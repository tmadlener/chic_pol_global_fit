#!/usr/bin/env python
"""
Script to make graphs as a function of ptm for differerent variables from a ptm
scan produced by run_ptM_scan
"""
import sys
import pandas as pd
import numpy as np

from collections import OrderedDict

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.data_handling import get_dataframe, apply_selections
from utils.misc_helpers import quantile
from utils.selection_functions import select_bin

from common_helpers import (
    get_var_name, get_variable_binning, frac_to_lam, identity
)

CHIC2_POL_SEL = select_bin('lth_chic2', -0.6, 1.0)
CHIC1_POL_SEL = select_bin('lth_chic1', -1./3., 1.0)


def process_input_variables(var_str):
    """Process the input variables into a usable format"""
    var_funcs = OrderedDict()
    for vstr in var_str.split(','):
        tmp = vstr.split(':')
        if len(tmp) == 2:
            func, var = tmp
            if not func in globals():
                logging.fatal('%s is not a valid (i.e. defined) transformation',
                              func)
                sys.exit(64)
            func = globals()[func]
        else:
            var = tmp[0]
            func = identity

        var_funcs[var] = func

    return var_funcs


def get_quantile_values(var_name, trans_f, quants=[0.16, 0.5, 0.84]):
    """Get the function that takes a dataframe as single argument and returns the
    desired quantiles of the variable after applying the passed
    transformation"""
    return lambda d: quantile(trans_f(d.loc[:, var_name].values), quants)



def make_graph(ptm_group_by, ptm_vals, var_name, trans_f):
    """Make a graph of the passed variable"""
    quantiles = ptm_group_by.apply(get_quantile_values(var_name, trans_f))

    central = np.ones_like(ptm_vals)
    low = np.zeros_like(ptm_vals)
    high = np.zeros_like(ptm_vals)

    for ipt in xrange(len(ptm_vals)):
        central[ipt] = quantiles.iloc[ipt][1] # quantiles is a pd.Series!
        low[ipt], high[ipt] = np.diff(quantiles.iloc[ipt])

    return r.TGraphAsymmErrors(len(ptm_vals), ptm_vals, central,
                               np.zeros_like(ptm_vals), np.zeros_like(ptm_vals),
                               low, high)



def main(args):
    """Main"""
    data = get_dataframe(args.scanfile)

    if args.physical_lambdas:
        logging.info('Removing rows with unphysical chic1 and chic2 lambdas')
        n_pre = data.shape[0]
        data = apply_selections(data, (CHIC1_POL_SEL, CHIC2_POL_SEL))
        logging.info('Rows before: {}, rows afterwards: {}'.format(n_pre, data.shape[0]))

    # Get a "helper" ptm binning, that is only used in the groupby to easily
    # identify the different ptM values
    ptm_binning, ptm_vals = get_variable_binning(data.ptM, get_vals=True)
    ptm_bins = data.groupby(pd.cut(data.ptM, ptm_binning))

    outfile = r.TFile.Open(args.outfile, 'recreate')

    for var, trans in args.variables.iteritems():
        graph = make_graph(ptm_bins, ptm_vals,  var, trans)
        graph.SetName('_'.join([get_var_name(var, trans.__name__), 'v', 'ptm']))
        graph.Write()


    outfile.Close()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make ptm dependent '
                                     'graphs for different variables from a ptm '
                                     'scan')
    parser.add_argument('scanfile', help='File containing the scan results from '
                        ' a run of run_ptM_scan')
    parser.add_argument('-o', '--outfile', help='File to which the TGraph '
                        'contour will be stored', default='ptm_dep_graphs.root')
    parser.add_argument('-v', '--variables', help='Comma separated list of '
                        'variables for which graphs should be produced. If a '
                        'transformation should be applied to the variable this '
                        'can be indicated by separating the transformation and '
                        'the variable with a colon: TRANSFORMATION:VARIABLE',
                        default='frac_to_lam:f_long_c1,frac_to_lam:f_long_c2',
                        type=process_input_variables)
    parser.add_argument('-pl', '--physical-lambdas', help='Clip the lambdas of '
                        'the chic1 and chic2 to their physically allowed ranges '
                        'before obtaining the bands', action='store_true',
                        default=False)

    clargs = parser.parse_args()
    main(clargs)
