#!/usr/bin/env python
"""
Script to obtain the 2d contour from a fit scan
"""
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.stats import chi2

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')


from utils.data_handling import get_dataframe, apply_selections
from utils.hist_utils import get_array, get_binning
from utils.selection_functions import select_bin

from common_helpers import (
    frac_to_lam, identity, get_var_name, contour_as_tgraph
)

def has_branch(tree, branches):
    """
    Check whether the tree has all the branches

    Args:
        tree (ROOT.TTree): TTree which should be checked
        branches (or list of str): Names of the branches for which the existence
            in the tree should be checked

    Returns:
        boolean: True if all branches exist in the TTree, False if one or more
            are not existing
    """
    tree_branches = [b.GetName() for b in tree.GetListOfBranches()]
    return all(b in tree_branches for b in branches)


def load_data(scanfile, var_x, var_y, tree='log_like_scan'):
    """Load the dataframe containing the scan results and do some cleanup"""
    # Only load the necessary variables
    load_vars = [var_x, var_y, 'llh']
    rfile = r.TFile.Open(scanfile)
    if has_branch(rfile.Get(tree), ['goodFit']):
        load_vars.append('goodFit')
    rfile.Close()

    data = get_dataframe(scanfile, treename=tree, columns=load_vars)

    n_full = data.shape[0]

    logging.info('Loaded %d rows of data', n_full)

    # Check if the data come from a fit scan
    if 'goodFit' in data.columns:
        # remove all rows where the fit did not converge and where there are
        # nan values for any of the fitted parameters
        data = apply_selections(data, lambda d: d.goodFit > 0)
        n_bad_fit = n_full - data.shape[0]

        logging.info('Removed %d because of no good fit', n_bad_fit)

    return data


def get_contour_graph(data, conf_level, var_x, var_y,
                      tf_x='identity', tf_y='identity',
                      bound_x=None, bound_y=None):
    """Get the two dimensional histogram from which the contour will be obtained"""
    trans_f_x = globals()[tf_x]
    trans_f_y = globals()[tf_y]

    llh_min = data.llh.min()
    llh_cond = lambda d: 2 * (d.llh - llh_min) < chi2.ppf(conf_level, 2)
    # remove the minimum value data point. In case it is outside the scanned
    # range, than this would "distort" the contour.
    rm_min = lambda d: d.llh != llh_min

    select_funcs = [llh_cond, rm_min]
    if bound_x is not None:
        select_funcs.append(select_bin(var_x, *bound_x))
    if bound_y is not None:
        select_funcs.append(select_bin(var_y, *bound_y))

    sel_data = apply_selections(data, select_funcs)

    filled = np.array(sel_data.loc[:, [var_x, var_y]])
    filled[:, 0] = trans_f_x(filled[:, 0])
    filled[:, 1] = trans_f_y(filled[:, 1])

    return contour_as_tgraph(filled)


def get_best_fit_graph(data, var_x, var_y,
                       tf_x='identity', tf_y='identity'):
    """Get the two-dimensional best fit point"""
    trans_f_x = globals()[tf_x]
    trans_f_y = globals()[tf_y]

    llh_min = data.llh.min()

    # In case the scan file has been merged from more than one scans, the
    # minimum will be present more than once
    min_data = data[data.llh == llh_min]
    if data.shape[1] > 0:
        min_data = min_data.iloc[0]

    return r.TGraph(1, np.array(trans_f_x(min_data.loc[var_x])),
                    np.array(trans_f_y(min_data.loc[var_y])))


def main(args):
    """Main"""
    llh_data = load_data(args.scanfile, args.variable_x, args.variable_y)

    outfile = r.TFile(args.outfile, 'recreate')
    outfile.cd()
    name_x = get_var_name(args.variable_x, args.transform_x)
    name_y = get_var_name(args.variable_y, args.transform_y)

    for conf_l in args.conf_levels:
        contour = get_contour_graph(llh_data, conf_l,
                                    args.variable_x, args.variable_y,
                                    args.transform_x, args.transform_y)


        contour.SetName('_'.join(['contour', name_x, name_y,
                                  str(conf_l).replace('.', 'p')]))
        contour.Write()


    best_fit = get_best_fit_graph(llh_data, args.variable_x, args.variable_y,
                                  args.transform_x, args.transform_y)
    best_fit.SetName('_'.join(['best', 'fit', name_x, name_y]))
    best_fit.Write()

    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Obtain and store the 2d '
                                     'contour using a fit scan as input')
    parser.add_argument('scanfile', help='File containing the TTree with the '
                        'results from the scan')
    parser.add_argument('-o', '--outfile', help='File to which the TGraph '
                        'contour will be stored', default='contour_graph.root')
    parser.add_argument('-vx', '--variable-x', help='Variable to use for the '
                        'horizontal direction', default='lambda_1')
    parser.add_argument('-vy', '--variable-y', help='Variable to use for the '
                        'vertical direction', default='lambda_2')
    parser.add_argument('-tx', '--transform-x', help='Transformation to apply to'
                        ' the variable on the horizontal axis (identity or '
                        'frac_to_lam)', default='identity')
    parser.add_argument('-ty', '--transform-y', help='Transformation to apply to'
                        ' the variable on the vertical axis (identity or '
                        'frac_to_lam)', default='identity')
    parser.add_argument('-cl', '--conf-levels', help='The CL levels (0, 1) for '
                        'which the contour should be obtained', default='0.68',
                        type=lambda s: [float(v) for v in s.split(',')])


    clargs = parser.parse_args()
    main(clargs)
