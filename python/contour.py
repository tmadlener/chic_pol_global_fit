#!/usr/bin/env python
"""
Script to obtain the 2d contour from a fit scan
"""
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.spatial import ConvexHull
from scipy.stats import chi2

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')


from utils.data_handling import get_dataframe, apply_selections
from utils.hist_utils import get_array, get_binning
from utils.selection_functions import select_bin


def frac_to_lam(frac):
    """Transform a longitudinal fraction into a lambda value"""
    return (1 - 3 * frac) / (1 + frac)


def identity(x):
    """Identity function"""
    return x


def load_data(scanfile):
    """Load the dataframe containing the scan results and do some cleanup"""
    data = get_dataframe(scanfile)

    n_full = data.shape[0]

    # remove all rows where the fit did not converge and where there are nan
    # values for any of the fitted parameters
    data = apply_selections(data, lambda d: d.goodFit == 1)
    n_bad_fit = n_full - data.shape[0]

    logging.info('Loaded %d rows and removed %d because of no good',
                 n_full, n_bad_fit)

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

    hull = ConvexHull(filled)

    # "close" the contour by appending the first point again at the end
    xcont = filled[hull.vertices, 0]
    xcont = np.append(xcont, np.array(xcont[0]))
    ycont = filled[hull.vertices, 1]
    ycont = np.append(ycont, np.array(ycont[0]))

    return r.TGraph(len(hull.vertices) + 1, xcont, ycont)


def get_var_name(var, trans):
    """Get the variable name for storing the contour"""
    if trans == 'identity':
        return var
    return var.replace('f_long', 'lambda')


def main(args):
    """Main"""
    llh_data = load_data(args.scanfile)

    contour = get_contour_graph(llh_data, args.conf_level,
                                args.variable_x, args.variable_y,
                                args.transform_x, args.transform_y)

    outfile = r.TFile(args.outfile, 'recreate')
    outfile.cd()
    name_x = get_var_name(args.variable_x, args.transform_x)
    name_y = get_var_name(args.variable_y, args.transform_y)

    contour.SetName('_'.join(['contour', name_x, name_y,
                              str(args.conf_level).replace('.', '0')]))
    contour.Write()
    outfile.Close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Obtain and store the 2d '
                                     'contour using a fit scan as input')
    parser.add_argument('scanfile', help='File containing the TTree with the '
                        'results from the scan')
    parser.add_argument('-o', '--outfile', help='File to which the TGraph '
                        'contour will be stored', default='contour_graph.root')
    parser.add_argument('-vx', '--variable_x', help='Variable to use for the '
                        'horizontal direction', default='lambda_1')
    parser.add_argument('-vy', '--variable_y', help='Variable to use for the '
                        'vertical direction', default='lambda_2')
    parser.add_argument('-tx', '--transform_x', help='Transformation to apply to'
                        ' the variable on the horizontal axis (identity or '
                        'frac_to_lam)', default='identity')
    parser.add_argument('-ty', '--transform_y', help='Transformation to apply to'
                        ' the variable on the vertical axis (identity or '
                        'frac_to_lam)', default='identity')
    parser.add_argument('-cl', '--conf-level', help='The CL level (0, 1) for '
                        'which the contour should be obtained', default=0.68,
                        type=float)


    clargs = parser.parse_args()
    main(clargs)
