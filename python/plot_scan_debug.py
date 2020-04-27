#!/usr/bin/env python
"""
Script to make some debug plots from the scanning
"""
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.data_handling import get_dataframe, apply_selections
from utils.hist_utils import (
    get_array, get_binning, from_array, hist2d, find_bin
)
from utils.misc_helpers import cond_mkdir
from utils.setup_plot_style import set_basic_style
from utils.plot_helpers import mkplot, _setup_canvas

from common_helpers import (
    frac_to_lam, identity, get_var_name, get_variable_binning
)

def to_bw_hist(hist):
    """Fill all filled bins with value 1 and all empty ones with 0"""
    arr = get_array(hist)
    # TODO: generalize and put into hist_utils
    binning = np.array([get_binning(hist, 0), get_binning(hist, 1)])

    arr = arr != 0

    return from_array(arr, binning)


def remove_color_bar(can, hist_idx=1):
    """Remove the color bar from the histogram"""
    hist = can.pltables[hist_idx]
    palette = hist.GetListOfFunctions().FindObject('palette')
    palette.SetX1NDC(1.2)
    palette.SetX2NDC(1.3)
    can.Modified()
    can.Update()


def make_plot_good_fit(data, var_x, var_y, tf_x, tf_y):
    """Make a plot showing whether the fit was successful or not for each point
    of the scan"""
    # Remove the minimum since that is in any case valid
    data = apply_selections(data, lambda d: d.llh != d.llh.min())

    trans_f_x = globals()[tf_x]
    trans_f_y = globals()[tf_y]

    x_vals = trans_f_x(data.loc[:, var_x])
    y_vals = trans_f_y(data.loc[:, var_y])

    x_binning = get_variable_binning(x_vals)
    y_binning = get_variable_binning(y_vals)

    good_fit = data.goodFit > 0

    hist = hist2d(x_vals[good_fit], y_vals[good_fit],
                  x_hist_sett=(len(x_binning) - 1, x_binning),
                  y_hist_sett=(len(y_binning) - 1, y_binning))

    can = mkplot(to_bw_hist(hist), drawOpt='colz',
                 xRange=[x_binning[0], x_binning[-1]],
                 yRange=[y_binning[0], y_binning[-1]],
                 xLabel=get_var_name(var_x, tf_x),
                 yLabel=get_var_name(var_y, tf_y))

    hist = can.pltables[1]
    hist.GetZaxis().SetRangeUser(0, 10)

    remove_color_bar(can)

    return can


def make_plot_min_chi2(data, var_x, var_y, tf_x, tf_y, gf_only=False):
    """Make a plot showing the min chi2 from the scan for each scan point, where
    the fit was OK"""
    data.loc[:, 'delta_chi2'] = 2 * (data.llh - data.llh.min())

    bin_data = apply_selections(data, lambda d: d.delta_chi2 != 0)
    # Cleanup the data frame by dropping duplicates which can stem from
    # overlapping scanning when done in parallel
    bin_data = bin_data.drop_duplicates()

    if gf_only:
        bin_data = apply_selections(bin_data, lambda d: d.goodFit > 0)

    # Get the binning
    trans_f_x = globals()[tf_x]
    trans_f_y = globals()[tf_y]
    x_vals = trans_f_x(bin_data.loc[:, var_x])
    y_vals = trans_f_y(bin_data.loc[:, var_y])

    x_binning = get_variable_binning(x_vals)
    y_binning = get_variable_binning(y_vals)

    arr = np.zeros((len(x_binning) - 1, len(y_binning) - 1))
    x_bin = find_bin(x_binning, x_vals)
    y_bin = find_bin(y_binning, y_vals)

    dchi2 = bin_data.delta_chi2.values
    arr[x_bin, y_bin] = dchi2[:]

    hist = from_array(arr, np.array([x_binning, y_binning]))

    can = _setup_canvas(None)
    can.SetRightMargin(0.12)

    mkplot(hist, can=can, drawOpt='colz',
           xRange=[x_binning[0], x_binning[-1]],
           yRange=[y_binning[0], y_binning[-1]],
           xLabel=get_var_name(var_x, tf_x),
           yLabel=get_var_name(var_y, tf_y))

    hist.GetZaxis().SetRangeUser(0, 25)

    can.Update()


    return can


def main(args):
    """Main"""
    set_basic_style()
    llh_data = get_dataframe(args.scanfile, 'log_like_scan')

    cond_mkdir(args.outdir)

    can = make_plot_good_fit(llh_data, args.variable_x, args.variable_y,
                             args.transform_x, args.transform_y)
    can.SaveAs('/'.join([args.outdir, 'scan_good_fit.pdf']))


    can = make_plot_min_chi2(llh_data, args.variable_x, args.variable_y,
                             args.transform_x, args.transform_y)
    can.SaveAs('/'.join([args.outdir, 'chi2_vals_all_fits.pdf']))

    can = make_plot_min_chi2(llh_data, args.variable_x, args.variable_y,
                             args.transform_x, args.transform_y,
                             gf_only=True)
    can.SaveAs('/'.join([args.outdir, 'chi2_vals_good_fits.pdf']))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make scan debug'
                                     ' plots')
    parser.add_argument('scanfile', help='File containing the results of the '
                        'scan')
    parser.add_argument('-o', '--outdir', help='Output directory into which all '
                        'the plots will be stored', default='.')
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



    clargs = parser.parse_args()
    main(clargs)
