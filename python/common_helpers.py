#!/usr/bin/env python
"""
Module containing some common helper functionality
"""
import numpy as np

from scipy.spatial import ConvexHull
from scipy.optimize import root_scalar

import ROOT as r

from utils.hist_utils import hist2d, get_array, get_binning
from utils.misc_helpers import quantile

def frac_to_lam(frac):
    """Transform a longitudinal fraction into a lambda value"""
    return (1 - 3 * frac) / (1 + frac)


def identity(x):
    """Identity function"""
    return x


def get_var_name(var, trans):
    """Get the variable name for storing the contour"""
    if trans == 'identity':
        return var
    return var.replace('f_long_c', 'lambda_')


def get_variable_binning(var_vals, get_vals=False):
    """Get a binning that places every data point into a separate bin. If get_vals
    is True, also return the unique values"""
    # Get the unique values, and split the intervals between them by half to
    # get the binning. For the first and last bin reuse the intervals of the
    # adjacent bins
    vals = np.sort(np.unique(var_vals))
    diffs = np.diff(vals)

    diffs = np.append(np.array([diffs[0]]), diffs)
    # Subtract half the interval to get to the bin borders from the central
    # values
    binning = vals - diffs / 2
    # Now have to use the full interval for the last bin
    binning = np.append(binning, np.array([binning[-1]] + diffs[-1]))

    if get_vals:
        return binning, vals

    return binning


def contour_as_tgraph(filled_ps):
    """Get the convex hull of the filled points and return a (closed TGraph)"""
    if not len(filled_ps):
        print('No filled points')
        return r.TGraph()

    hull = ConvexHull(filled_ps)

    # "close" the contour by appending the first point again at the end
    xcont = filled_ps[hull.vertices, 0]
    xcont = np.append(xcont, np.array(xcont[0]))
    ycont = filled_ps[hull.vertices, 1]
    ycont = np.append(ycont, np.array(ycont[0]))

    return r.TGraph(len(hull.vertices) + 1, xcont, ycont)


def get_2d_hist(data, varx, vary, nbinsx=100, nbinsy=100):
    """Get the 2d distribution of vary vs varx"""
    # Get rid of very large outliers
    xmin, xmax = quantile(data.loc[:, varx].values, [0.0001, 0.9999])
    ymin, ymax = quantile(data.loc[:, vary].values, [0.0001, 0.9999])
    return hist2d(data.loc[:, varx], data.loc[:, vary], minx=xmin, maxx=xmax,
                  miny=ymin, maxy=ymax, nbinsx=nbinsx, nbinsy=nbinsy)


def get_coverage_contour(hist, coverage=0.683):
    """
    Get the contour from the passed histogram that surpasses the specified coerage
    """
    vals = get_array(hist)
    sum_vals = np.sum(vals)

    def _coverage(level):
        """Calculate the coverage corresponding to the passed level"""
        return np.sum(vals * (vals >= level)) / sum_vals

    cov_level = root_scalar(lambda x: _coverage(x) - coverage,
                            bracket=[np.min(vals), np.max(vals)])

    filled = vals >= cov_level.root

    x_vals, y_vals = get_binning(hist, 'X'), get_binning(hist, 'Y')
    # get the bin centers
    x_vals = 0.5 * (x_vals[1:] + x_vals[:-1])
    y_vals = 0.5 * (y_vals[1:] + y_vals[:-1])

    filled_coords = []
    for ix, xv in enumerate(x_vals):
        for iy, yv in enumerate(y_vals):
            if filled[ix, iy]:
                filled_coords.append([xv, yv])

    return contour_as_tgraph(np.array(filled_coords))
