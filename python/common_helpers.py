#!/usr/bin/env python
"""
Module containing some common helper functionality
"""
import numpy as np
from scipy.spatial import ConvexHull
from collections import OrderedDict
from functools import partial

import ROOT as r

from utils.hist_utils import hist2d, get_array, get_binning, find_bin
from utils.misc_helpers import quantile
from utils.data_handling import apply_selections
from utils.selection_functions import select_bin

class RootResult(object):
    """Helper class to mimic the result from root_scalar"""
    def __init__(self, val, success=True):
        self.root = val
        self.converged = success


def secant(func, bracket, eps=1e-4, maxsteps=100, nsteps=0):
    """Simple root finding for cases where scipy is not recent enough and
    root_scalar is not present
    """
    x0, x1 = bracket
    fx0 = func(x0)
    fx1 = func(x1)
    x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0)

    if np.abs(func(x2)) < eps:
        return RootResult(x2)

    if nsteps >= maxsteps:
        print('Could not find root for function within {} steps. Value = {}'.format(maxsteps, func(x2)))
        return RootResult(x2, False)

    return secant(func, [x0, x1], eps, maxsteps, nsteps + 1)

try:
    from scipy.optimize import root_scalar
except ImportError: # scipy not recent enough
    root_scalar = partial(secant, eps=1e-3)


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

    # do some pre-processing to start from a slightly better bracket for the
    # secant method
    dec_cov = np.array([_coverage(0.05 * i * np.max(vals)) for i in xrange(21)])
    q_bin = find_bin(dec_cov, np.array([coverage]))
    search_brack = [q_bin * 0.05 * np.max(vals), (q_bin + 1) * 0.05 * np.max(vals)]

    cov_level = root_scalar(lambda x: _coverage(x) - coverage,
                            bracket=search_brack)

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


LABELS = {
    'lth': '#lambda_{#vartheta}',
    'dlth': '#Delta#lambda',
    'jpsi': 'J/#psi',
    'psip': '#psi(2S)',
    'chic1': '#chi_{c1}',
    'chic2': '#chi_{c2}',
    'chic': '#chi_{c}'
}


def get_label(var):
    """Get a label for the passed variable"""
    parts = var.split('_')
    if parts[0] == 'r':
        return '{} #rightarrow {} feed down'.format(LABELS.get(parts[1], parts[1]),
                                                    LABELS.get(parts[2], parts[2]))

    if parts[0] == 'lth':
        if len(parts) == 2:
            return '{}^{{{}}}'.format(LABELS['lth'], LABELS.get(parts[1], parts[1]))
        return '{}^{{{}}} from {}'.format(LABELS['lth'], LABELS.get(parts[1], parts[1]),
                                          LABELS.get(parts[2], parts[2]))

    if parts[0] == 'dlth':
        return '{}({},{})'.format(LABELS['dlth'], LABELS.get(parts[1], parts[1]),
                                  LABELS.get(parts[2], parts[2]))

    return var


def get_range(var):
    """Get the plotting range for a variable"""
    parts = var.split('_')
    if parts[0] == 'r':
        if parts[1] == 'chic2':
            return [0, 0.15]
        return [0, 0.3]

    if parts[0] == 'lth':
        if 'psi' in parts[1]:
            return [-0.75, 0.75]
        return [-1.5, 2.0]

    if parts[0] == 'dlth':
        return [-0.75, 0.75]

    return [None, None]


def parse_inputs(inputs, enforce_keys=False):
    files = OrderedDict()
    for inp in inputs:
        try:
            key, fn = inp.split(':')
        except ValueError as e:
            # When there is only one input make it possible to pass that
            # without a legend identifier unless they keys are enforced
            if enforce_keys:
                raise(e)
            key, fn = 'no_leg', inp

        files[key] = r.TFile.Open(fn)

    return files



LTH_CHI1_SEL = select_bin('lth_chic1', -1./3., 1)
LTH_CHI2_SEL = select_bin('lth_chic2', -0.6, 1)

def select_phys_data(dfr, plot_vars, do_selection=True):
    """Select the data to plot according to the variables that are plotted. If
    lth_chicJ is in the plot_vars no selection will be applied to it, if it is
    projected over, it will be restricted to the physically allowed domain"""
    if not do_selection:
        return dfr

    if not isinstance(plot_vars, (list, tuple)):
        plot_vars = [plot_vars]

    selections = []
    if 'lth_chic1' not in plot_vars:
        selections.append(LTH_CHI1_SEL)
    if 'lth_chic2' not in plot_vars:
        selections.append(LTH_CHI2_SEL)

    selections = selections if selections else None

    return apply_selections(dfr, selections)
