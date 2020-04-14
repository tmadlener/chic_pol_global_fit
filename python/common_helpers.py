#!/usr/bin/env python
"""
Module containing some common helper functionality
"""

import numpy as np


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
    else:
        return binning

