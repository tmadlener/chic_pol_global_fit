#!/usr/bin/env python
"""
Module containing some common helper functionality
"""

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
