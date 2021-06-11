#!/usr/bin/env python3
"""
Script to make a basic fit result summary pdf
"""

import os

from argparse import Namespace

from utils.reporting import open_tex_file, create_figure
from utils.data_handling import get_dataframe
from utils.misc_helpers import fmt_float

import make_fit_result_plots as plots
from make_fit_par_table import create_table, get_vals_uncers, print_val

PRETTY_PAR = {
    # LDMEs
    'l_3S1_8_c0': r'$\ldme{\chi_{c0}}{\Swave{3}{1}{8}}$',
    'l_3P0_1_c0': r'$\ldme{\chi_{c0}}{\Pwave{0}{1}}$',

    'l_3S1_1_jpsi': r'$\ldme{J/\psi}{\Swave{3}{1}{1}}$',
    'l_3S1_8_jpsi': r'$\ldme{J/\psi}{\Swave{3}{1}{8}}$', # derived
    'l_3PJ_8_jpsi': r'$\ldme{J/\psi}{\Pwave{J}{8}}$', # derived

    'l_3S1_1_psip': r'$\ldme{\psi(2S)}{\Swave{3}{1}{1}}$',
    'l_3S1_8_psip': r'$\ldme{\psi(2S)}{\Swave{3}{1}{8}}$', # derived
    'l_3PJ_8_psip': r'$\ldme{\psi(2S)}{\Pwave{J}{8}}$', # derived

    'l_1S0_8_jpsi': r'$\ldme{J/\psi}{\Swave{1}{0}{8}}$',
    'l_1S0_8_psip': r'$\ldme{\psi(2S)}{\Swave{1}{0}{8}}$',

    # LDME ratios
    'l_r_3PJ_8_1S0_8_jpsi': r'$\ldmeratio{\Pwave{J}{8}}{J/\psi}$',
    'l_r_3S1_8_1S0_8_jpsi': r'$\ldmeratio{\Swave{3}{1}{8}}{J/\psi}$',

    # double ratios
    'l_rr_3PJ_8_1S0_8_psip_jpsi': r'$\ldmedr{\Pwave{J}{8}}$',
    'l_rr_3S1_8_1S0_8_psip_jpsi': r'$\ldmedr{\Swave{3}{1}{8}}$',

    # costh ratio norms
    'norm_costh_1': r'$n_{1}$',
    'norm_costh_2': r'$n_{2}$',
    'norm_costh_3': r'$n_{3}$',

    # Nuisance params
    'br_psip_dp': r'\br{\psi(2S)}{J/\psi \pi\pi}',
    'br_psip_mm': r'\br{\psi(2S)}{\mu\mu}',
    'br_psip_c2': r'\br{\psi(2S)}{\chi_{c2}}',
    'br_psip_c1': r'\br{\psi(2S)}{\chi_{c1}}',
    'br_psip_jpsi': r'\br{\psi(2S)}{J/\psi}',
    'br_c2_jpsi': r'\br{\chi_{c2}}{J/\psi}',
    'br_c1_jpsi': r'\br{\chi_{c1}}{J/\psi}',
    'br_jpsi_mm': r'\br{J/\psi}{\mu\mu}',
    'L_ATLAS': r'\lumi{ATLAS}',
    'L_CMS': r'\lumi{ATLAS}'
}

DERIV_PARS = (
    'l_3S1_8_jpsi',
    'l_3PJ_8_jpsi',
    'l_3S1_8_psip',
    'l_3PJ_8_psip',
)

PREAMBLE = r'''\documentclass[a4paper, 11pt]{scrartcl}

\usepackage{graphicx}
\usepackage{subfig}
\usepackage[margin=2cm]{geometry}\usepackage{tabulary}
\usepackage{multirow}
\usepackage{amsmath}

\newcommand{\lumi}[1]{$\mathcal{L}_{\textrm{#1}}$}
\newcommand{\br}[2]{$\mathcal{B}(#1 \to #2)$}
\newcommand{\Swave}[3]{{}^{#1}S_{#2}^{[#3]}}
\newcommand{\Pwave}[2]{{}^{3}P_{#1}^{[#2]}}
\newcommand{\ldme}[2]{\mathcal{O}^{#1}(#2)}
\newcommand{\ldmeratio}[2]{\mathcal{R}_{#2}(#1, \Swave{1}{0}{8})}
\newcommand{\ldmedr}[1]{\mathcal{RR}(#1, \Swave{1}{0}{8})}
'''

SYMBOL_CAPTION = r'''
$\mathcal{R}_{\mathcal{Q}}(X, Y) \equiv \ldme{\mathcal{Q}}{X} / \ldme{\mathcal{Q}}{Y}$.
$\mathcal{RR}(X, Y) \equiv \mathcal{R}_{\psi(2S)}(X, Y) / \mathcal{R}_{J/\psi}(X, Y)$
'''

PLOTS_TO_USE = (
    'combined_cs.pdf',
    'psi_pol.pdf',
    'chic_ratio_cs.pdf',
    'costh_ratio_ptm_3p29.pdf',
    'costh_ratio_ptm_4p64.pdf',
    'costh_ratio_ptm_7p10.pdf',
)

PLOT_DIR = 'plots_for_report'

def add_fit_plots(plot_dir):
    """Add the plots from the plot_dir relative to the result dir"""
    def _label(pdf_name):
        return pdf_name.replace('.pdf', '').replace('_', r'\_')

    return create_figure({_label(p): f'./{plot_dir}/{p}' for p in PLOTS_TO_USE})

def pretty_label(par):
    """Get a pretty label if one exists"""
    if par in PRETTY_PAR: return PRETTY_PAR[par]
    return par.replace('_', r'\_')

def add_fit_param_table(res_file):
    """Add the fit parameter table from the fit result file (relative to the result dir)"""
    par_indices = get_dataframe(res_file, 'parameter_indices')
    par_indices = {
        p: i for p, i in zip(par_indices.columns, par_indices.values[0])
    }

    vals_uncers = get_vals_uncers(res_file, par_indices)
    rows = [r'parameter & best-fit value \\ \hline']
    for par, (val, err) in vals_uncers.items():
        rows.append(
            fr'{pretty_label(par)} & ${fmt_float(val)} \pm {fmt_float(err)}$ \\'
        )

    return create_table(rows, 'l | c', caption='Fit parameter values. ' + SYMBOL_CAPTION)

def add_deriv_par_table(data):
    """Add the table with the derived parameters"""
    rows = [r'parameter & value \\ \hline']

    for par in DERIV_PARS:
        rows.append(fr'{pretty_label(par)} & {print_val(data, par)} \\')

    return create_table(rows, 'l | c', caption='Derived parameter values at $p_{T}/M = 5$')

def main(args):
    """Main"""
    graphfile = os.path.join(args.resultdir, 'fit_graphs_and_models_nrqcd_global_fit.root')

    # Produce the plots
    plot_args = Namespace(graphfile=graphfile, outdir=f'{args.resultdir}/{PLOT_DIR}')
    plots.main(plot_args)

    # Scan data loading
    scan_data = get_dataframe(f'{args.resultdir}/param_scan_ptm5.root')

    with open_tex_file(f'{args.resultdir}/fit_results_report.tex', PREAMBLE) as texfile:
        texfile.write(add_fit_plots(PLOT_DIR))
        texfile.write('\n')

        texfile.write(add_fit_param_table(f'{args.resultdir}/fit_results_nrqcd_global_fit.root'))
        texfile.write('\n')

        texfile.write(add_deriv_par_table(scan_data))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to generate an overview fit result pdf')
    parser.add_argument('resultdir', help='Directory containing the fit results')

    clargs = parser.parse_args()
    main(clargs)
