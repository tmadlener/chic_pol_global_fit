#!/usr/bin/env python
"""
Script to make a basic fit parameter table
"""
import numpy as np

from collections import OrderedDict
from copy import deepcopy

from utils.data_handling import get_dataframe
from utils.reporting import open_tex_file
from utils.misc_helpers import quantile

from common_helpers import select_phys_data

PAR_TO_IDX = {
    "sigma_psip": 0,
    "sigma_chic2": 1,
    "sigma_chic1": 2,
    "sigma_jpsi": 3,

    "f_long_psi": 4,
    "f_long_c1": 5,
    "f_long_c2": 6,

    "gamma": 7,
    "beta_total_psi": 8,
    "beta_long_psi": 9,
    "beta_total_c1": 10,
    "beta_long_c1": 11,
    "beta_total_c2": 12,
    "beta_long_c2": 13,

    "br_psip_dp": 14,
    "br_psip_mm": 15,
    "br_psip_c2": 16,
    "br_psip_c1": 17,
    "br_psip_jpsi": 18,
    "br_c2_jpsi": 19,
    "br_c1_jpsi": 20,
    "br_jpsi_mm": 21,

    "L_CMS": 22,
    "L_ATLAS": 23,

    "norm_costh_1": 24,
    "norm_costh_2": 25,
    "norm_costh_3": 26,
}

# Parameters that appear after the betas in the indices
REST_PARS = [
    "br_psip_dp", "br_psip_mm", "br_psip_c2", "br_psip_c1",
    "br_psip_jpsi", "br_c2_jpsi", "br_c1_jpsi", "br_jpsi_mm",
    "L_CMS", "L_ATLAS",
    "norm_costh_1", "norm_costh_2", "norm_costh_3"
]

PRETTY_PAR = OrderedDict()

PRETTY_PAR['sigma_psip'] = r'$\sigma^{\psi(2S)}$'
PRETTY_PAR['sigma_chic2'] = r'$\sigma^{\chi_{c2}}$'
PRETTY_PAR['sigma_chic1'] = r'$\sigma^{\chi_{c1}}$'
PRETTY_PAR['sigma_jpsi'] = r'$\sigma^{J/\psi}$'

PRETTY_PAR['f_long_psi'] = r'$f_{l}^{\psi}$'
PRETTY_PAR['f_long_c1'] = r'$f_{l}^{\chi_{c1}}$'
PRETTY_PAR['f_long_c2'] = r'$f_{l}^{\chi_{c2}}$'

PRETTY_PAR['gamma'] = r'$\gamma$'
PRETTY_PAR['beta_total_psi'] = r'$\beta_{t}^{\psi}$'
PRETTY_PAR['beta_total_c1'] = r'$\beta_{t}^{\chi_{c1}}$'
PRETTY_PAR['beta_total_c2'] = r'$\beta_{t}^{\chi_{c2}}$'
PRETTY_PAR['beta_long_psi'] = r'$\beta_{l}^{\psi}$'
PRETTY_PAR['beta_long_c1'] = r'$\beta_{l}^{\chi_{c1}}$'
PRETTY_PAR['beta_long_c2'] = r'$\beta_{l}^{\chi_{c2}}$'

PRETTY_PAR['norm_costh_1'] = r'$n_{1}$'
PRETTY_PAR['norm_costh_2'] = r'$n_{2}$'
PRETTY_PAR['norm_costh_3'] = r'$n_{3}$'

PRETTY_PAR['br_psip_dp'] = r'\br{\psi(2S)}{J/\psi \pi\pi}'
PRETTY_PAR['br_psip_mm'] = r'\br{\psi(2S)}{\mu\mu}'
PRETTY_PAR['br_psip_c2'] = r'\br{\psi(2S)}{\chi_{c2}}'
PRETTY_PAR['br_psip_c1'] = r'\br{\psi(2S)}{\chi_{c1}}'
PRETTY_PAR['br_psip_jpsi'] = r'\br{\psi(2S)}{J/\psi}'
PRETTY_PAR['br_c2_jpsi'] = r'\br{\chi_{c2}}{J/\psi}'
PRETTY_PAR['br_c1_jpsi'] = r'\br{\chi_{c1}}{J/\psi}'
PRETTY_PAR['br_jpsi_mm'] = r'\br{J/\psi}{\mu\mu}'
PRETTY_PAR['L_ATLAS'] = r'\lumi{ATLAS}'
PRETTY_PAR['L_CMS'] = r'\lumi{CMS}'


DERIV_PAR = {
    'chic2_chic1_cs_br': r'$\mathcal{B}\sigma_{\chi_{c2}} / \mathcal{B}\sigma_{\chi_{c1}}$',
    'r_chic1_jpsi': r'$\chi_{c1} \to J/\psi$',
    'r_chic2_jpsi': r'$\chi_{c2} \to J/\psi$',
    'r_chic_jpsi': r'$\chi_{cJ} (J=1,2) \to J/\psi$',
    'r_psip_jpsi': r'$\psi(2S) \to J/\psi$',
    'r_jpsi_direct': r'direct $J/\psi$',
    'r_psip_chic1': r'$\psi(2S) \to \chi_{c1}$',
    'r_psip_chic2': r'$\psi(2S) \to \chi_{c2}$',
    'lth_jpsi': r'$J/\psi$',
    'lth_chic1': r'$\chi_{c1}$',
    'lth_chic2': r'$\chi_{c2}$',
    'lth_psip': r'$\psi(2S)$',
    'lth_jpsi_chic': r'$J/\psi$ from $\chi_{cJ}$ (J=1,2)',
}



PREAMBLE = r'''\documentclass[a4paper, 11pt]{scrartcl}

\usepackage{tabulary}
\usepackage{multirow}
\usepackage{amsmath}

\newcommand{\lumi}[1]{$\mathcal{L}_{\textrm{#1}}$}
\newcommand{\br}[2]{$\mathcal{B}(#1 \to #2)$}

'''

def get_par_to_idx(fit_option=1):
    """
    Get the table of parameter names to indices in the parameter value vector
    """
    if fit_option == 0:
        return PAR_TO_IDX

    table = deepcopy(PAR_TO_IDX)

    if fit_option == 1:
        table['beta_long_c1'] = 10
        table['beta_total_c1'] = 8
        table['beta_long_c2'] = 11
        table['beta_total_c2'] = 8

        for par in REST_PARS:
            table[par] -= 2

        return table

    if fit_option == 2:
        table['beta_long_psi'] = 8
        table['beta_total_c1'] = 8
        table['beta_total_c2'] = 8
        table['beta_long_c1'] = 8
        table['beta_long_c2'] = 8

        for par in REST_PARS:
            table[par] -= 5

        return table


def get_vals_uncers(datafile, par_to_idx):
    """Get the values and uncertainties of all parameters"""
    data = get_dataframe(datafile, 'fit_result')
    pars = data.parameters[0]
    n_pars = np.max([int(i) for i in par_to_idx.values()]) + 1
    if len(pars) != n_pars:
        print('Parameters and par_to_idx do not have the same length!')
        return

    cov = np.reshape(data.cov_matrix[0], (n_pars, n_pars))

    vals = {}
    for par, idx in par_to_idx.items():
        vals[par] = (pars[idx], np.sqrt(cov[idx, idx]))

    return vals


def create_table(rows, table_format, caption=None, label=None):
    """Create a table from the passed rows"""
    table = []
    table.append(r'\begin{table}')
    table.append(r'\centering')
    table.append(r'\begin{{tabulary}}{{1.0\linewidth}}{{{}}}'.format(table_format))

    table += rows

    table.append(r'\end{tabulary}')

    if caption is not None:
        table.append(r'\caption{{{}}}'.format(caption))
    if label is not None:
        table.append(r'\label{{{}}}'.format(label))

    table.append(r'\end{table}')
    return '\n'.join(table)


def create_fit_par_table(v_free_b, v_e_t_cs, v_equal_b):
    """Get the table string"""
    rows = [r' & free & default & constrained \\ \hline']
    fmt_vals = lambda v: r'${:.3f} \pm {:.3f}$'.format(*v)

    for par, ltx in PRETTY_PAR.iteritems():
        rows.append(
            r'{} & {} & {} & {} \\'.format(ltx, fmt_vals(v_free_b[par]),
                                           fmt_vals(v_e_t_cs[par]),
                                           fmt_vals(v_equal_b[par])))

    return create_table(rows, 'l|c|c|c', caption='Best fit parameters and uncertainties')


def print_val(data, var, percent=False):
    data = select_phys_data(data, var)
    quants = quantile(data.loc[:, var].values, [0.16, 0.5, 0.84])
    err_lo, err_hi = np.diff(quants)

    asym = np.abs(err_lo - err_hi) > 1e-3

    if percent:
        fmt_str_asym = r'${:.1f}^{{+{:.1f}}}_{{-{:.1f}}}$'
        fmt_str_sym = r'${:.1f} \pm {:.1f}$'
        quants *= 100
        err_lo, err_hi = np.diff(quants)
    else:
        fmt_str_asym = r'${:.3e}^{{+{:.3e}}}_{{-{:.3e}}}$'
        fmt_str_sym = r'${:.3e} \pm {:.3e}$'

    if asym:
        return fmt_str_asym.format(quants[1], err_hi, err_lo)
    else:
        return fmt_str_sym.format(quants[1], err_lo)


def create_deriv_par_table(data_fb, data_def, data_eb):
    """
    Make the tables with the derived parameters
    """
    rows = [r' & $6\beta$ & $4\beta$ & $1\beta$ \\ \hline']

    fmt_row = lambda n, p0, p1, p2: r'{} & {} & {} & {} \\'.format(n, p0, p1, p2)

    par = 'chic2_chic1_cs_br'
    rows.append(
        fmt_row(DERIV_PAR[par], print_val(data_fb, par),
                print_val(data_def, par), print_val(data_eb, par))
    )
    rows.append(r'\hline\hline')
    rows.append(r'\multicolumn{4}{c}{feed down fractions [\%]} \\ \hline')

    for par in ['r_chic1_jpsi', 'r_chic2_jpsi', 'r_chic_jpsi', 'r_psip_jpsi',
                'r_jpsi_direct', 'r_psip_chic1', 'r_psip_chic2']:
        rows.append(
            fmt_row(DERIV_PAR[par], print_val(data_fb, par, True),
                    print_val(data_def, par, True), print_val(data_eb, par, True))
        )

    rows.append(r'\hline\hline')
    rows.append(r'\multicolumn{4}{c}{$\lambda_{\vartheta}$} \\ \hline')

    for par in ['lth_jpsi', 'lth_psip', 'lth_chic1', 'lth_chic2', 'lth_jpsi_chic']:
        rows.append(
            fmt_row(DERIV_PAR[par], print_val(data_fb, par),
                    print_val(data_def, par), print_val(data_eb, par))
        )

    return create_table(rows, 'l c c c')


def main(args):
    """Main"""
    vals_free_b = get_vals_uncers('/'.join([args.resultdir, 'fit_results_free_betas.root']),
                                  get_par_to_idx(0))
    vals_default = get_vals_uncers('/'.join([args.resultdir, 'fit_results_equal_total_cs.root']),
                                   get_par_to_idx(1))
    vals_equal_b = get_vals_uncers('/'.join([args.resultdir, 'fit_results_equal_betas.root']),
                                   get_par_to_idx(2))

    data_default = get_dataframe('/'.join([args.resultdir, 'scan_ptm_5_2500k_equal_total_cs.root']))
    data_equal_b = get_dataframe('/'.join([args.resultdir, 'scan_ptm_5_2500k_equal_betas.root']))
    data_free_b = get_dataframe('/'.join([args.resultdir, 'scan_ptm_5_2500k_free_betas.root']))

    jpsi_dir = lambda d: 1 - (d.r_chic2_jpsi + d.r_chic1_jpsi + d.r_psip_jpsi)

    data_default.loc[:, 'r_jpsi_direct'] = jpsi_dir(data_default)
    data_equal_b.loc[:, 'r_jpsi_direct'] = jpsi_dir(data_equal_b)
    data_free_b.loc[:, 'r_jpsi_direct'] = jpsi_dir(data_free_b)

    with open_tex_file(args.outfile, PREAMBLE) as texf:
        texf.write(
            create_fit_par_table(vals_free_b, vals_default, vals_equal_b)
        )

        texf.write('\n\n\n')

        texf.write(
            create_deriv_par_table(data_free_b, data_default, data_equal_b)
        )



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script to generate an overview'
                                     ' table over the three scenarios after a '
                                     'run of run_all_options.sh')
    parser.add_argument('resultdir', help='Directon into which run_all_options.sh'
                        ' has stored its results')
    parser.add_argument('-o', '--outfile', help='Output .tex file',
                        default='fit_parameters_table.tex')

    clargs = parser.parse_args()
    main(clargs)
