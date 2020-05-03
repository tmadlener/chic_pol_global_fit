#!/usr/bin/env python
"""
Script to make a correlation graph from a scan produced by run_ptM_scan assuming
that all points have the same ptm
"""
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.data_handling import get_dataframe
from utils.misc_helpers import flatten

from common_helpers import get_2d_hist, get_coverage_contour, select_phys_data


def process_input_variables(var_str):
    """Process the input variables"""
    var_pairs = []
    for vstr in var_str.split(','):
        var_pairs.append(vstr.split('|'))

    return var_pairs


def make_pair_correlation(data, varx, vary, phys_lam):
    """Fit the 2D distribution with a 2d normal distribution and get the
    1 sigma ellipse"""
    print('Processing pair: {}, {}'.format(varx, vary))
    data = select_phys_data(data, [varx, vary], phys_lam)
    dist_2d = get_2d_hist(data, varx, vary)

    cov_graphs = []
    for cov, sigma in [(0.683, 1), (0.955, 2), (0.997, 3)]:
        cov_graph = get_coverage_contour(dist_2d, cov)
        cov_graph.SetName('corr_2d_{}_{}_{}sigma'.format(varx, vary, sigma))

        cov_graphs.append(cov_graph)

    return cov_graphs


def main(args):
    """Main"""
    data = get_dataframe(args.scanfile, columns=list(flatten(args.variables)))

    outfile = r.TFile.Open(args.outfile, 'recreate')

    for varx, vary in args.variables:
        corr_graphs = make_pair_correlation(data, varx, vary, args.physical_lambdas)
        for cg in corr_graphs:
            cg.Write()

    outfile.Close()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to produce 2d '
                                     'correlation graphs from a run of '
                                     'run_ptM_scan')
    parser.add_argument('scanfile', help='File containing the scan results')
    parser.add_argument('-o', '--outfile', help='File to which the TGraphs will '
                        'be stored', default='correlation_graphs.root')
    parser.add_argument('-v', '--variables', help='Comma separated list of '
                        'variable pairs of the form VAR1|VAR2, for which the '
                        'graphs should be produced.',
                        default='r_chic1_jpsi|r_chic2_jpsi',
                        type=process_input_variables)
    parser.add_argument('-pl', '--physical-lambdas', help='Clip the lambdas of '
                        'the chic1 and chic2 to their physically allowed ranges '
                        'before obtaining the bands', action='store_true',
                        default=False)


    clargs = parser.parse_args()
    main(clargs)
