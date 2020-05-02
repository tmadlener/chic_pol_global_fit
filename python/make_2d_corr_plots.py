#!/usr/bin/env python
"""
Script for making the 2d correlation graph plots
"""
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.setup_plot_style import set_basic_style
from utils.misc_helpers import cond_mkdir

from common_helpers import get_label, get_range, parse_inputs

PLOT_PAIRS = (
    ('lth_jpsi', 'lth_chic1'),
    ('lth_jpsi', 'lth_chic2'),

    ('lth_chic1', 'lth_chic2'),

    ('lth_psip', 'lth_chic1'),
    ('lth_psip', 'lth_chic2'),

    ('dlth_jpsi_psip', 'lth_chic1'),
    ('dlth_jpsi_psip', 'lth_chic2'),

    ('r_chic1_jpsi', 'r_chic2_jpsi'),
    ('r_chic1_jpsi', 'lth_chic1'),
    ('r_chic2_jpsi', 'lth_chic2')
)

def flip_axis(graph):
    """Flip the x & and y coordinates of the passed graph"""
    return r.TGraph(graph.GetN(), graph.GetY(), graph.GetX())


def make_2d_contour_plot(input_files, var_pair):
    """Make one plot of all input_files"""
    varx, vary = var_pair
    graphs = [
        f.Get('corr_2d_{}_{}_1sigma'.format(varx, vary)) for f in input_files.values()
    ]

    # These are possibly stored in reversed order so the axis have to be
    # flipped
    if var_pair == ('r_chic1_jpsi', 'r_chic2_jpsi'):
        # try to get them the other way around
        graphs_2 = [
            f.Get('corr_2d_{}_{}_1sigma'.format(vary, varx)) for f in input_files.values()
        ]
        # Check which graphs are not present and replace them with the flipped axis version
        for i, g in enumerate(graphs):
            if not g:
                graphs[i] = flip_axis(graphs_2[i])


    if len(graphs) == 1:
        leg = None
        leg_entries = None
    else:
        leg = setup_legend(0.75, 0.8, 0.92, 0.8 + len(graphs) * 0.035)
        leg_entries = input_files.keys()

    can = mkplot(graphs, drawOpt='L', attr=default_attributes(line=1),
                 leg=leg, legEntries=leg_entries, legOpt='L',
                 yLabel=get_label(vary), yRange=get_range(vary),
                 xLabel=get_label(varx), xRange=get_range(varx))

    return can


def main(args):
    """Main"""
    set_basic_style()

    input_files = parse_inputs(args.input)
    outdir = '/'.join(['.', args.outdir])
    cond_mkdir(outdir)

    for var_pair in PLOT_PAIRS:
        can = make_2d_contour_plot(input_files, var_pair)
        save_name = outdir + '/' + 'corr_2d_{}_{}.pdf'.format(*var_pair)
        can.SaveAs(save_name)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make 2d correlation '
                                     'plots')
    parser.add_argument('-i', '--input', help='Legend entry and input root file '
                        'separaed by a colon for multiple inputs or just the '
                        'input root file for a single input', action='append')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'plots are stored', default='plots_2d_corr')


    clargs = parser.parse_args()
    main(clargs)
