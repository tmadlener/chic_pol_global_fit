#!/usr/bin/env python
"""
Script for making 1d band plots as a function of pT/M
"""
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes
from utils.setup_plot_style import set_basic_style
from utils.misc_helpers import cond_mkdir

PTM_RANGE = [2, 10]
PTM_LABEL = 'p_{T}/M'

BAND_ATTR = default_attributes(open_markers=False)
for att in BAND_ATTR:
    att['fillalpha'] = (att['color'], 0.25)


PLOT_BANDS = {
    ('lth_psip', 'lth_jpsi_chic'): {
        'yRange': [-1.5, 2.0], 'yLabel': '#lambda_{#vartheta}^{J/#psi}',
        'legPos': (0.2, 0.2, 0.5, 0.28),
        'legEntries': ['direct', 'from #chi_{cJ}, J=1,2'],
        'pdfname': 'polarization_bands_jpsi.pdf'
    },

    ('lth_chic1', 'lth_chic2'): {
        'yRange': [-1.5, 2.0], 'yLabel': '#lambda_{#vartheta}',
        'legPos': (0.8, 0.86, 0.9, 0.94),
        'legEntries': ['#chi_{c1}', '#chi_{c2}'],
        'pdfname': 'polarization_bands_chi.pdf'
    },
    ('r_chic1_jpsi', 'r_chic2_jpsi', 'r_psip_jpsi', 'r_chic_jpsi'): {
        'yRange': [0, 0.45], 'yLabel': 'X #rightarrow J/#psi feed down',
        'legPos': (0.2, 0.78, 0.4, 0.94),
        'legEntries': ['#chi_{c1}', '#chi_{c2}', '#psi(2S)', '#chi_{c1} + #chi_{c2}'],
        'pdfname': 'feed_down_bands.pdf'
    }
}


def make_1d_band_plot(gfile, variables, **kwargs):
    """Make 1d band plot"""
    graphs = [
        gfile.Get('{}_v_ptm'.format(v)) for v in variables
    ]

    return mkplot(graphs, attr=BAND_ATTR, drawOpt='LE3',
                  xRange=PTM_RANGE, xLabel=PTM_LABEL, legOpt='FL',
                  **kwargs)


def main(args):
    """Main"""
    set_basic_style()

    outdir = '/'.join(['.', args.outdir])
    cond_mkdir(outdir)

    graphfile = r.TFile.Open(args.graphfile)

    for variables, kwargs in PLOT_BANDS.iteritems():
        can = make_1d_band_plot(graphfile, variables, **kwargs)
        can.SaveAs('/'.join([outdir, kwargs['pdfname']]))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for making pT/M band '
                                     'plots')
    parser.add_argument('graphfile', help='File containing the graphs of the '
                        'banbds')
    parser.add_argument('-o', '--outdir', help='Output directory into which the '
                        'plots are stored', default='plots_1d_band')


    clargs = parser.parse_args()
    main(clargs)
