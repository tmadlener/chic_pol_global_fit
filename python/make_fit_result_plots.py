#!/usr/bin/env python
"""
Script that plots all data to best fit model comparisons
"""
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, setup_legend, default_colors
from utils.misc_helpers import cond_mkdir
from utils.setup_plot_style import set_basic_style
from utils.data_handling import list_obj
from utils.graph_utils import scale_graph


COL = default_colors()

DATA_ATTR = [{'color': 1, 'marker': 20, 'size': 0.75},
             {'color': 1, 'marker': 25, 'size': 0.75}]

CS_RANGE = [2, 35]
CS_LABEL = 'd#sigma/dp_{T} (nb/GeV)'
POL_LABEL = '#lambda_{#vartheta}'
PTM_LABEL = 'p_{T}/M'
CTH_LABEL = '|cos#vartheta^{HX}|'
POL_RANGE = [0, 1]

B_CHIC2_JPSI = 19.0e-2
B_CHIC1_JPSI = 34.3e-2

def psi2S_cs(gfile):
    """psi(2S) cross section plot"""
    graphs = [gfile.Get(n) for n in ['psi2S_CMS_cs', 'psi2S_ATLAS_cs']]
    leg = setup_legend(0.6, 0.7, 0.88, 0.82)
    can = mkplot(graphs, drawOpt='PE', logy=True, attr=DATA_ATTR,
                 leg=leg, legEntries=['#psi(2S) {}'.format(e) for e in 'CMS', 'ATLAS'],
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yLabel=CS_LABEL)

    mkplot(gfile.Get('psi_cs_direct'), can=can, drawOpt='Lsame',
           leg=leg, legEntries=['best fit'], legOpt='L')

    return can


def jpsi_cs(gfile):
    """jpsi cross section plot"""
    leg = setup_legend(0.6, 0.7, 0.88, 0.78)
    can = mkplot(gfile.Get('jpsi_CMS_cs'), drawOpt='PE', logy=True, attr=DATA_ATTR,
                 leg=leg, legEntries=['J/#psi CMS'],
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yLabel=CS_LABEL)

    mkplot(gfile.Get('jpsi_cs_full'), can=can, drawOpt='Lsame',
           leg=leg, legEntries=['best fit'], legOpt='L')

    return can


def chic_cs(gfile):
    """chic cross sections plot"""
    leg = setup_legend(0.6, 0.7, 0.88, 0.86)
    can = mkplot([gfile.Get(n) for n in ['chic1_ATLAS_cs', 'chic2_ATLAS_cs']],
                 drawOpt='PE', logy=True, attr=DATA_ATTR,
                 leg=leg, legEntries=['#chi_{c1} ATLAS', '#chi_{c2} ATLAS'],
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yLabel=CS_LABEL,
                 yRange=[1e-3, 40])

    mkplot([gfile.Get(n) for n in ['chic1_cs_full', 'chic2_cs_full']],
           drawOpt='Lsame', can=can,
           leg=leg, legEntries=['best fit', 'best fit'], legOpt='L')

    return can


def chic_ratio_cs(gfile):
    """chic cross section ratio plot"""
    leg = setup_legend(0.6, 0.7, 0.88, 0.78)
    can = mkplot(scale_graph(gfile.Get('chic_ratio_CMS_cs'), B_CHIC2_JPSI / B_CHIC1_JPSI),
                 drawOpt='PE', attr=DATA_ATTR,
                 leg=leg, legEntries=['B #chi_{c2} / B #chi_{c1} CMS'],
                 yRange=[0, 1],
                 xRange=[2, 10], xLabel=PTM_LABEL, yLabel='ratio')

    mkplot(scale_graph(r.TGraph(gfile.Get('chic_ratio_cs_full')), B_CHIC2_JPSI / B_CHIC1_JPSI),
           can=can, drawOpt='Lsame',
           leg=leg, legEntries=['model'], legOpt='L')

    return can


def psi_pol(gfile):
    """psi(2S) and J/psi polarization plot"""
    leg = setup_legend(0.6, 0.2, 0.88, 0.36)
    can = mkplot([gfile.Get(n) for n in ['psi2S_CMS_pol', 'jpsi_CMS_pol']],
                 drawOpt='PE', attr=DATA_ATTR,
                 leg=leg, legEntries=['#psi(2S) CMS', 'J/#psi CMS'],
                 xRange=[2, 30], xLabel=PTM_LABEL, yLabel=POL_LABEL, yRange=[-1, 1])

    mkplot([gfile.Get(n) for n in ['psi_pol_direct', 'jpsi_pol_full']],
           can=can, drawOpt='Lsame', leg=leg, legOpt='L',
           legEntries=['#psi(2S) model', 'J/#psi model'])

    return can



def costh_ratios(gfile):
    """Make the costh ratio plots"""
    graphs = [gfile.Get(n) for n in list_obj(gfile, 'TGraph', 'chic_CMS_pol_ptM')]

    cans = []
    for graph in graphs:
        ptm_s = graph.GetName().split('_')[-1]
        ptm = float(ptm_s.replace('p', '.'))

        leg = setup_legend(0.55, 0.7, 0.88, 0.78)
        can = mkplot(graph, drawOpt='PE', attr=DATA_ATTR,
                     leg=leg, legEntries=['CMS @ p_{{T}} / M = {:.2f}'.format(ptm)],
                     xRange=[0, 1], xLabel=CTH_LABEL,
                     yRange=[0.2, 0.6], yLabel='#chi_{c2} / #chi_{c1}')

        model = gfile.Get('chic_pol_ptM_{}'.format(ptm_s))
        mkplot(model, can=can, drawOpt='Lsame',
               leg=leg, legEntries=['model'], legOpt='L')


        cans.append((ptm_s, can))

    return cans


def combined_cs(gfile):
    """Figure combining all cross sections into one figure"""
    MSIZE = 0.75
    graph_attrs = {
        'jpsi_CMS_cs': {'color': COL[0], 'marker': 24, 'size': MSIZE},
        'psi2S_CMS_cs': {'color': COL[0], 'marker': 25, 'size': MSIZE},
        'psi2S_ATLAS_cs': {'color': COL[1], 'marker': 25, 'size': MSIZE},
        'chic1_ATLAS_cs': {'color': COL[1], 'marker': 22, 'size': MSIZE},
        'chic2_ATLAS_cs': {'color': COL[1], 'marker': 23, 'size': MSIZE}
    }

    can = mkplot([gfile.Get(g) for g in graph_attrs.keys()], attr=graph_attrs.values(),
                 drawOpt='PE', xRange=[2, 50], yRange=[1e-5, 30], logy=True, logx=True,
                 xLabel=PTM_LABEL, yLabel=CS_LABEL)

    mkplot([gfile.Get(f) for f in  ['psi_cs_direct', 'jpsi_cs_full', 'chic1_cs_full', 'chic2_cs_full']],
           can=can, drawOpt='L same', attr=[{'color': 1, 'width': 1, 'line': 1}])


    leg_graphs = [r.TGraph(1, np.array([10000]), np.array([10000])) for _ in range(4)]
    mkplot(leg_graphs, can=can, drawOpt='P same', legPos=(0.2, 0.2, 0.4, 0.36),
           legEntries=['J/#psi', '#psi(2S)', '#chi_{c1}', '#chi_{c2}'], legOpt='P',
           attr=[graph_attrs[n] for n in ['jpsi_CMS_cs', 'psi2S_CMS_cs',
                                          'chic1_ATLAS_cs', 'chic2_ATLAS_cs']])

    return can


def main(args):
    """Main"""
    set_basic_style()
    gfile = r.TFile.Open(args.graphfile)


    cond_mkdir(args.outdir)

    for pfunc in [psi2S_cs, jpsi_cs, chic_cs, chic_ratio_cs, psi_pol, combined_cs]:
        can = pfunc(gfile)
        can.SaveAs('/'.join([args.outdir, pfunc.__name__ + '.pdf']))

    for ptm, can in costh_ratios(gfile):
        can.SaveAs('/'.join([args.outdir, 'costh_ratio_ptm_{}.pdf'.format(ptm)]))




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot best fit models'
                                     ' on top of (corrected) data graphs')
    parser.add_argument('graphfile', help='File containing the graphs and models')
    parser.add_argument('-o', '--outdir', help='directory into which the plots '
                        'will be stored', default='.')

    clargs = parser.parse_args()
    main(clargs)
