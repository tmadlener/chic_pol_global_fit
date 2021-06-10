#!/usr/bin/env python
"""
Script that plots all data to best fit model comparisons
"""
import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import (
    mkplot, setup_legend, default_colors, _setup_canvas
)
from utils.misc_helpers import cond_mkdir, get_vals_from_rwbuffer
from utils.setup_plot_style import set_basic_style
from utils.data_handling import list_obj
from utils.graph_utils import scale_graph, pull_graph, divide_func, shift_graph


COL = default_colors()

DATA_ATTR = [{'color': 1, 'marker': 20, 'size': 0.75},
             {'color': 1, 'marker': 25, 'size': 0.75},
             {'color': 1, 'marker': 26, 'size': 0.75}]

CS_RANGE = [2, 35]
CS_LABEL = 'd#sigma/dp_{T} (nb/GeV)'
POL_LABEL = '#lambda_{#vartheta}'
PTM_LABEL = 'p_{T}/M'
CTH_LABEL = '|cos#vartheta^{HX}|'
POL_RANGE = [0, 1]

B_CHIC2_JPSI = 19.0e-2
B_CHIC1_JPSI = 34.3e-2

# define colors for J/psi and psi(2S) polarization plot
JPSI_COL = r.TColor.GetFreeColorIndex()
TCOLOR_JPSI = r.TColor(JPSI_COL, 13./255, 154./255, 0, "")
PSI_COL = r.TColor.GetFreeColorIndex()
TCOLOR_PSI = r.TColor(PSI_COL, 178./255, 33./255, 179./255, "")


def _pull_plot(graphs, model, **kwargs):
    """Make the pull plot"""

    pull_graphs = [pull_graph(g, model) for g in graphs]

    can = mkplot([r.TLine(CS_RANGE[0], v, CS_RANGE[1], v) for v in [0, 2, -2]],
                 drawOpt='L', attr=[{'color': 12, 'line': 1, 'width': 2},
                                    {'color': 12, 'line': 7, 'width': 2},
                                    {'color': 12, 'line': 7, 'width': 2}],
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yRange=[-4, 4], logx=True,
    yLabel='(data - model) / data uncertainty')

    mkplot(pull_graphs, drawOpt='P same', can=can, **kwargs)

    return can


def _rel_diff_plot(graphs, model, **kwargs):
    """Make a plot with the relative differences"""
    # first assign the model values as uncertainties to the graphs
    delta_graphs = []
    for g in graphs:
        x_vals = get_vals_from_rwbuffer(g.GetX(), g.GetN())
        y_vals = get_vals_from_rwbuffer(g.GetY(), g.GetN())
        fy_vals = np.array([model.Eval(x) for x in x_vals])
        dy_vals = (y_vals - fy_vals)
        delta_graphs.append(divide_func(
            shift_graph(g, -fy_vals), model
        ))
        # delta_graphs.append(type(g)(len(x_vals), x_vals, dy_vals, *get_errors(g)))

    can = mkplot(delta_graphs, drawOpt='PE',
                 xRange=CS_RANGE, xLabel='p_{T}/M', logx=True,
                 yRange=[-0.5, 0.5],
                 yLabel='(data - model) / model', **kwargs)

    mkplot(r.TLine(CS_RANGE[0], 0, CS_RANGE[1], 0), can=can, drawOpt='L same',
           attr=[{'color': 12, 'line': 1, 'width': 2}])

    return can


def _get_present(gfile, gnames, lentries):
    """Get the graphs and leg_entries that are present in the file"""
    graphs, leg_entries = [], []
    for name, entry in zip(gnames, lentries):
        g = gfile.Get(name)
        if g:
            graphs.append(g)
            leg_entries.append(entry)

    return graphs, leg_entries


def psi2S_cs(gfile):
    """psi(2S) cross section plot"""
    names = ['psi2S_CMS_cs', 'psi2S_ATLAS_cs']
    leg_entries = ['#psi(2S) #rightarrow #mu#mu CMS',
                             '#psi(2S) #rightarrow J/#psi #pi#pi ATLAS']

    graphs, leg_entries = _get_present(gfile, names, leg_entries)


    leg = setup_legend(0.5, 0.7, 0.88, 0.86)
    can = mkplot(graphs, drawOpt='PE', logy=True, logx=True, attr=DATA_ATTR,
                 leg=leg,
                 legEntries=leg_entries,
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yLabel=CS_LABEL)

    mkplot(gfile.Get('psip_cs_direct'), can=can, drawOpt='Lsame',
           leg=leg, legEntries=['best fit'], legOpt='L')

    return can


def psi2S_cs_pulls(gfile):
    """psi(2S) cross section pulls"""
    names = ['psi2S_CMS_cs', 'psi2S_ATLAS_cs']
    leg_entries = ['#psi(2S) #rightarrow #mu#mu CMS',
                             '#psi(2S) #rightarrow J/#psi #pi#pi ATLAS']

    graphs, leg_entries = _get_present(gfile, names, leg_entries)

    model = gfile.Get('psip_cs_direct')

    return _pull_plot(graphs, model, legPos=(0.2, 0.8, 0.4, 0.92), legOpt='P',
                      legEntries=leg_entries)


def psi_2S_rel_diff(gfile):
    """psi(2S) relative diff plot"""
    names = ['psi2S_CMS_cs', 'psi2S_ATLAS_cs']
    leg_entries = ['#psi(2S) #rightarrow #mu#mu CMS',
                             '#psi(2S) #rightarrow J/#psi #pi#pi ATLAS']

    graphs, leg_entries = _get_present(gfile, names, leg_entries)
    model = gfile.Get('psip_cs_direct')

    return _rel_diff_plot(graphs, model, legEntries=leg_entries, legOpt='PE',
                          legPos=(0.2, 0.2, 0.4, 0.32))


def jpsi_cs(gfile):
    """jpsi cross section plot"""
    graphs, leg_entries = _get_present(gfile, ['jpsi_CMS_cs'],
                                       ['J/#psi CMS'])

    leg = setup_legend(0.6, 0.7, 0.88, 0.82)
    can = mkplot(graphs,
                 drawOpt='PE', logy=True, logx=True, attr=DATA_ATTR,
                 leg=leg, legEntries=leg_entries,
                 xRange=CS_RANGE, xLabel=PTM_LABEL, yLabel=CS_LABEL)

    mkplot(gfile.Get('jpsi_cs_full'), can=can, drawOpt='Lsame',
           leg=leg, legEntries=['best fit'], legOpt='L')

    return can


def jpsi_cs_pulls(gfile):
    """jpsi cross section pulls"""
    graphs, leg_entries = _get_present(gfile, ['jpsi_CMS_cs'],
                                       ['CMS'])
    model = gfile.Get('jpsi_cs_full')

    return _pull_plot(graphs, model, legPos=(0.2, 0.8, 0.4, 0.92), legOpt='P',
                      legEntries=leg_entries)

def jpsi_rel_diff(gfile):
    """jpsi relative difference"""
    """jpsi cross section pulls"""
    graphs, leg_entries = _get_present(gfile, ['jpsi_CMS_cs'],
                                       ['CMS'])
    model = gfile.Get('jpsi_cs_full')

    return _rel_diff_plot(graphs, model, legEntries=leg_entries, legOpt='P',
                          legPos=(0.2, 0.2, 0.4, 0.32))


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

    can = mkplot(scale_graph(r.TGraph(gfile.Get('chic_ratio_cs_full')), B_CHIC2_JPSI / B_CHIC1_JPSI),
                 drawOpt='L', attr=[{'color': 1, 'line': 1, 'width': 2}],
                 yRange=[0, 1], xRange=[2, 10], xLabel=PTM_LABEL, yLabel='ratio',
                 leg=leg, legEntries=['model'], legOpt='L')


    mkplot(scale_graph(gfile.Get('chic_ratio_CMS_cs'), B_CHIC2_JPSI / B_CHIC1_JPSI),
           drawOpt='PE same', can=can,
           attr=[{'color': COL[0], 'marker': 20, 'size': 1.8}],
           leg=leg, legEntries=['B #chi_{c2} / B #chi_{c1} CMS'])


    return can


def psi_pol(gfile):
    """psi(2S) and J/psi polarization plot"""
    leg = setup_legend(0.6, 0.2, 0.88, 0.36)


    can = mkplot([gfile.Get(n) for n in ['psip_pol_direct', 'jpsi_pol_full']],
                 drawOpt='L', colors=[PSI_COL, JPSI_COL],
                 leg=leg, legEntries=['#psi(2S) model', 'J/#psi model'], legOpt='L',
                 xRange=[2, 30], xLabel=PTM_LABEL, yLabel=POL_LABEL, yRange=[-1, 1])

    mkplot([gfile.Get(n) for n in ['psi2S_CMS_pol', 'jpsi_CMS_pol']],
                 attr=[{'color': PSI_COL, 'marker': 20, 'size': 0.75},
                       {'color': JPSI_COL, 'marker': 21, 'size': 0.75}],
           can=can, drawOpt='PE same', leg=leg,
           legEntries=['#psi(2S) CMS', 'J/#psi CMS'])

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


def combined_cs(gfile, only=None):
    """Figure combining all cross sections into one figure"""
    MSIZE = 0.75
    graph_attrs = {
        'jpsi_CMS_cs': {'color': COL[0], 'marker': 24, 'size': MSIZE},
        'psi2S_CMS_cs': {'color': COL[0], 'marker': 25, 'size': MSIZE},
        'psi2S_ATLAS_cs': {'color': COL[1], 'marker': 26, 'size': MSIZE},
        'chic1_ATLAS_cs': {'color': COL[1], 'marker': 22, 'size': MSIZE},
        'chic2_ATLAS_cs': {'color': COL[1], 'marker': 23, 'size': MSIZE},
    }


    if only is not None:
        # No need for leg_entries here
        graphs, _ = _get_present(gfile, [g for g in graph_attrs.keys() if only in g],
                                 [g for g in graph_attrs.keys() if only in g])
        attrs = [graph_attrs[g] for g in graph_attrs.keys() if only in g]
    else:
        graphs, _ = _get_present(gfile, graph_attrs.keys(), graph_attrs.keys())
        attrs = list(graph_attrs.values())

    can = mkplot(graphs, attr=attrs,
                 drawOpt='PE', xRange=[2, 50], yRange=[1e-5, 30], logy=True, logx=True,
                 xLabel=PTM_LABEL, yLabel=CS_LABEL)

    if only is None:
        mkplot([gfile.Get(f) for f in  ['psip_cs_direct', 'jpsi_cs_full', 'chic1_cs_full', 'chic2_cs_full']],
               can=can, drawOpt='L same', attr=[{'color': 1, 'width': 1, 'line': 1}])


    leg_graphs = [r.TGraph(1, np.array([10000.]), np.array([10000.])) for _ in range(5)]
    mkplot(leg_graphs, can=can, drawOpt='P same', legPos=(0.2, 0.2, 0.4, 0.40),
           legEntries=['J/#psi', '#psi(2S) #rightarrow #mu#mu', '#psi(2S) #rightarrow J/#psi #pi#pi',
                       '#chi_{c1}', '#chi_{c2}'], legOpt='P',
           attr=[graph_attrs[n] for n in ['jpsi_CMS_cs', 'psi2S_CMS_cs', 'psi2S_ATLAS_cs',
                                          'chic1_ATLAS_cs', 'chic2_ATLAS_cs']])

    return can


def main(args):
    """Main"""
    set_basic_style()
    gfile = r.TFile.Open(args.graphfile)


    cond_mkdir(args.outdir)

    for pfunc in [psi2S_cs, jpsi_cs, chic_cs, chic_ratio_cs, psi_pol, combined_cs,
                  psi2S_cs_pulls, jpsi_cs_pulls, psi_2S_rel_diff, jpsi_rel_diff]:
        can = pfunc(gfile)
        can.SaveAs('/'.join([args.outdir, pfunc.__name__ + '.pdf']))

    for ptm, can in costh_ratios(gfile):
        can.SaveAs('/'.join([args.outdir, 'costh_ratio_ptm_{}.pdf'.format(ptm)]))

    can = combined_cs(gfile, 'ATLAS')
    can.SaveAs('/'.join([args.outdir, 'combined_cs_ATLAS']) + '.pdf')

    can = combined_cs(gfile, 'CMS')
    can.SaveAs('/'.join([args.outdir, 'combined_cs_CMS']) + '.pdf')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot best fit models'
                                     ' on top of (corrected) data graphs')
    parser.add_argument('graphfile', help='File containing the graphs and models')
    parser.add_argument('-o', '--outdir', help='directory into which the plots '
                        'will be stored', default='.')

    clargs = parser.parse_args()
    main(clargs)
