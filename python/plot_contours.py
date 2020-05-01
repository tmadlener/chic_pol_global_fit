#!/usr/bin/env python
"""
Script to make plots of contours (and best fit values)
"""
from collections import OrderedDict

import numpy as np

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.setup_plot_style import set_TDR_style

from common_helpers import parse_inputs

# define attributes for up to four different confidence levels per input
CONT_ATTRS = []
CONT_ATTRS.append(default_attributes(line=1, width=2))
CONT_ATTRS.append(default_attributes(line=7, width=2))
CONT_ATTRS.append(default_attributes(line=2, width=2))
CONT_ATTRS.append(default_attributes(line=5, width=2))

# attributes for best fit
BF_ATTRS = default_attributes(size=1.5, open_markers=False)
for att in BF_ATTRS:
    att['marker'] = 34

# legend attributes for the contours and the best fit
LEG_C_ATTR = [
    {'color': 12, 'width': 2, 'line': 1},
    {'color': 12, 'width': 2, 'line': 7},
    {'color': 12, 'width': 2, 'line': 2},
    {'color': 12, 'width': 2, 'line': 5}
]
LEG_BF_ATTR = [{'color': 12, 'marker': 34, 'size': 1.5}]


def get_contours(files, conf_level, var_1, var_2):
    """Get the contours of the corresponding conf_level from all files"""
    cont_name = '_'.join(['contour', var_1, var_2,
                          str(conf_level).replace('.', 'p')])

    return conf_level * 100, [f.Get(cont_name) for f in files.values()]


def get_best_fits(files, var_1, var_2):
    """Get the best fit points from all files"""
    name = '_'.join(['best', 'fit', var_1, var_2])
    return [f.Get(name) for f in files.values()]


def get_all_contours(files, conf_levels, var_1, var_2):
    """Get all the contours for all confidence levels from all the files in a
    list of lists form"""
    contours = OrderedDict()

    for conf_l in conf_levels:
        cl, conts = get_contours(files, conf_l, var_1, var_2)
        contours[cl] = conts

    return contours


def make_plot(contours, best_fits, leg_keys, **kwargs):
    """Make the plot and return the canvas"""
    conts = contours.values()
    conf_levels = contours.keys()

    leg_coords = (0.2, 0.6, 0.4, 0.8)
    if len(conf_levels) > 1:
        leg_coords = (0.2, 0.6, 0.4, 0.84)

    leg = setup_legend(*leg_coords)


    can = mkplot(conts[0], drawOpt='L', leg=leg, legEntries=leg_keys,
                 attr=CONT_ATTRS[0], legOpt='L', **kwargs)

    if len(conts) > 1:
        leg_graphs = [r.TGraph(1, np.array([-1000]), np.array([-10000]))
                      for _ in conf_levels]
        mkplot(leg_graphs, drawOpt='Lsame', can=can, attr=LEG_C_ATTR,
               leg=leg, legEntries=['{:.1f}% CL'.format(v) for v in conf_levels],
               legOpt='L')

        for i in range(1, len(contours)):
            mkplot(conts[i], can=can, drawOpt='Lsame', attr=CONT_ATTRS[i])

    # best fit points
    mkplot(best_fits, drawOpt='Psame', can=can, attr=BF_ATTRS)
    leg_graph = r.TGraph(1, np.array([-1000]), np.array([-10000]))
    mkplot(leg_graph, drawOpt='Psame', can=can, leg=leg, legEntries=['best fit'],
           attr=LEG_BF_ATTR, legOpt='P')


    return can


def main(args):
    """Main"""
    set_TDR_style()
    files = parse_inputs(args.input, enforce_keys=True)

    var_x = 'lambda_1'
    var_y = 'lambda_2'

    contours = get_all_contours(files, args.conf_levels, var_x, var_y)
    best_fits = get_best_fits(files, var_x, var_y)

    can = make_plot(contours, best_fits, files.keys(),
                    xRange=[-2, 2], yRange=[-2, 2],
                    xLabel='#lambda_{#vartheta}^{#chi_{c1}}',
                    yLabel='#lambda_{#vartheta}^{#chi_{c2}}')

    can.SaveAs(args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot the fit results')
    parser.add_argument('-i', '--input', help='Legend entry and input root file '
                        'separated by a colon for each result that should be on '
                        'the plot', action='append')
    parser.add_argument('-cl', '--conf-levels', help='The CL levels (0, 1) for '
                        'which the contour should be plotted', default='0.68',
                        type=lambda s: [float(v) for v in s.split(',')])
    parser.add_argument('-o', '--outfile', help='Name of the produced output '
                        'file')


    clargs = parser.parse_args()
    main(clargs)
