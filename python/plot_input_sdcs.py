#!/usr/bin/env python3

import os
import sys
sys.path.append(os.path.expandvars('${CHIB_CHIC_POLFW_DIR}/python'))
import glob
from copy import deepcopy

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

# ROOT can't figure out the scalar_op template implementation yet to be able to
# directly do this in python. So easiest solution: Simply define the
# functionality inside a function that we can then call from python
AUXILIARY_CPP_FUNCS = r'''
sdc::SDC lth_from_flong(const sdc::SDC& fLong) {
  return (1 - 3 * fLong) / (1 + fLong);
}

sdc::SDC scale_sdc(double factor, sdc::SDC& sdc) {
  return factor * sdc;
}
'''

# Load sdcs and some reading functionality
r.gInterpreter.AddIncludePath(f'{THIS_DIR}/../interface')
r.gInterpreter.LoadFile(f'{THIS_DIR}/../src/read_sdcs.h')
r.gInterpreter.LoadFile(f'{THIS_DIR}/../interface/nrqcd_helpers.h')
r.gInterpreter.LoadText(AUXILIARY_CPP_FUNCS)
from ROOT import sdc
from ROOT import readPsiSDCs, readChic1SDCs, readChic2SDCs, lth_from_flong, scale_sdc
XRANGE = (5, 120) # Same as in PLB 773 (2017) 476
X_AXIS_UNIT = 'p_{T}'

# Load plotting helper
from utils.plot_helpers import mkplot, default_attributes, setup_legend
from utils.misc_helpers import cond_mkdir
from utils.setup_plot_style import set_basic_style

# Consistent attributes for plots involving multiple states
ATTRS = {
    '3S1_1': [default_attributes(size=0.5)[0]],
    '3S1_8': [default_attributes(size=0.5)[1]],
    '3PJ_8': [default_attributes(size=0.5)[2]],
    '1S0_8': [default_attributes(size=0.5)[3]],
    '3P1_1': [default_attributes(size=0.5)[4]],
    '3P2_1': [default_attributes(size=0.5)[4]],
}

# Scale factors for comparability with the Shao SDCs from PLB 773 (2017) 476
RAP_SCALE = 1 / 2.4 * 1.5
SCALE_FACTORS = {
    '3S1_1': RAP_SCALE * 6,
    '3S1_8': RAP_SCALE,
    '3PJ_8': RAP_SCALE,
    '1S0_8': RAP_SCALE,
    '3P1_1': RAP_SCALE * 6,
    '3P2_1': RAP_SCALE * 6,
}

# Attributes for plotting comparison plots of total and longitudinal SDCs
ATTRS_LT = [
    deepcopy(default_attributes(size=0.5)[0]),
    deepcopy(default_attributes(size=0.5)[1]),
    default_attributes(size=0.5, open_markers=False)[0],
    default_attributes(size=0.5, open_markers=False)[1],
]

def get_order(order):
    """Translate the string into the enum from c++"""
    if order == 'LP' or order == 'LP+NLO':
        return sdc.SDCType.LP_NLO
    if order == 'NLO':
        return sdc.SDCType.NLO
    if order == 'LO':
        return sdc.SDCType.LO

    raise argparse.ArgumentTypeError("type has to be one of 'LP', 'NLO' or 'LO'")

def get_order_str(order):
    """The other way around"""
    if order == sdc.SDCType.LP_NLO: return 'LP_NLO'
    if order == sdc.SDCType.NLO: return 'NLO'
    if order == sdc.SDCType.LO: return 'LO'

    return ''


def plot_all_sdcs(sdcdir, order, outdir):
    """Plot all sdcs that are present in the sdcdir"""
    all_sdcs = glob.glob(f'{sdcdir}/*.txt')
    for fn in all_sdcs:
        name = os.path.basename(fn).replace('.txt', '')
        s = sdc.read_from_file(fn, order)
        s.scale_support(1 / r.ccbarMassSDCs) # rescale to pT/M
        can = mkplot(s.asTGraph(), drawOpt='P',
                     xRange=XRANGE, xLabel=X_AXIS_UNIT,
                     yLabel=name)
        can.SaveAs(f'{outdir}/{name}.pdf')

        graphs = [g for g in s.asAlwaysPositiveGraphs() if g.GetN() > 0]
        can = mkplot(graphs,
                     drawOpt='P', logy=True,
                     xRange=XRANGE, xLabel=X_AXIS_UNIT,
                     legOpt='P', legPos=(0.75, 0.8, 0.88, 0.88), legEntries=('original', 'flipped'),
                     ydscale=0.1, yLabel=f'{name} (negative values flipped)')
        can.SaveAs(f'{outdir}/{name}_log_positive.pdf')


def _plot_state_sdcs(sdcdir, order, outdir, base_name, state_enum, states, read_f):
    """Plot the total SDCs and the lth of all intermediate states for a given
    quarkonium state"""
    sdcs = read_f(sdcdir, order)

    def _plot_total_sdc(state, **kwargs):
        index = getattr(state_enum, f's{state}')
        sdc_tot = sdcs.tot[index]
        sdc_tot = scale_sdc(SCALE_FACTORS[state], sdc_tot)
        graphs = [g for g in sdc_tot.asAlwaysPositiveGraphs() if g.GetN() > 0]
        can = mkplot(graphs[0], attr=ATTRS[state], legEntries=[state], legOpt='L', **kwargs)
        if len(graphs) > 1:
            can = mkplot(graphs[1], can=can, drawOpt=kwargs.get('drawOpt', 'C') + 'same',
                         attr=ATTRS[state])
        return can

    leg = setup_legend(0.75, 0.75, 0.92, 0.91)
    can = _plot_total_sdc(states[0], drawOpt='C', leg=leg,
                           xRange=XRANGE, logx=True, xLabel=X_AXIS_UNIT,
                           logy=True, yLabel='total SDCs', yRange=[1e-3, 0.5e5])

    for state in states[1:]:
        can = _plot_total_sdc(state, drawOpt='C same', leg=leg, can=can)

    can.SaveAs(f'{outdir}/{base_name}_total_SDCs_{get_order_str(order)}.pdf')

    def _plot_lth(state, **kwargs):
        index = getattr(state_enum, f's{state}')
        f_long = sdcs.lng[index] / sdcs.tot[index]
        lth = lth_from_flong(f_long)
        return mkplot(lth.asTGraph(),
                      attr=ATTRS[state], legEntries=[state], legOpt='L', **kwargs)

    leg = setup_legend(0.75, 0.8, 0.92, 0.94)
    can = _plot_lth(states[0], drawOpt='C', leg=leg,
                    xRange=XRANGE, xLabel=X_AXIS_UNIT, logx=True,
                    yRange=[-2, 2], yLabel='#lambda_{#vartheta}^{#psi}')
    for state in states[1:]:
        can = _plot_lth(state, drawOpt='C same', leg=leg, can=can)

    can.SaveAs(f'{outdir}/{base_name}_lth_SDCs_{get_order_str(order)}.pdf')

    def _plot_tot_long_signs(state):
        index = getattr(state_enum, f's{state}')
        sdc_long, sdc_tot = sdcs.lng[index], sdcs.tot[index]
        leg = setup_legend(0.75, 0.8, 0.93, 0.94)
        can = mkplot([g for g in sdc_tot.asAlwaysPositiveGraphs() if g.GetN() > 0],
                     drawOpt='P', logy=True, attr=ATTRS_LT[:2],
                     xRange=XRANGE, xLabel=X_AXIS_UNIT,
                     legOpt='P', leg=leg, legEntries=('tot +', 'tot -'),
                     ydscale=0.1, yLabel=f'{state} SDCs')

        return mkplot([g for g in sdc_long.asAlwaysPositiveGraphs() if g.GetN() > 0],
                      can=can, drawOpt='P same', attr=ATTRS_LT[2:],
                      legOpt='P', leg=leg, legEntries=('long +', 'long -'))

    for state in states:
        can = _plot_tot_long_signs(state)
        can.SaveAs(f'{outdir}/{base_name}_tot_long_SDCs_{state}_{get_order_str(order)}.pdf')


def plot_psi_sdcs(sdcdir, order, outdir):
    """Plot the relevant psi sdcs"""
    _plot_state_sdcs(sdcdir, order, outdir, 'psi', r.PsiSDCs,
                     ['3S1_8', '3PJ_8', '1S0_8', '3S1_1'], readPsiSDCs)


def plot_chic1_sdcs(sdcdir, order, outdir):
    """Plot the relevant chic1 sdcs"""
    _plot_state_sdcs(sdcdir, order, outdir, 'chic1', r.Chic1SDCs, ['3P1_1', '3S1_8'], readChic1SDCs)


def plot_chic2_sdcs(sdcdir, order, outdir):
    """Plot the relevant chic1 sdcs"""
    _plot_state_sdcs(sdcdir, order, outdir, 'chic2', r.Chic2SDCs, ['3P2_1', '3S1_8'], readChic2SDCs)


def main(args):
    """Main"""
    set_basic_style()
    sdc_dir = args.inputdir
    cond_mkdir(args.outdir)
    if args.plot_all:
        plot_all_sdcs(sdc_dir, args.type, args.outdir)

    plot_psi_sdcs(sdc_dir, args.type, args.outdir)
    plot_chic1_sdcs(sdc_dir, args.type, args.outdir)
    plot_chic2_sdcs(sdc_dir, args.type, args.outdir)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to plot input SDCs that are used in the fit')
    parser.add_argument('-i', '--inputdir', help='Directory where the SDCs can be found',
                        default=f'{THIS_DIR}/../../SDCs_Chung')
    parser.add_argument('-o', '--outdir', help='Output directory where the plots will be stored',
                        default='.')
    parser.add_argument('-t', '--type', help='Type (or better order) of the SDCs (LP=LP+NLO, NLO, LO)',
                        default='LP', type=get_order)
    parser.add_argument('-a', '--plot-all', help='Plot all SDCs (i.e. the "raw" input ones) not just the ones that are used in the end',
                        action='store_true', default=False)

    clargs = parser.parse_args()
    main(clargs)
