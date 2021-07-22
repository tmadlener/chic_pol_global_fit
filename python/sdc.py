#!/usr/bin/env python3
"""
Module to make the sdc interface available via PyROOT bindings
"""

import os
THIS_DIR = os.path.dirname(os.path.realpath(__file__))

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch()

r.gInterpreter.AddIncludePath(f'{THIS_DIR}/../interface')
r.gInterpreter.LoadFile(f'{THIS_DIR}/../interface/sdc.h')

from ROOT import sdc
