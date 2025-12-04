#!/usr/bin/env python
# Copyright 2025, the PyITM Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

"""Routines to perform statistical comparisons."""

import numpy as np

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_normalization(vals, doReport = False):
    norm = 1.0
    m = np.mean(np.abs(vals))
    n = int(np.floor(np.log10(m)))
    if (np.abs(n) >= 2):
        norm = 10**n
    if (doReport):
        print('Normalization Factor : ', norm)
    return norm

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_mean_diff(vals1, vals2, doReport = False):
    m = np.mean(vals1 - vals2)
    if (doReport):
        print(' Mean Diff (v1-v2) : ', m)
    return m

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_normMean_mean_diff(vals1, vals2, doReport = False):
    m = np.mean(vals1 - vals2)/np.mean(vals2)
    if (doReport):
        print(' Norm (to mean v2) Mean Diff (v1-v2) : ', m)
    return m

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_normStd_mean_diff(vals1, vals2, doReport = False):
    m = np.mean(vals1 - vals2)/np.std(vals2)
    if (doReport):
        print(' Norm (to std v2) Mean Diff (v1-v2) : ', m)
    return m

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_rms_diff(vals1, vals2, doReport = False):
    m = np.sqrt(np.mean((vals1 - vals2)**2))
    if (doReport):
        print(' RMS Diff : ', m)
    return m

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_normMean_rms_diff(vals1, vals2, doReport = False):
    m = calc_rms_diff(vals1, vals2)
    m = m / np.mean(vals2)
    if (doReport):
        print(' Norm (to mean v2) RMS Diff : ', m)
    return m

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_normStd_rms_diff(vals1, vals2, doReport = False):
    m = calc_rms_diff(vals1, vals2)
    m = m / np.std(vals2)
    if (doReport):
        print(' Norm (to std) RMS Diff : ', m)
    return m

