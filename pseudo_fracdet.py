#!/usr/bin/env python
"""
Pseudo-fracdet map
"""
__author__ = "Sidney Mau"

import os
import glob
import yaml
import numpy as np
import healpy as hp
import fitsio as fits


with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['data']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']


############################################################

infiles = glob.glob ('{}/*.fits'.format(datadir))

nside = 2048
pix = []
for infile in infiles:
    print('loading {}'.format(infile))
    data = fits.read(infile, columns=['RA','DEC'])
    p = hp.ang2pix(nside, data['RA'], data['DEC'], lonlat=True)
    #pix.append(p)
    pix.append(np.unique(p))

print('Constructing map')
pix = np.concatenate(pix)
pix = np.unique(pix)
coverage_map = np.zeros(hp.nside2npix(nside))
coverage_map[pix] = 1

print('Writing output')
result = '{}_pseudo_fracdet.fits.gz'.format(survey)
hp.write_map(result, coverage_map)
