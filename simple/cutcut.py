#!/usr/bin/env python

import numpy as np
import healpy as hp
import ugali.utils.healpix
import ugali.utils.projector
import glob
import argparse
import yaml
import subprocess

import percent

p = argparse.ArgumentParser()
p.add_argument('delta_x', type=float)
p.add_argument('bin_edge', type=float)
args = p.parse_args()


with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.full_load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']

cutcut_dir = '{}/cutcut_arrays_dx={}_binedge={}'.format(datadir, round(args.delta_x,3), round(args.bin_edge,1))
subprocess.call('mkdir -p {}'.format(cutcut_dir).split())
#infiles = glob.glob('{}/*.fits'.format(datadir))
#
#print('Pixelizing...')
#pix_nside = [] # Equatorial coordinates, RING ordering scheme
#for infile in infiles:
#    pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))
with open('{}/pix_nside.txt'.format(datadir)) as f:
    pix_nside = map(int, f.read().splitlines())

delta_x = args.delta_x # 0.003 seems ideal
area = delta_x**2
bins = np.arange(-args.bin_edge, args.bin_edge + 1.e-10, delta_x)
centers = 0.5 * (bins[0: -1] + bins[1:])

yy, xx = np.meshgrid(centers, centers)
xxflat, yyflat = xx.flatten(), yy.flatten()

i = 0
for pix_nside_select in pix_nside:
    ra, dec = ugali.utils.healpix.pixToAng(nside, pix_nside_select)
    proj = ugali.utils.projector.Projector(ra, dec)

    rara, decdec = proj.imageToSphere(xxflat, yyflat)
    cutcut = (ugali.utils.healpix.angToPix(nside, rara, decdec) == pix_nside_select).reshape(xx.shape)
    np.savez_compressed('{}/cutcut_array_{}'.format(cutcut_dir, pix_nside_select), cutcut)

    i += 1
    percent.bar(i, len(pix_nside))

