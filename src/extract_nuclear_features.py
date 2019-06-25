#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190624
# 
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import re
import sys

from ggc.args import check_threads
from joblib import Parallel, delayed
from skimage.io import imread
from tqdm import tqdm

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
Extract nuclear features from nuclear boxes, generating using the -u flag
of the gpseq_anim script. Use also --nuclear-sel with no parameters to extract
all nuclear boxes (useful if new selection should be performed).
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('rootdir', type = str, help = '''
	Path to root directory with nuclear boxes.''')
parser.add_argument('outdir', type = str, help = '''
	Path to output directory where the output should be written to.''')

# Add arguments with default value
parser.add_argument('-P', '--prefix', type = str, help = '''
	Prefix for output table name. Defaults to output directory name.''')
parser.add_argument('-a', '--aspect', type = float, nargs = 3,
	help = """Physical size of Z, Y and X voxel sides in nm. Default: 300.0 216.6 216.6""",
	metavar = ('Z', 'Y', 'X'), default = [300., 216.6, 216.6])
parser.add_argument('-t', '--threads', metavar = 'nthreads', type = int,
	default = 1, help = """Number of threads to be used for parallelization.""")

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

assert os.path.isdir(args.rootdir)
assert os.path.isdir(args.outdir)

if type(None) == type(args.prefix):
	args.prefix = os.basename(args.outdir)

args.threads = check_threads(args.threads)

# FUNCTIONS ====================================================================

def get_nucleus_feats(signature, condition):
	sid, nid = [int(x) for x in re.match(r's([0-9]+)n([0-9]+)', signature).groups()]

	mask = gp.tools.image.add_top_bottom_slides(
		imread(os.path.join(args.rootdir, f'{signature}{condition}.tif')))
	idna = gp.tools.image.add_top_bottom_slides(
		imread(os.path.join(args.rootdir, f'{signature}{condition}.dna.tif')))
	isig = gp.tools.image.add_top_bottom_slides(
		imread(os.path.join(args.rootdir, f'{signature}{condition}.sig.tif')))

	nARadius_um = (mask.max(0).sum()*pxArea/np.pi)**.5
	nARadius_vx = (mask.max(0).sum()/np.pi)**.5
	nRadius_um = (mask.sum()*vxVolume*3/4/np.pi)**(1/3)
	nRadius_vx = (mask.sum()*3/4/np.pi)**(1/3)
	nSurface = gp.tools.image.calc_surface(mask, args.aspect)

	nData = {
		'condition' : [condition],						# Dataset
		'sid' : [sid], 'nid' : [nid],					# Series/Nucleus ID
		'area_px' : mask.max(0).sum(),					# Z-projection area in px
		'area_um2' : mask.max(0).sum()*pxArea,			# Z-projection area in um2
		'a2radius_px' : [nARadius_um],					# Radius in um (from area)
		'a2radius_um' : [nARadius_vx],					# Radius in px (from area)
		'volume_vx' : [mask.sum()],						# Volume in vx
		'volume_um3' : [mask.sum()*vxVolume],			# Volume in um3
		'v2radius_px' : [nRadius_um],					# Radius in um (from volume)
		'v2radius_um' : [nRadius_vx],					# Radius in vx (from volume)
		'dna_sum' : [np.nansum(idna[mask != 0])],		# DNA intensity integral
		'dna_mean' : [np.nanmean(idna[mask != 0])],		# DNA intensity mean
		'sig_sum' : [np.nansum(isig[mask != 0])],		# Signal intensity integral
		'sig_mean' : [np.nanmean(isig[mask != 0])],		# Signal intensity mean
		'surface_um2' : [nSurface],						# Surface in um2
		'sphericity' : 4*np.pi*(nRadius_um**2)/nSurface	# Sphericity
	}

	return(pd.DataFrame.from_dict(nData))

# RUN ==========================================================================

flist = os.listdir(args.rootdir)

nuclei_list = set([f.split(".")[0] for f in flist])
nuclei_dict = {}
for n in nuclei_list:
	signature, condition = n.split("iSK")
	condition = "iSK" + condition
	if condition not in nuclei_dict.keys():
		nuclei_dict[condition] = [signature]
	else:
		nuclei_dict[condition].append(signature)

pxArea = np.prod(args.aspect[1:])*1e-6
vxVolume = np.prod(args.aspect)*1e-9 # (um3)

data = []
for condition in tqdm(list(nuclei_dict.keys()), desc = "Condition"):
	cData = Parallel(n_jobs = args.threads, verbose = 11)(
	    delayed(get_nucleus_feats)(signature, condition)
	    for signature in nuclei_dict[condition])
	data.append(pd.concat(cData))
outData = pd.concat(data).reset_index(drop = True)

outData.to_csv(os.path.join(args.outdir, f"{args.prefix}.nuclear_features.tsv"),
	sep = "\t",index = False, header = True)

# END ==========================================================================

################################################################################
