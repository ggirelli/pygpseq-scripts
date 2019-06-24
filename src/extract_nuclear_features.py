#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190624
# Project: GPSeq-YFISH
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
parser.add_argument('-p', '--prefix', type = str, help = '''
	Prefix for output table name. Defaults to output directory name.''')
parser.add_argument('-a', '--aspect', type = float, nargs = 3,
	help = """Physical size of Z, Y and X voxel sides in nm. Default: 300.0 216.6 216.6""",
	metavar = ('Z', 'Y', 'X'), default = [300., 216.6, 216.6])

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

# FUNCTIONS ====================================================================

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
	cData = []
	for signature in tqdm(nuclei_dict[condition], desc = "Nuclei"):
		sid, nid = [int(x) for x in re.match(r's([0-9]+)n([0-9]+)',"s12n23").groups()]

		mask = gp.tools.image.add_top_bottom_slides(
			imread(os.path.join(args.rootdir, f'{signature}{condition}.tif')))
		idna = gp.tools.image.add_top_bottom_slides(
			imread(os.path.join(args.rootdir, f'{signature}{condition}.dna.tif')))
		isig = gp.tools.image.add_top_bottom_slides(
			imread(os.path.join(args.rootdir, f'{signature}{condition}.sig.tif')))

		nRadius = (mask.sum()*vxVolume*3/4/np.pi)**(1/3)
		nSurface = gp.tools.image.calc_surface(mask, args.aspect)

		nData = {
			'condition' : [condition],						# Dataset
			'sid' : [sid], 'nid' : [nid],					# Series/Nucleus ID
			'area_px' : mask.max(0).sum(),					# Z-projection area in px
			'area_um2' : mask.max(0).sum()*pxArea,			# Z-projection area in um2
			'volume_vx' : [mask.sum()],						# Volume in vx
			'volume_um3' : [mask.sum()*vxVolume],			# Volume in um3
			'radius_um' : [nRadius],						# Radius in um
			'dna_sum' : [np.nansum(idna[mask != 0])],		# DNA intensity integral
			'dna_mean' : [np.nanmean(idna[mask != 0])],		# DNA intensity mean
			'sig_sum' : [np.nansum(isig[mask != 0])],		# Signal intensity integral
			'sig_mean' : [np.nanmean(isig[mask != 0])],		# Signal intensity mean
			'surface_um2' : [nSurface],						# Surface in um2
			'sphericity' : 4*np.pi*(nRadius**2)/nSurface	# Sphericity
		}

		cData.append(pd.DataFrame.from_dict(nData))

	data.append(pd.concat(cData))

outData = pd.concat(data).reset_index(drop = True)

outData.to_csv(os.path.join(args.outdir, f"{args.prefix}.nuclear_features.tsv"),
	sep = "\t",index = False, header = True)

# END ==========================================================================

################################################################################
