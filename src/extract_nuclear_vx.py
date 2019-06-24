#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20181114
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

from ggc.args import check_threads
from joblib import Parallel, delayed
from scipy.ndimage.measurements import center_of_mass
from skimage.io import imread
from tqdm import tqdm

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
Extract single nucleus voxel data from nuclear boxes, generating using the -u flag
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
parser.add_argument('-t', '--threads', metavar = 'nthreads', type = int,
	default = 1, help = """Number of threads to be used for parallelization.""")

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

if type(None) == type(args.prefix):
	args.prefix = os.basename(args.outdir)

args.threads = check_threads(args.threads)

# FUNCTIONS ====================================================================

def buildVxNuclearProfile(n, condition):
	import commonFunctions as cF
	from scipy.interpolate import RegularGridInterpolator

	sample = "%s%s" % (n, condition)

	# Read images
	mask = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.tif'))
	idna = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.dna.tif'))
	isig = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.sig.tif'))

	# Identify nuclear center of mass
	ncom = np.array(center_of_mass(mask)).astype("i")

	empty = np.zeros(mask.shape)
	empty[ncom[0], ncom[1], ncom[2]] = 1
	empty -= 1
	empty *= -1

	laminD = gp.tools.distance.calc_lamina_distance(mask, aspect)
	centrD = gp.tools.distance.calc_lamina_distance(empty, aspect)

	mask_flat = mask.reshape([np.prod(mask.shape)])
	mask_flat = np.where(mask_flat.tolist())[0].tolist()

	# Prepare output
	DTYPE_NUCLEAR_DATA = [
		('dna', 'u8'),
		('signal', 'u8'),
		('ratio', 'f8'),
		("lamina_dist", 'f8'),
		("centr_dist", 'f8'),
		("norm_lamina_dist", 'f8'),
	]
	data = np.zeros(len(mask_flat), dtype = DTYPE_NUCLEAR_DATA)

	# Flatten data for export
	data['dna'] = gp.tools.vector.flatten_and_select(idna, mask_flat)
	data['signal'] = gp.tools.vector.flatten_and_select(isig, mask_flat)
	data['lamina_dist'] = gp.tools.vector.flatten_and_select(laminD, mask_flat)
	data['centr_dist'] = gp.tools.vector.flatten_and_select(centrD, mask_flat)

	data['norm_lamina_dist'] = gp.tools.vector.flatten_and_select(
		laminD / (laminD + centrD), mask_flat)
	data['ratio'] = gp.tools.vector.flatten_and_select(
		isig / idna, mask_flat)

	df = pd.DataFrame(data)
	df['c'] = condition
	df['sn'] = n
	df.to_csv(os.path.join(args.outdir, f"{sample}.vx.tsv"),
		sep = "\t", index = False)

# RUN ==========================================================================

# Identify nuclei
sampleList = list(set([x.split(".")[0] for x in os.listdir(args.rootdir)]))

regexp = re.compile(r"^(.*)(%s.*)$" % prefix)
nd = {}
for n in tqdm(sampleList, desc = "Dividing in conditions"):
	m = re.match(regexp, n).groups()
	if m[1] in nd.keys():
		nd[m[1]].append(m[0])
	else:
		nd[m[1]] = [m[0]]

for condition in nd.keys():
	Parallel(n_jobs = ncores, verbose = 11)(
	    delayed(buildVxNuclearProfile)(n, condition)
	    for n in nd[condition])

# END ==========================================================================

################################################################################
