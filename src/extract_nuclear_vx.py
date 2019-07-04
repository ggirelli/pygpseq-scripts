#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20181114
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

parser = argparse.ArgumentParser(description = '''
Extract single nucleus voxel data from nuclear boxes, generating using the -u flag
of the gpseq_anim script. Use also --nuclear-sel with no parameters to extract
all nuclear boxes (useful if new selection should be performed).
''', formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('prefix', type = str, help = '''
	User prefix (usually in the format "iFL").''')
parser.add_argument('rootdir', type = str, help = '''
	Path to root directory with nuclear boxes.''')
parser.add_argument('outdir', type = str, help = '''
	Path to output directory where the output should be written to.''')

parser.add_argument('-a', '--aspect', type = float, nargs = 3,
	help = """Physical size of Z, Y and X voxel sides in nm. Default: 300.0 216.6 216.6""",
	metavar = ('Z', 'Y', 'X'), default = [300., 216.6, 216.6])
parser.add_argument('-t', '--threads', metavar = 'nthreads', type = int,
	default = 1, help = """Number of threads to be used for parallelization.""")
parser.add_argument('--dist-type', type = str,
	help = """Method for lamina distance calculation/normalization.
	Default: '%s'""" % (gp.const.LD_ARG_LABELS[gp.const.LD_DEFAULT]),
	choices = list(gp.const.LD_ARG_LABELS),
	default = gp.const.LD_ARG_LABELS[gp.const.LD_DEFAULT])
parser.add_argument('--mid-type', type = str,
	help = """Method for mid-section selection. Default: '%s'""" % (
		gp.const.MID_SEC_ARG_LABELS[gp.const.MID_SEC_DEFAULT]),
	choices = list(gp.const.MID_SEC_ARG_LABELS),
	default = gp.const.MID_SEC_ARG_LABELS[gp.const.MID_SEC_DEFAULT])

parser.add_argument('-2', action = 'store_const', dest = 'only_from_mid',
	const = True, default = False, help = '''Extract only voxels from the nucleus mid section.
	Check the --mid-type option for how to change section selection method.''')

version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

args = parser.parse_args()

assert os.path.isdir(args.rootdir)
assert os.path.isdir(args.outdir)

if type(None) == type(args.prefix):
	args.prefix = os.basename(args.outdir)

args.threads = check_threads(args.threads)
args.mid_type = gp.const.MID_SEC_ARG_LABELS.index(args.mid_type)
args.dist_type = gp.const.LD_ARG_LABELS.index(args.dist_type)

print(f'''
 # {sys.argv[0]} v{version}
  
     Prefix : {args.prefix}
       Root : {args.rootdir}
     Output : {args.outdir}
     Aspect : {args.aspect}
   Distance : {args.dist_type}
 Midsection : {args.mid_type}
   #threads : {args.threads}
2D analysis : {args.only_from_mid}
''')

# FUNCTIONS ====================================================================

def buildVxNuclearProfile(n, condition):
	sample = "%s%s" % (n, condition)

	mask = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.tif'))
	idna = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.dna.tif'))
	isig = gp.tools.image.add_top_bottom_slides(imread(f'{args.rootdir}/{sample}.sig.tif'))

	if args.only_from_mid:
		largestID = gp.tools.image.get_mid_section_idx(idna, mask, args.mid_type)
		mask = mask[largestID]
		idna = idna[largestID]
		isig = isig[largestID]

	ncom = np.array(center_of_mass(mask)).astype("i")

	laminD, centrD = gp.tools.distance.calc_nuclear_distances(
		args.dist_type, mask, args.aspect)

	mask_flat = mask.reshape([np.prod(mask.shape)])
	mask_flat = np.where(mask_flat.tolist())[0].tolist()

	DTYPE_NUCLEAR_DATA = [
		('dna', 'u8'),
		('signal', 'u8'),
		('ratio', 'f8'),
		("lamina_dist", 'f8'),
		("centr_dist", 'f8'),
		("norm_lamina_dist", 'f8'),
	]
	data = np.zeros(len(mask_flat), dtype = DTYPE_NUCLEAR_DATA)

	data['dna'] = gp.tools.vector.flatten_and_select(idna, mask_flat)
	data['signal'] = gp.tools.vector.flatten_and_select(isig, mask_flat)
	data['lamina_dist'] = gp.tools.vector.flatten_and_select(laminD, mask_flat)
	data['centr_dist'] = gp.tools.vector.flatten_and_select(centrD, mask_flat)

	data['norm_lamina_dist'] = gp.tools.vector.flatten_and_select(
		laminD / (laminD + centrD), mask_flat)
	np.seterr(divide='ignore', invalid='ignore')
	data['ratio'] = gp.tools.vector.flatten_and_select(isig / idna, mask_flat)

	df = pd.DataFrame(data)
	df['c'] = condition
	df['sn'] = n
	df.to_csv(os.path.join(args.outdir, f"{sample}.vx.tsv"),
		sep = "\t", index = False)

# RUN ==========================================================================

sampleList = list(set([x.split(".")[0] for x in os.listdir(args.rootdir)]))
regexp = re.compile(r"^(.*)(%s.*)$" % args.prefix)

nd = {}
for n in tqdm(sampleList, desc = "Dividing in conditions"):
	m = re.match(regexp, n).groups()
	if m[1] in nd.keys():
		nd[m[1]].append(m[0])
	else:
		nd[m[1]] = [m[0]]

cid = 0
for condition in sorted(nd.keys()):
	Parallel(n_jobs = args.threads, verbose = 0)(
	    delayed(buildVxNuclearProfile)(n, condition)
	    for n in tqdm(nd[condition],
	    	desc = f'{cid+1}/{len(nd)} {condition} [n.threads={args.threads}]'))
	cid += 1

# END ==========================================================================

################################################################################
