#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190620
# Project: GPSeq-YFISH
# 
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import sys
from tqdm import tqdm

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser( description = '''
Generate mean and median condition profiles, from lamina to center, starting
from previously extract nuclear voxel data.
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('prefix', type = str, help = '''
	User prefix (usually in the format "iFL").''')
parser.add_argument('rootdir', type = str, help = '''
	Path to root directory with nuclear voxel content.''')

# Add arguments with default value
parser.add_argument('-n', '--args.n', type = int,
	help = """Number of bins from lamina to center.""", default = 200)

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

# FUNCTIONS ====================================================================

def mkStatProfile(data, k):
	d = {}
	if 0 != len(xbin[k]):
		binData[k+'_mean'] = [np.nanmean(xbin[k])]
		binData[k+'_median'] = [np.nanmedian(xbin[k])]
	else:
		binData[k+'_mean'] = [np.nan]
		binData[k+'_median'] = [np.nan]
	return(d)

# RUN ==========================================================================

flist = os.listdir(args.rootdir)
meta = {}
for fname in flist:
	eid = args.prefix + fname.split("_")[0].split(args.prefix)[1]
	if eid not in meta.keys():
		meta[eid] = [fname]
	else:
		meta[eid].append(fname)

breaks = np.linspace(0, 1, args.n)
for eid in meta.keys():
	bins = [{"mid":np.mean(breaks[i:(i+1)]),"dna":[],"sig":[],"rat":[]}
		for i in range(args.n)]

	for fname in tqdm(meta[eid], desc = "Nucleus"):
		with open(os.path.join(args.rootdir, fname), "r") as IH:
			drop = next(IH)
			for line in tqdm(IH, desc = "Reading file"):
				data = line.strip().split("\t")
				bid = np.where((breaks<float(data[5])))[0][-1]
				bins[bid]['dna'].append(float(data[0]))
				bins[bid]['sig'].append(float(data[1]))
				bins[bid]['rat'].append(float(data[2]))

	data = []
	for xbin in bins:
		binData = {"mid":[xbin['mid']]}
		for k in ['dna', 'sig', 'rat']:
			binData.update(mkStatProfile(xbin[k], k))
		data.append(pd.DataFrame.from_dict(binData))

	data = pd.concat(data).reset_index(drop = True)
	data['eid'] = eid

	data.to_csv(os.path.join(os.path.dirname(args.rootdir),
			f'{eid}.condition.profiles.tsv'),
		sep = '\t', index = False)

# END ==========================================================================

################################################################################
