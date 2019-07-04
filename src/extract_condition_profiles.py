#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190620
# 
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import sys

from ggc.args import check_threads
from joblib import Parallel, delayed
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
parser.add_argument('-n', '--nbins', type = int, help = """
	Number of bins from lamina to center. Default: 200""", default = 200)
parser.add_argument('--selected', type = str, help = """
	Path to table of selected nuclei. Mandatory columns: condition, sid, nid""")
parser.add_argument('-S', '--suffix', type = str, help = """
	Suffix for output files.""", default = "")
parser.add_argument('-O', '--outdir', type = str, help = """
	Path to output directory where the output should be written to.""", default = "")
parser.add_argument('-t', '--threads', metavar = 'nthreads', type = int,
	default = 1, help = """Number of threads to be used for parallelization.""")

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

assert os.path.isdir(args.rootdir)
if type(None) == type(args.outdir):
	args.outdir = os.path.dirname(args.rootdir)
else:
	os.path.dirname(args.outdir)

if 0 != len(args.suffix):
	if not args.suffix.startswith("."):
		args.suffix = f".{args.suffix }"

args.threads = check_threads(args.threads)

print(f'''
 # {sys.argv[0]} v{version}

   Prefix : {args.prefix}
     Root : {args.rootdir}
 Selected : {args.selected}
   Suffix : "{args.suffix}"
   Output : {args.outdir}
    #bins : {args.nbins}
 #threads : {args.threads}
''')

# FUNCTIONS ====================================================================

def mkStatProfile(data, k):
	d = {}
	if 0 != len(data):
		d[k+'_mean'] = [np.nanmean(data)]
		d[k+'_median'] = [np.nanmedian(data)]
	else:
		d[k+'_mean'] = [np.nan]
		d[k+'_median'] = [np.nan]
	return(d)

def get_nucleus_bins(fname, breaks, args):
	bins = [{"mid":np.mean(breaks[i:(i+1)]),"dna":[],"sig":[],"rat":[]}
		for i in range(args.nbins)]

	fname = f"{fname}.vx.tsv"
	with open(os.path.join(args.rootdir, fname), "r") as IH:
		drop = next(IH)
		for line in IH:
			data = line.strip().split("\t")

			x = float(data[5])
			for bid in range(args.nbins):
				if breaks[bid] >= x:
					break

			bins[bid]['dna'].append(float(data[0]))
			bins[bid]['sig'].append(float(data[1]))
			bins[bid]['rat'].append(float(data[2]))

	return(bins)

# RUN ==========================================================================

flist = os.listdir(args.rootdir)

meta = {}
for fname in flist:
	if fname.endswith(".vx.tsv"):
		eid = args.prefix + fname.split("_")[0].split(args.prefix)[1]
		if eid not in meta.keys():
			meta[eid] = [fname.split(".")[0]]
		else:
			meta[eid].append(fname.split(".")[0])

selectedNuclei = set()
if type(None) != type(args.selected):
	assert os.path.isfile(args.selected)
	nTable = pd.read_csv(args.selected, sep = "\t")
	reqCols = ("condition", "sid", "nid")

	assert all([x in nTable.columns for x in reqCols])
	for i in range(nTable.shape[0]):
		n = nTable.loc[i]
		signature = f's{n["sid"]}n{n["nid"]}{n["condition"]}'
		selectedNuclei.add(signature)

	assert 0 != np.sum([len(meta[x]) for x in meta.keys()])
	for eid in meta.keys():
		meta[eid] = [n for n in meta[eid] if n in selectedNuclei]

assert 0 != np.sum([len(meta[x]) for x in meta.keys()])

breaks = np.linspace(0, 1, args.nbins+1)
allData = []
eidn = 0
for eid in sorted(meta.keys()):
	bins = [{"mid":np.mean(breaks[i:(i+1)]),"dna":[],"sig":[],"rat":[]}
		for i in range(args.nbins)]

	binList = Parallel(n_jobs = args.threads, verbose = 0)(
	    delayed(get_nucleus_bins)(fname, breaks, args)
	    for fname in tqdm(meta[eid],
	    	desc = f'{eidn+1}/{len(meta)} {eid} [n.threads={args.threads}]'))

	while 0 < len(binList):
		binSet = binList.pop()
		bid = 0
		while 0 < len(binSet):
			single_bin = binSet.pop()
			bins[bid]['dna'].extend(single_bin['dna'])
			bins[bid]['sig'].extend(single_bin['sig'])
			bins[bid]['rat'].extend(single_bin['rat'])
			bid += 1

	for bid in range(len(bins)):
		binData = {"mid":[bins[bid]['mid']],"eid":eid}
		for k in ['dna', 'sig', 'rat']:
			binData.update(mkStatProfile(bins[bid][k], k))
		bins[bid] = pd.DataFrame.from_dict(binData)

	bins = pd.concat(bins).reset_index(drop = True)
	bins.to_csv(os.path.join(args.outdir,
			f'{eid}.condition.profiles{args.suffix}.tsv'),
		sep = '\t', index = False)
	allData.append(bins)
	eidn += 1

print("Merging and writing...")
pd.concat(allData).reset_index(drop = True).to_csv(
	os.path.join(args.outdir, f'condition.profiles{args.suffix}.tsv'),
	sep = '\t', index = False, header = True)

# END ==========================================================================

################################################################################
