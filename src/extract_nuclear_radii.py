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
import math
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import random
import re
import sys

from ggc.args import check_threads
from joblib import Parallel, delayed
from skimage.measure import marching_cubes_lewiner
from scipy.ndimage.measurements import center_of_mass
from scipy.interpolate import RegularGridInterpolator
from skimage.io import imread
from tqdm import tqdm


# PARAMETERS ===================================================================

parser = argparse.ArgumentParser(description = '''
Extract single nucleus single radii profiles from nuclear boxes, generated using
the -u flag of the gpseq_anim script. Use also --nuclear-sel with no parameters
to extract all nuclear boxes (useful if new selection should be performed).
''', formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('prefix', type = str, help = '''
	User prefix (usually in the format "iFL").''')
parser.add_argument('rootdir', type = str, help = '''
	Path to root directory with nuclear boxes.''')
parser.add_argument('outdir', type = str, help = '''
	Path to output directory where the output should be written to.''')

parser.add_argument('-n', '--nradii', metavar = 'nradii', type = int,
	default = 200, help = """Number of radii to extract. Default: 200""")
parser.add_argument('-b', '--npoints', metavar = 'npoints', type = int,
	default = 100, help = """Number of points to sample on the radii. Default: 100""")
parser.add_argument('-N', '--maxnuclei', metavar = 'maxNuclei', type = int,
	default = 500, help = """Maximum number nuclei to extract radii from. Default: 500""")
parser.add_argument('-S', '--seed', metavar = 'nthreads', type = int,
	help = """Seed for pseudorandom number generator.""")

parser.add_argument('--selected', type = str, help = """
	Path to table of selected nuclei. Mandatory columns: condition, sid, nid""")

parser.add_argument('-a', '--aspect', type = float, nargs = 3,
	help = """Physical size of Z, Y and X voxel sides in nm. Default: 300.0 216.6 216.6""",
	metavar = ('Z', 'Y', 'X'), default = [300., 216.6, 216.6])
parser.add_argument('-t', '--threads', metavar = 'nthreads', type = int,
	default = 1, help = """Number of threads to be used for parallelization.""")
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

if args.only_from_mid:
	sys.exit("2D analysis not yet implemented.")

print(f'''
 # {sys.argv[0]} v{version}
  
     Prefix : {args.prefix}
       Root : {args.rootdir}
   Selected : {args.selected}
     Output : {args.outdir}
     Aspect : {args.aspect}
    # Radii : {args.nradii}
   # Points : {args.npoints}
   # Nuclei : {args.maxnuclei}

  # threads : {args.threads}
       Seed : {args.seed}

 Midsection : {args.mid_type}
2D analysis : {args.only_from_mid}
''')

# FUNCTIONS ====================================================================

def extract_nuclear_radii(n, condition):
	sample = "%s%s" % (n, condition)

	# Read images
	mask = gp.tools.image.add_top_bottom_slides(
		imread("%s/%s.tif" % (args.rootdir, sample)))
	idna = gp.tools.image.add_top_bottom_slides(
		imread("%s/%s.dna.tif" % (args.rootdir, sample)))
	isig = gp.tools.image.add_top_bottom_slides(
		imread("%s/%s.sig.tif" % (args.rootdir, sample)))

	# Identify nuclear center of mass
	ncom = np.round(np.array(center_of_mass(mask)), 4)

	# Spread points homogeneously on a sphere
	points = np.round(np.array(
		fibonacci_sphere(samples = args.nradii, seed = args.seed)), 4)
	points += ncom

	# Identify 3D contour vertices (triangular mesh)
	vertsReal, faces, ns, vs = marching_cubes_lewiner(mask, 0.0, args.aspect)

	# Convert nm to vx
	vertsVX = vertsReal / args.aspect

	# Build lines from nuclear center of mass through dots on sphere surface
	lines = []
	for i in range(points.shape[0]):
		if any(0 == points[i] - ncom):
			continue
		my, qy, mz, qz = get3DLineThroughTwoPoints_byX(
			ncom, points[i])
		line3d_byX = build3Dline_byX(my, qy, mz, qz)
		lines.append({
			"coeff" : (my, qy, mz, qz),	# By X
			"fun" : line3d_byX,			# By X
			"p1" : ncom,
			"p2" : points[i]
		})

	# For each line, identify the intersection point with the 3D contour
	intersectionPointList = []
	for lineID in range(len(lines)):
		line = lines[lineID]
		# Find the 3D contour triangle with the closest baricenter
		faceDistList = []
		baricenterList = []
		for fi in range(faces.shape[0]):
			face = faces[fi, :]
			# Skip if not a triangle (e.g., overlapping points)
			if 3 != len(np.unique(vertsVX[face], axis = 0)):
				continue
			baricenter = np.mean(vertsVX[face], axis = 0).tolist()
			baricenter.append(fi)
			baricenterList.append(np.array(baricenter))
		baricenterTable = np.vstack(baricenterList)
		for bi in range(baricenterTable.shape[0]):
			faceDistList.append([getDistBetweenPointAnd3DLineThroughTwoPoints(
				baricenterTable[bi, 0:3], line['p1'], line['p2']), baricenterTable[bi, 3]])
		faceDistList = np.array(faceDistList)
		# Find 2 triangles that are the closest to a line
		faceIDs = faceDistList[np.argsort(
			faceDistList[:, 0])[:2], 1].astype('i')
		triangles = faces[faceIDs]
		# Find triangle closest to a point on the sphere
		# (used to draw the line)
		c1 = np.mean(vertsVX[triangles[0]], 0)
		c2 = np.mean(vertsVX[triangles[1]], 0)
		distTriangle1 = getEuDistBetweenTwoPoints(c1, line['p2'])
		distTriangle2 = getEuDistBetweenTwoPoints(c2, line['p2'])
		if distTriangle1 <= distTriangle2:
			triangle = vertsVX[triangles[0]]
		else:
			triangle = vertsVX[triangles[1]]
		# Find intersection point
		intersectionPoint = getPlaneLineIntersection(
			line['p1'], line['p2'],
			triangle[0, :], triangle[1, :], triangle[2, :])
		intersectionPointList.append(intersectionPoint)

	# Sample space from ncom to each intersection point
	sampledPointList = []
	for i in range(len(intersectionPointList)):
		intersectionPoint = intersectionPointList[i]
		line3D = build3Dline_vector(ncom, intersectionPoint)
		ds = np.linspace(0, 1, args.npoints)
		for d in ds:
			coords = line3D(d).tolist()
			coords.extend(intersectionPoint.tolist())
			coords.append(i)
			sampledPointList.append(np.array(coords))

	sampledPoints = np.array(sampledPointList)

	# For each point, interpolate the signal and calculate lamina distance --
	# Assemble all data into a table, for plotting purposes
	gridShape = np.array(idna.shape)
	regularGrid = [np.linspace(0, gridShape[i]-1, gridShape[i])
		for i in range(len(gridShape))]

	dna_gradient = RegularGridInterpolator(regularGrid, idna,
		method = "linear", bounds_error = False, fill_value = 0)
	sig_gradient = RegularGridInterpolator(regularGrid, isig,
		method = "linear", bounds_error = False, fill_value = 0)

	data = []
	for i in range(sampledPoints.shape[0]):
		pCoords = sampledPoints[i, :3].tolist()
		pSurf = sampledPoints[i, 3:6].tolist()
		lineID = sampledPoints[i, 6]
		dnaValue = (dna_gradient(pCoords))
		sigValue = (sig_gradient(pCoords))
		dCenter = getDistanceBetweenTwoPoints(ncom, pCoords)
		dLamina = getDistanceBetweenTwoPoints(pSurf, pCoords)
		row = pCoords
		row.extend([lineID, dnaValue, sigValue, dCenter, dLamina])
		data.append(row)

	data = pd.DataFrame(np.array(data))
	data.columns = ["z", "y", "x", "line",
		"dna", "signal", "centr_dist", "lamina_dist"]

	data['ratio'] = data['signal'] / data['dna']
	data['norm_lamina_dist'] = (data['centr_dist'] + data['lamina_dist'])
	data['norm_lamina_dist'] = data['lamina_dist'] / data['norm_lamina_dist']

	data['n'] = n
	data['condition'] = condition

	return data

# General ----------------------------------------------------------------------

def vec2scal(v):
	return(np.sqrt(np.sum(v**2)))

def getEuDistBetweenTwoPoints(p1, p2):
	return(np.sqrt(np.sum(np.power(p1 - p2, 2))))

# 2D ---------------------------------------------------------------------------

def getLineThroughTwoPoints(x0, y0, x1, y1):
	'''Calculate slope and intercept of a straight line cutting through
	two points.'''
	slope = np.round((y1-y0)/(x1-x0), 6)
	b = np.round(y0 - x0*(y1-y0)/(x1-x0), 6)
	return((slope, b))

def getLineThroughPointWithAngle(x0, y0, angle, isRadiants = False):
	'''Calculate slope intercept of a straight line cutting through a point
	with given angle.'''
	angle = np.deg2rad(angle) if not isRadiants else angle
	x1 = x0 + np.cos(angle)
	y1 = y0 + np.sin(angle)
	return(getLineThroughTwoPoints(x0, y0, x1, y1))

def getAngleOfLineFromSlope(m, inRadiants = False):
	'''Calculate angle from line slope. Carefule about the sign.'''
	angle = np.arctan(m)
	angle = np.deg2rad(angle) if inRadiants else angle
	return(angle)

def getAngleOfLineThroughTwoPoints(x0, y0, x1, y1,
	inRadiants = False, correctTo360 = True):
	'''Calculate angle of line cutting through two points.
	To use the first point as center, and the X axis as 0, use correctTo360.'''
	angle = getAngleOfLineFromSlope((y1-y0)/(x1-x0), inRadiants)
	angle = np.rad2deg(angle)
	if correctTo360:
		if x1 - x0 >= 0 :
			if y1 - y0 <=0:
				angle += 360
		else:
			if y1 - y0 >=0:
				angle = 180 + angle
			else:
				angle = 180 + angle
		#angle += 90 if y1 - y0 >=0 else 270
		#angle = (angle + 180) % 360
	angle = np.deg2rad(angle) if inRadiants else angle
	return(angle)

def getClosestPoints(com, x, y, pLower = None, pHigher = None, shift = 0):
	'''Identify the closest points to a line drawn through a point (com).'''
	if type(None) == type(pLower):
		pLower = (np.nan, np.nan, -np.inf)
	if type(None) == type(pHigher):
		pHigher = (np.nan, np.nan, np.inf)
	contourAngle = getAngleOfLineThroughTwoPoints(com[0], com[1], x, y)
	if contourAngle >= (angle + shift):
		if contourAngle <= pHigher[2]:
			pHigher = (x, y, contourAngle)
	else:
		if contourAngle >= pLower[2]:
			pLower = (x, y, contourAngle)
	return((pLower, pHigher))

def getIntersectionOfTwoLines(slope1, b1, slope2, b2):
	'''Find intersection of two straight lines.'''
	x = (b2 - b1) / (slope1 - slope2)
	y = slope1 * x + b1
	return((x, y))

def sampleLineBetweenTwoPoints(x0, y0, x1, y1, n = 200):
	'''Get coordinates of homogeneously spread points on a line, drawn between
	two points.'''
	if x0 == x1 or y0 == y1:
		xs = np.linspace(x0, x1, n)
		ys = np.linspace(y0, y1, n)
	else:
		slope, b = getLineThroughTwoPoints(x0, y0, x1, y1)
		xs = np.linspace(x0, x1, n)
		ys = np.array([slope*x+b for x in xs])
	return(np.vstack([xs, ys]).T)

def abline(slope, intercept, xShift = 0, yShift = 0):
	"""Plot a line from slope and intercept"""
	# From https://stackoverflow.com/a/43811762/1593536
	axes = plt.gca()
	x_vals = np.array(axes.get_xlim())
	y_vals = intercept + slope * x_vals
	plt.plot(x_vals+xShift, y_vals+yShift, '--')

def getDistanceBetweenTwoPoints(p1, p2):
	'''Euclidean distance between two points.'''
	p1 = np.array(p1)
	p2 = np.array(p2)
	return((np.sum((p1 - p2)**2))**.5)

# 3D ---------------------------------------------------------------------------

def fibonacci_sphere(samples = 1, randomize = True, seed = None):
	# From https://stackoverflow.com/a/26127012/1593536
	random.seed(seed)
	rnd = 1.
	if randomize:
		rnd = random.random() * samples
	points = []
	offset = 2./samples
	increment = math.pi * (3. - math.sqrt(5.));
	for i in range(samples):
		y = ((i * offset) - 1) + (offset / 2);
		r = math.sqrt(1 - pow(y,2))
		phi = ((i + rnd) % samples) * increment
		x = math.cos(phi) * r
		z = math.sin(phi) * r
		points.append([x,y,z])
	return points

def get3DLineThroughTwoPoints_byX(p1, p2):
	# Z: slice
	# Y: col
	# X: row
	z1, y1, x1 = tuple(p1.tolist())
	z2, y2, x2 = tuple(p2.tolist())
	my = (y2 - y1) / (x2 - x1)
	qy = y1 - x1 * (y2 - y1) / (x2 - x1)
	mz = (z2 - z1) / (x2 - x1)
	qz = z1 - x1 * (z2 - z1) / (x2 - x1)
	return((my, qy, mz, qz))

def build3Dline_byX(my, qy, mz, qz):
	return(eval("lambda x: np.round(np.array((x*%f+%f, x*%f+%f, x)), 6)" % (
		mz, qz, my, qy)))

def getDistBetweenPointAnd3DLineThroughTwoPoints(p0, p1, p2):
	# distance between p0 and the line through p1 and p2
	# With cross product
	return vec2scal(np.cross((p0 - p1), (p0 - p2))) / vec2scal(p2 - p1)

def findTriangleCutByLine(line, verts, faces):
	dists = []
	for fi in range(faces.shape[0]):
		f = faces[fi, :]
		if 3 != len(np.unique(verts2[f], axis = 0)):
			continue
		fCenter = np.mean(verts2[f], axis = 0)
		faceDist = getDistBetweenPointAnd3DLineThroughTwoPoints(fCenter, line['p1'], line['p2'])
		dists.append(np.array([faceDist, fi]))
	dists = np.array(dists)
	# Find 2 triangles that are the closest to a line
	triangles = faces[dists[np.argsort(dists[:, 0])[:2], 1].astype('i')]
	# Find triangle closest to a point on the sphere (used to draw the line)
	c1 = np.mean(verts2[triangles[0]], 0)
	c2 = np.mean(verts2[triangles[1]], 0)
	if getEuDistBetweenTwoPoints(c1, line['p2']) <= getEuDistBetweenTwoPoints(c2, line['p2']):
		triangle = verts2[triangles[0]]
	else:
		triangle = verts2[triangles[1]]
	return(triangle)

def getPlaneLineIntersection(P0, P1, A, B, C):
	# A, B, C are points that identify the plane
	# P0 and P1 are points that identify the line
	n = np.cross(A - C, B - C)
	l = P1 - P0
	d = np.dot(A - P0, n) / np.dot(l, n)
	intersectionPoint = d * l + P0
	return(intersectionPoint)

def build3Dline_vector(p1, p2):
	# Z: slice
	# Y: col
	# X: row
	z1, y1, x1 = tuple(p1.tolist())
	z2, y2, x2 = tuple(p2.tolist())
	funString = "lambda d: d*(np.array([%f,%f,%f])-np.array([%f,%f,%f]))" % (
		z2, y2, x2, z1, y1, x1) + "+np.array([%f,%f,%f])" % (z1, y1, x1)
	return(eval(funString))

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

selectedNuclei = {}
if type(None) != type(args.selected):
	assert os.path.isfile(args.selected)
	nTable = pd.read_csv(args.selected, sep = "\t")
	reqCols = ("condition", "sid", "nid")

	assert all([x in nTable.columns for x in reqCols])
	for i in range(nTable.shape[0]):
		n = nTable.loc[i]
		signature = f's{n["sid"]}n{n["nid"]}'
		if n['condition'] not in selectedNuclei.keys():
			selectedNuclei[n['condition']] = set()
		selectedNuclei[n['condition']].add(signature)

	assert 0 != np.sum([len(nd[x]) for x in nd.keys()])
	for eid in nd.keys():
		nd[eid] = [n for n in nd[eid] if n in selectedNuclei[eid]]

assert 0 != np.sum([len(nd[x]) for x in nd.keys()])

data = []
cid = 0
for condition in sorted(nd.keys()):
	cdata = Parallel(n_jobs = args.threads, verbose = 0)(
	    delayed(extract_nuclear_radii)(n, condition)
	    for n in tqdm(nd[condition][:args.maxnuclei],
	    	desc = f'{cid+1}/{len(nd)} {condition} [n.threads={args.threads}]'))
	print(condition)
	pd.concat(cdata).to_csv(os.path.join(args.outdir,
		f"{condition}.nuclear_radii.tsv"), index = False, sep = "\t")
	data.extend(cdata)
	cid += 1
pd.concat(data).to_csv(os.path.join(args.outdir,
	f"nuclear_radii.tsv"), index = False, sep = "\t")

# END ==========================================================================

################################################################################
