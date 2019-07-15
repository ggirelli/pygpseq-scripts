#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190704
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(require(argparser))
suppressMessages(require(data.table))
suppressMessages(require(pbapply))
suppressMessages(require(rootSolve))

setDTthreads(1)

# INPUT ========================================================================

scriptName = "extract_profile_descriptors.R"
parser = arg_parser(paste0(
"Calculate "), name = scriptName)

parser = add_argument(parser, arg = 'inpath',
	help = 'Path to nuclei table.')

parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)
parser = add_argument(parser, arg = '--idcols', type = class(""),
	help = 'Columns for profile identification. Default: eid, sn.',
	default = c("eid", "sn"), nargs = Inf)
parser = add_argument(parser, arg = '--profcol', type = class(""),
	help = 'Column with profile value. Default: sig_median.',
	default = "sig_median", nargs = 1)

version_flag = "0.0.1"
parser = add_argument(parser, arg = '--version', short = '-V',
	help = 'Print version and quit.', flag = T)
args = commandArgs(trailingOnly=TRUE)
if ( "--version" %in% args ) {
	cat(sprintf("%s v%s\n", scriptName, version_flag))
	quit()
}

p = parse_args(parser)
attach(p['' != names(p)], warn.conflicts = F)

# FUNCTION =====================================================================

label_neighbours = function(x, xs, f) {
	xDiff = xs-x
	if ( !any(xDiff == 0) ) {
		bIDs = c(last(which(xDiff < 0)), which(xDiff > 0)[1])
	} else {
		match = which(xDiff == 0)
		bIDs = c(max(0, match-1), min(length(xs), match+1))
	}
	bYs = f(xs[bIDs])
	label = rep(0, length(bYs))
	label[bYs > 0] = "+"
	label[bYs < 0] = "-"
	return(paste(label, collapse = ""))
}

extract_descriptors = function(nd, sigCol = "sig_median") {
	regrC = coefficients(lm(unlist(nd[, ..sigCol]) ~ poly(nd$mid, 5, raw = T)))
	f <- eval(parse(text = paste0("function(x) (", regrC[1], "+x*", regrC[2],
			"+(x**2)*", regrC[3], "+(x**3)*", regrC[4], "+(x**4)*", regrC[5],
			"+(x**5)*", regrC[6], ")")))
	g <- function(x) {}
	body(g) <- D(body(f), 'x')
	h <- function(x) {}
	body(h) <- D(body(g), 'x')

	groots = uniroot.all(g, c(0,1))
	gdata = NULL
	if ( 0 != length(groots) ) {
		gdata = data.table(x = groots, y = 0,
			type = "1st Peak", col = sigCol)
		gdata$label = unlist(lapply(gdata$x, label_neighbours, nd$mid, g))
		gdata = gdata[label == "+-"]
	}

	hroots = uniroot.all(h, c(0,1))
	hdata = NULL
	if ( 0 != length(hroots) ) {
		hdata = data.table(x = hroots, y = 0,
			type = "1st Inflection", col = sigCol)
		hdata$label = unlist(lapply(hdata$x, label_neighbours, nd$mid, h))
		hdata = hdata[label == "-+"]
	}

	rData = rbindlist(list(gdata[1, ], hdata[1, ]))
	if ( 0 == length(rData) ) return(NULL)
	rData$eid = nd[1, eid]
	rData$sn = nd[1, sn]

	return(rData)
}

# RUN ==========================================================================

cat("Reading input...\n")
nProfiles = fread(inpath)
groupIDs = apply(nProfiles[, ..idcols], 1, paste, collapse = "-")

cat("Going through profiles...\n")
pboptions(type = "timer")
pointData = rbindlist(pblapply(split(nProfiles, groupIDs),
	extract_descriptors, sigCol = profcol, cl = threads))

cat("Writing output...\n")
write.table(pointData, file.path(dirname(inpath), sprintf("%s.%s.points.tsv",
		tools::file_path_sans_ext(basename(inpath)), profcol)),
	quote = F, row.names = F, col.names = T, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
