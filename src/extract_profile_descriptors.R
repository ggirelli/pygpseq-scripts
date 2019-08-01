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
# parser = arg_parser(paste0(
# "Calculate "), name = scriptName)

# parser = add_argument(parser, arg = 'inpath',
# 	help = 'Path to nuclei table.')

# parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
# 	help = 'Number of threads for parallelization.', default = 1, nargs = 1)
# parser = add_argument(parser, arg = '--idcols', type = class(""),
# 	help = 'Columns for profile identification. Default: eid, sn.', nargs = Inf)
# parser = add_argument(parser, arg = '--profcol', type = class(""),
# 	help = 'Column with profile value. Default: sig_median.',
# 	default = "sig_median", nargs = 1)
# parser = add_argument(parser, arg = '--distcol', type = class(""),
# 	help = 'Column with distance from lamina value. Default: mid.',
# 	default = "mid", nargs = 1)
# parser = add_argument(parser, arg = '--xthr', short = '-x', type = class(0),
# 	help = 'Set threshold on the profile x axis for fitting and plotting.',
# 	default = c(0, 1), nargs = 2)

# version_flag = "0.0.1"
# parser = add_argument(parser, arg = '--version', short = '-V',
# 	help = 'Print version and quit.', flag = T)
# args = commandArgs(trailingOnly=TRUE)
# if ( "--version" %in% args ) {
# 	cat(sprintf("%s v%s\n", scriptName, version_flag))
# 	quit()
# }

# p = parse_args(parser)
# attach(p['' != names(p)], warn.conflicts = F)

# if ( all(is.na(idcols)) ) idcols = c("eid", "sn")

inpath = "/mnt/data/RADIANT/YFISH/iTK295_303/out3d_allNuclei/extra/3d_trajectories/nuclear_radii.tsv"
threads = 20
idcols = c("condition", "n", "line")
profcol = "signal"
distcol = "norm_lamina_dist"
xthr = c(0,1)

# FUNCTION =====================================================================

label_neighbours = function(x, xs, f) {
	xDiff = xs-x
	if ( !any(xDiff == 0) ) {
		bIDs = c(last(which(xDiff < 0)), which(xDiff > 0)[1])
		if ( 1 == length(bIDs) ) return(NA)
		if ( 1 != abs(diff(bIDs)) )
			bIDs = c(which(xDiff < 0)[1], last(which(xDiff > 0)))
	} else {
		match = which(xDiff == 0)
		bIDs = c(max(0, match-1), min(length(xs), match+1))
	}
	if ( 1 == length(bIDs) ) return(NA)
	bYs = f(xs[bIDs])
	label = rep(0, length(bYs))
	label[bYs > 0] = "+"
	label[bYs < 0] = "-"
	return(paste(label, collapse = ""))
}

extract_descriptors = function(nd, sigCol = "sig_median", distCol = "mid") {
	nd = nd[get(distCol) <= xthr[2] & get(distCol) >= xthr[1]]
	#nd[, c(sigCol) := .(as.numeric(gsub("\\[(.*)\\]", "\\1", get(sigCol))))]
	#nd[, c(distCol) := .(as.numeric(gsub("\\[(.*)\\]", "\\1", get(distCol))))]
	regrC = coefficients(lm(unlist(unlist(nd[, ..sigCol])) ~ poly(unlist(nd[, ..distCol]), 5, raw = T)))
	f <- eval(parse(text = paste0("function(x) (", regrC[1], "+x*", regrC[2],
			"+(x**2)*", regrC[3], "+(x**3)*", regrC[4], "+(x**4)*", regrC[5],
			"+(x**5)*", regrC[6], ")")))
	g <- function(x) {}
	body(g) <- D(body(f), 'x')
	h <- function(x) {}
	body(h) <- D(body(g), 'x')

	# layout(1:3)
	# xx = seq(0, 1, by = .001)
	# plot(xx, f(xx), type = "l")
	# plot(xx, g(xx), type = "l", col = 2); abline(h = 0, lty = 2)
	# plot(xx, h(xx), type = "l", col = 2); abline(h = 0, lty = 2)

	groots = uniroot.all(g, c(0,1))
	gdata = NULL
	if ( 0 != length(groots) ) {
		gdata = data.table(x = groots, y = f(groots),
			type = "1st Peak", col = sigCol)
		gdata$label = unlist(lapply(gdata$x, label_neighbours, nd[, get(distCol)], g))
		gdata = gdata[label == "+-"]
	}

	hroots = uniroot.all(h, c(0,1))
	hdata = NULL
	if ( 0 != length(hroots) ) {
		hdata = data.table(x = hroots, y = f(hroots),
			type = "1st Inflection", col = sigCol)
		hdata$label = unlist(lapply(hdata$x, label_neighbours, nd[, get(distCol)], h))
		hdata = hdata[label == "-+"]
	}

	rData = rbindlist(list(gdata[1, ], hdata[1, ]))
	if ( 0 == length(rData) ) return(NULL)

	rData = rbind(rData, data.table(x = NA,
		y = median(unlist(nd[, ..sigCol]), na.rm = T),
		type = "Median intensity",
		col = sigCol, label = NA))
	for (colName in idcols) {
		rData[, c(colName) := .(nd[1, get(colName)])]
	}

	return(rData)
}

# RUN ==========================================================================

cat("Reading input...\n")
nProfiles = fread(inpath, showProgress = T)
groupIDs = apply(nProfiles[, ..idcols], 1, paste, collapse = "-")
nProfiles = split(nProfiles, groupIDs)
remove("groupIDs")

cat("Going through profiles...\n")
pboptions(type = "timer")
pointData = rbindlist(pblapply(nProfiles, extract_descriptors,
	sigCol = profcol, distCol = distcol, cl = threads))

cat("Writing output...\n")
write.table(pointData, file.path(dirname(inpath), sprintf("%s.%s.points.tsv",
		tools::file_path_sans_ext(basename(inpath)), profcol)),
	quote = F, row.names = F, col.names = T, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
