#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190703
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(ggrepel))
suppressMessages(require(LaCroixColoR))
suppressMessages(require(rootSolve))
suppressMessages(require(RColorBrewer))
suppressMessages(require(viridis))

theme_set(theme_cowplot())

# PARAMS =======================================================================

scriptName = "plot_condition_profiles.R"
parser = arg_parser(paste0(
"Generate plots on condition profiles, maxima peaks, inflection points, and 
relative AUC. Requires extracted condition profiles. The extra folder should
contain the following subfolders: nuclear_selection, mid_profiles, 3d_profiles;
and the following files: meta.tsv.
"), name = scriptName)

parser = add_argument(parser, arg = 'dataset',
	help = 'Dataset ID, e.g., iXX1_2')
parser = add_argument(parser, arg = 'rootDir',
	help = 'Path to extra folder.')

parser = add_argument(parser, arg = '--timepoints', short = '-t',
	type = class(""), nargs = 1,
	help = 'Path to timepoint meta tsv table with levels and time (min) columns.')
parser = add_argument(parser, arg = '--xthr', short = '-x', type = class(0),
	help = 'Set threshold on the profile x axis for fitting and plotting.',
	default = c(0, 1), nargs = 2)
parser = add_argument(parser, arg = '--outpath', short = "-O", type = class(""),
	help = 'Path to output folder. Default to root dir.')

p = parse_args(parser)
attach(p['' != names(p)])

if ( is.na(timepoints) ) {
	cLevels = c(paste0(c(1, 2, 5, 10, 15, 30, 45), "min"),
		paste0(c(1, 2, 6), "h"), "ON", "LongON")
	cTime = c(1, 2, 5, 10, 15, 30, 45, c(1, 2, 6)*60, 16*60, 72*60)
	cMeta = data.table(
		levels = factor(cLevels, levels = cLevels),
		time = cTime, key = "levels"
	)
} else {
	stopifnot(file.exists(timepoints))
	cMeta = fread(timepoints, col.names = c("levels", "time"), key = "levels")
}
cMeta$colors = lacroix_palette(type = "paired", 8)[1:nrow(cMeta)]

if ( is.na(outpath) ) outpath = rootDir

# FUNCTIONS ====================================================================

labelNeighbours = function(x, xs, f) {
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

save_and_plot <- function(x, bname, width, height,
	dpi=300, use.cowplot=FALSE, ncol=1, nrow=1,
	base_height=4, base_width=NULL, base_aspect_ratio=1,
	plot=FALSE, formats = c("png", "pdf")){
  # Function to save the plot (and separately its data)
  # to file and show the plot in the notebook
  if( !use.cowplot ){
	if ( "png" %in% formats) {
		png(filename=file.path(paste0(bname, ".png")),
			units="in", height=height, width=width, res=dpi)
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "pdf" %in% formats) {
		cairo_pdf(filename=file.path(paste0(bname, "_cairo_pdf.pdf")),
			onefile = TRUE, height=height, width=width, family="Helvetica",
			pointsize=8, antialias="none")
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "eps" %in% formats) {
		cairo_ps(filename=file.path(paste0(bname, "_cairo_ps.eps")),
			onefile = TRUE, height=height, width=width, family="Helvetica",
			pointsize=8, antialias="none")
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
		postscript(file=file.path(paste0(bname, "_postscript.eps")),
			onefile = TRUE, paper="special", height=height, width=width,
			family="Helvetica", pointsize=8, horizontal=FALSE)
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
  }else{
	if ( "png" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, ".png")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			dpi=dpi)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "pdf" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, "_cairo_pdf.pdf")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device=cairo_pdf)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "eps" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, "_cairo_ps.eps")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device=cairo_ps)
		save_plot(x, filename=file.path(paste0(bname, "_postscript.eps")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device="ps")  
		while(length(dev.list())>0) invisible(dev.off())
	}
  }
  if( plot ) print(x)
}

# RUN ==========================================================================

meta = fread(file.path(rootDir, "meta.tsv"), key = "condition")
stopifnot(all(c("condition", "timeLab", "label") %in% names(meta)))
meta[, label := as.character(label)]
meta[is.na(label), label := ""]

nf = fread(file.path(rootDir,
	sprintf("nuclear_selection/%s.nuclear_features.tsv", dataset)))
nf$condition = unlist(lapply(nf$condition,
	function(x) unlist(strsplit(x, "_", fixed = T))[1]))
setkeyv(nf, c("condition", "sid", "nid"))
nf = nf[meta]
nf[, info := paste0(condition, "-", timeLab, "-", label)]

setkeyv(nf, "timeLab")
nf = nf[cMeta][!is.na(condition)]
nfMeta = nf[, .(color = unique(colors)), by = info]
nfColors = nfMeta$color
names(nfColors) = nfMeta$info

p = ggplot(melt(nf, id.vars = c("condition", "sid", "nid",
		"timeLab", "label", "info", "colors")),
		aes(x = value, color = info)
	) + geom_density(
	) + facet_wrap(~variable, scales = "free"
	) + guides(color = guide_legend(title = "Condition")
	) + theme(legend.position = "top",
		axis.text.x = element_text(angle = 45, hjust = 1)
	) + scale_color_manual(breaks = nfMeta$info, values = nfColors
	) + ggtitle(dataset
	)
save_and_plot(p, file.path(outpath,
	sprintf("%s.nuclear_features.density", dataset)),
	format = "png", 15, 9)

l = lapply(c("mid", "3d"), function(atype) {
	pr = fread(file.path(rootDir,
		sprintf("%s_profiles/condition.profiles.%s.%s.tsv",
		atype, dataset, atype)))

	pr$eid = unlist(lapply(pr$eid,
		function(x) unlist(strsplit(x, ".", fixed = T))[1]))
	setkeyv(pr, "eid")
	pr = pr[meta]

	pr[, timeLab := factor(timeLab, levels = cMeta$levels)]
	pr[, info := reorder(
		paste0(eid, "-", timeLab, "-", label), as.numeric(timeLab))]
	setkeyv(pr, "timeLab")
	pr = pr[cMeta][!is.na(eid)]
	prMeta = pr[, .(color = unique(colors)), by = c("info", "label")]
	prColors = prMeta$color
	names(prColors) = prMeta$info

	pr = pr[mid <= xthr[2] & mid >= xthr[1]]

	p = ggplot(pr, aes(x = mid, y = sig_median, color = info)
		) + geom_point() + geom_line(
		) + geom_smooth(method = "lm", formula = "y ~ poly(x, 5, raw = T)",
			fill = NA, linetype = "dashed", color = "#323232"
		) + xlab("Normalized lamina distance (a.u.)"
		) + ylab("Median Signal channel intensity (a.u.)"
		) + guides(color = F
		) + ggtitle(sprintf("%s - %s", dataset, atype)
		) + facet_wrap(~info
		) + scale_color_manual(breaks = prMeta$info, values = prColors
		) + xlim(0, 1
		)
	save_and_plot(p, file.path(outpath,
		sprintf("%s.%s.profiles.model", dataset, atype)),
		format = "png", 15, 9)

	p = ggplot(pr, aes(x = mid, y = sig_median, color = info)
		) + geom_point() + geom_line(
		) + geom_smooth(method = "lm", formula = "y ~ poly(x, 5, raw = T)",
			fill = NA, linetype = "dashed", color = "#323232"
		) + xlab("Normalized lamina distance (a.u.)"
		) + ylab("Median Signal channel intensity (a.u.)"
		) + guides(color = F
		) + ggtitle(sprintf("%s - %s", dataset, atype)
		) + facet_wrap(~info, scales = "free"
		) + scale_color_manual(breaks = prMeta$info, values = prColors
		) + xlim(0, 1
		)
	save_and_plot(p, file.path(outpath,
		sprintf("%s.%s.profiles.model.free", dataset, atype)),
		format = "png", 15, 9)

	p = ggplot(pr, aes(x = mid, y = sig_median, color = info)
		) + geom_smooth(method = "lm", formula = "y ~ poly(x, 5, raw = T)",
			fill = NA
		) + xlab("Normalized lamina distance (a.u.)"
		) + ylab("Median Signal channel intensity (a.u.)"
		) + guides(color = guide_legend(title = "Condition", nrow = 2)
		) + theme(legend.position = "top"
		) + ggtitle(sprintf("%s - %s", dataset, atype)
		) + scale_color_manual(breaks = prMeta$info, values = prColors
		) + xlim(0, 1
		)
	save_and_plot(p, file.path(outpath,
		sprintf("%s.%s.profiles.model.overlay", dataset, atype)),
		format = "png", 9, 9)

	pointData = rbindlist(by(pr, pr$info, function(ct) {
		regrC = coefficients(lm(ct$sig_median ~ poly(ct$mid, 5, raw = T)))
		f <- eval(parse(text = paste0("function(x) (",
			regrC[1], "+x*", regrC[2], "+(x**2)*", regrC[3], "+(x**3)*",
			regrC[4], "+(x**4)*", regrC[5], "+(x**5)*", regrC[6], ")")))
		g <- function(x) {}
		body(g) <- D(body(f), 'x')
		h <- function(x) {}
		body(h) <- D(body(g), 'x')

		groots = uniroot.all(g, c(0,1))
		gdata = NULL
		if ( 0 != length(groots) ) {
			gdata = data.table(x = groots, y = 0, v = "1st Peak", c = "signal")
			gdata$label = unlist(lapply(gdata$x, labelNeighbours, ct$mid, g))
			gdata = gdata[label == "+-"]
		}

		hroots = uniroot.all(h, c(0,1))
		hdata = NULL
		if ( 0 != length(hroots) ) {
			hdata = data.table(x = hroots, y = 0,
				v = "1st Inflection", c = "signal")
			hdata$label = unlist(lapply(hdata$x, labelNeighbours, ct$mid, h))
			hdata = hdata[label == "-+"]
		}

		rData = rbindlist(list(gdata[1, ], hdata[1, ]))
		rData$e = ct[1, eid]
		rData$l = ct[1, label]
		rData$t = ct[1, time]

		return(rData)
	}))
	write.table(pointData, file.path(outpath,
		sprintf("%s.%s.point_data.tsv", dataset, atype)),
		quote = F, row.names = F, col.names = T, sep = "\t")

	p = ggplot(pointData[!is.na(c)], aes(y = x, x = t, group = v, color = v)
		) + geom_point() + geom_line(
		) + ylab("Normalized distance from lamina"
		) + xlab("Time point") + ylim(0, 1
		) + guides(color = guide_legend("Point", title.position = "left")
		) + theme(legend.position = "top"
		) + scale_color_brewer(palette = "Set1"
		) + ggtitle(sprintf("%s - %s", dataset, atype)
		) + geom_label_repel(aes(label = l),
			segment.color = '#989898', segment.size = .25,
			color = "#323232", size = 3,
			point.padding = .5, force = 20
		)
	save_and_plot(p, file.path(outpath,
		sprintf("%s.%s.profiles.points", dataset, atype)),
		format = "png", 12, 6)
})

warnings()

# END ==========================================================================

################################################################################
