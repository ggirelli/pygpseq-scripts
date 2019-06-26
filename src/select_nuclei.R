#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190624
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

require(argparser)
require(cowplot)
require(data.table)
require(ggplot2)
require(parallel)
require(viridis)

# INPUT ========================================================================

# Create arguent parser
scriptName = "nucleiSelectPopulation.R"
parser = arg_parser(paste0(
"Subselect nuclei based on volume (vx) and DNA stain intensity integral.
Specifically, fits sum of Gaussians to the distribution, then select an interval
of +-k*sigma around the mean of the first Gaussian.
"), name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'inpath',
	help = 'Path to nuclei table.')

# Define elective arguments
parser = add_argument(parser, arg = '--ksigma', short = '-k', type = class(0),
	help = 'Sigma constant for interval definition.', default = 3, nargs = 1)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)
parser = add_argument(parser, arg = '--outpath', short = "-O", type = class(""),
	help = 'Path to output folder. Default to inpath parent folder.')

# Version argument
version_flag = "0.0.1"
parser = add_argument(parser, arg = '--version', short = '-V',
	help = 'Print version and quit.', flag = T)

args = commandArgs(trailingOnly=TRUE)
if ( "--version" %in% args ) {
	cat(sprintf("%s v%s\n", scriptName, version_flag))
	quit()
}

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)], warn.conflicts = F)

if ( !file.exists(inpath) )
	stop(sprintf("nuclei table not found: %s", inpath))
if ( is.na(outpath) ) outpath = dirname(inpath)
if ( !dir.exists(outpath) )
	stop(sprintf("output folder not found: %s", outpath))

# FUNCTION =====================================================================

cross_x = function(xs, ys) {
	idx = which(unlist(lapply(1:(length(ys)-1), FUN = function(i) {
		ys[i] * ys[i+1] < 0
	})))

	P1 = data.frame(x = xs[idx], y = ys[idx])
	P2 = data.frame(x = xs[idx+1], y = ys[idx+1])

	m = (P2$y - P1$y) / (P2$x - P1$x)
	i = P2$y - m * P2$x

	types = unlist(lapply(P1$y, FUN = function(x) {
		if( x >= 0 ) {
			return("maxima")
		} else {
			return("minima")
		}
	}))

	return(data.frame(i = idx + 1, x = -i/m, type = types))
}

fwhm = function(m, xs, ys) {
	hm = ys[xs == m] / 2

	crossings = cross_x(xs, ys - hm)

	fir = crossings$x[crossings$x <= m]
	fir = fir[length(fir)]
	sec = crossings$x[crossings$x >= m][1]

	return(2 * min(m - fir, sec - m))
}

calc_feature_range = function(data, column, k_sigma = 3, plot = F) {

	# Calculate starting points for fitting

	ddens = density(unlist(data[, ..column]))
	ddens = data.frame(x = ddens$x, density = ddens$y, stringsAsFactors = F)

	peaks = cross_x(ddens$x[-1], diff(ddens$density))
	maxima = peaks[peaks$type == "maxima",]
	maxima = maxima[order(ddens$density[maxima$i], decreasing = T),]
	peak1 = maxima[1,]
	peak2 = maxima[maxima$x > peak1$x,][1,]

	f1 = fwhm(ddens$x[peak1$i], ddens$x, ddens$density)
	sigma1 = f1 / sqrt(8 * log(2))
	f2 = fwhm(ddens$x[peak2$i], ddens$x, ddens$density)
	sigma2 = f2 / sqrt(8 * log(2))

	C1 = ddens$density[peak1$i]
	C2 = ddens$density[peak2$i]
	mean1 = peak1$x
	mean2 = peak2$x

	# Fit Sum of Gaussians
	
	model = tryCatch({
		nls(ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2)) +
			C2 * exp(-(ddens$x-mean2)**2/(2 * sigma2**2))), data=ddens,
			start=list(C1=C1, mean1=mean1, sigma1=sigma1,
			         C2=C2, mean2=mean2, sigma2=sigma2), algorithm="port") 
	}, error = function(e) {
		nls(ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2))),
			data=ddens, start=list(C1=C1, mean1=mean1, sigma1=sigma1),
			algorithm="port") 
	})
	cfs = as.data.frame(t(coef(model)), stringsAsFactors = T)
	ddens$G1 = cfs$C1*exp(-(ddens$x-cfs$mean1)**2/(2*cfs$sigma1**2))
	

	if ("C2" %in% names(cfs)) {
		# Enforce 2nd Gaussian to be AFTER and LOWER than the 1st
		# by reverting to single Gaussian fitting
		if (mean2 < (mean1 + 2 * sigma1) || C2 > C1) {
			model = nls(
				ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2))),
				data=ddens, start=list(C1=C1, mean1=mean1, sigma1=sigma1),
				algorithm="port") 
			cfs = as.data.frame(t(coef(model)), stringsAsFactors = T)
		}
	}
	if ("C2" %in% names(cfs))
		ddens$G2 = cfs$C2*exp(-(ddens$x-cfs$mean2)**2/(2*cfs$sigma2**2))

	return(list(
		range = c(mean1-k_sigma*sigma1, mean1+k_sigma*sigma1),
		density = ddens
	))
}

select_dataset_population = function(dataset, nt) {
	dt = nt[condition == dataset,]

	dna_sum_cfr = calc_feature_range(dt, "dna_sum", ksigma)
	volume_vx_cfr = calc_feature_range(dt, "volume_vx", ksigma)

	dt[, selected := "kept"]
	dt[dt$dna_sum < dna_sum_cfr$range[1] | dt$dna_sum > dna_sum_cfr$range[2], selected := "discarded"]
	dt[dt$volume_vx < volume_vx_cfr$range[1] | dt$volume_vx > volume_vx_cfr$range[2], selected := "discarded"]

	color = list(
		distrib = "#7F7F7F",
		range = "#b2182b",
		G1 = "#00770A",
		G2 = "#BF48CB",
		selected = "#D91C1C",
		discarded = "#2D50DC"
	)

	theme_set(theme_cowplot())

	p1 = ggplot(dt, aes(x = volume_vx, y = dna_sum, color = factor(selected))
		) + geom_point(size = .5, alpha = .5
		) + geom_vline(xintercept = volume_vx_cfr$range,
			color = color$range, linetype = "dashed"
		) + geom_hline(yintercept = dna_sum_cfr$range,
			color = color$range, linetype = "dashed"
		) + guides(color = F
		) + scale_color_manual(values = c(color$selected, color$discarded)
		) + xlab("Volume (vx)") + ylab("DNA stain intensity integral (a.u.)"
		) + ggtitle(sprintf("%s: %d/%d nuclei selected", dataset, dt[selected == "kept", .N], nrow(dt))
		) + coord_fixed(
			(max(dt$volume_vx)-min(dt$volume_vx))/(max(dt$dna_sum)-min(dt$dna_sum))
		)

	p2 = axis_canvas(p1, axis = "y", coord_flip = TRUE) + geom_line(
			data = dna_sum_cfr$density, aes(x = x, y = density),
			color = color$distrib
		) + geom_line(data = dna_sum_cfr$density, aes(x = x, y = G1),
			color = color$G1, alpha = 0.75, size = .5
		) + geom_vline(xintercept = dna_sum_cfr$range,
			color = color$range, linetype = "dashed",
			alpha = 0.75, size = .5
		) + coord_flip(
		)
	if ( "G2" %in% colnames(dna_sum_cfr$density) )
		p2 = p2 + geom_line(data = dna_sum_cfr$density, aes(x = x, y = G2),
			color = color$G2, alpha = 0.75, size = .5)

	p3 = axis_canvas(p1, axis = "x") + geom_line(
			data = volume_vx_cfr$density, aes(x = x, y = density),
			color = color$distrib
		) + geom_vline(xintercept = volume_vx_cfr$range,
			color = color$range, linetype = "dashed"
		) + geom_line(data = volume_vx_cfr$density, aes(x = x, y = G1),
			color = color$G1, alpha = 0.75, size = .5
		)
	if ( "G2" %in% colnames(volume_vx_cfr$density) )
		p3 = p3 + geom_line(data = volume_vx_cfr$density, aes(x = x, y = G2),
			color = color$G2, alpha = 0.75, size = .5)

	p12 <- insert_yaxis_grob(p1, p2, grid::unit(.2, "null"), position = "right")
	p123 <- insert_xaxis_grob(p12, p3, grid::unit(.2, "null"), position = "top")
	p = ggdraw(p123)

	data = data.table(
		dataset = dataset,
		volume_vx_min = volume_vx_cfr$range[1],
		volume_vx_max = volume_vx_cfr$range[2],
		dna_sum_min = dna_sum_cfr$range[1],
		dna_sum_max = dna_sum_cfr$range[2],
		nselected = dt[selected == "kept", .N],
		nnuclei = nrow(dt)
	)

	return(list(
		nuclei_table = dt,
		data = data,
		plot = p
	))
}

# RUN ==========================================================================

nt = fread(inpath)

reqCols = c("condition", "sid", "nid", "volume_vx", "volume_um3",
	"a2radius_px", "a2radius_um", "v2radius_px", "v2radius_um",
	"dna_sum", "dna_mean", "sig_sum", "sig_mean", "surface_um2", "sphericity")
stopifnot(all(reqCols %in% names(nt)))

cat("Nuclei counts:\n")
preCounts = nt[, .(total = .N), by = condition]
print(preCounts)

# Filter per dataset -----------------------------------------------------------

data = mclapply(sort(unique(nt$condition)),
	FUN = select_dataset_population, nt, mc.cores = threads)

outName = sprintf("%s.selected.%dsigma",
	tools::file_path_sans_ext(basename(inpath)), ksigma)

pdf(file.path(outpath, sprintf("%s.pdf", outName)),
	width = 8, height = 9)
p = lapply(data, FUN = function(x) {
	theme_set(theme_cowplot())
	print(x$plot)
})
graphics.off()

png(file.path(outpath, sprintf("%s.png", outName)), width = 1600, height = 800)
theme_set(theme_cowplot())
plot_grid(plotlist = lapply(data, function(x) x$plot))
graphics.off()

out = rbindlist(lapply(data, function(x) x$nuclei_table[selected == "kept"]))

postCounts = out[, .(kept = .N), by = condition]
setkeyv(preCounts, "condition")
setkeyv(postCounts, "condition")

logCounts = preCounts[postCounts]
logCounts[, removed := total-kept]
logCounts[, percKept := sprintf("%.2f%%", kept/total*100, 2)]
logCounts[, kSigma := ksigma]

cat("Selected nuclei:\n")
print(logCounts)

write.table(logCounts, file.path(outpath, sprintf("%s.log", outName)),
	quote = F, row.names = F, col.names = T, sep = "\t")
write.table(out, file.path(outpath, sprintf("%s.tsv", outName)),
	quote = F, row.names = F, col.names = T, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
