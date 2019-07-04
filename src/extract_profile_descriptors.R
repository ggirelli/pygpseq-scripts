#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190704
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

scriptName = "extract_profile_descriptors.R"
parser = arg_parser(paste0(
""), name = scriptName)

# parser = add_argument(parser, arg = 'inpath',
# 	help = 'Path to nuclei table.')

# parser = add_argument(parser, arg = '--ksigma', short = '-k', type = class(0),
# 	help = 'Sigma constant for interval definition.', default = 3, nargs = 1)
# parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
# 	help = 'Number of threads for parallelization.', default = 1, nargs = 1)
# parser = add_argument(parser, arg = '--outpath', short = "-O", type = class(""),
# 	help = 'Path to output folder. Default to inpath parent folder.')
# parser = add_argument(parser, arg = '--date-meta', short = "-d", type = class(""),
# 	help = paste0('Path to session meta table with "dataset", ',
# 		'"session" (id) and acqusition "date" columns.'))

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

# RUN ==========================================================================

# END --------------------------------------------------------------------------

################################################################################
