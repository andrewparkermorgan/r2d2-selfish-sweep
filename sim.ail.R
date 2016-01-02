#! /usr/bin/env Rscript

library(simcross, lib.loc = "~/lib/R")
library(plyr)
library(optparse)

source("~/lib/util/R/simcross.R")

opts <- list( make_option( c("-n","--nruns"), type = "integer", default = 1,
						   help = "number of simulation runs" ),
			  make_option( c("-b","--batch"), type = "integer", default = 1,
			  			 help = "batch number to be added to output" ),
			  make_option( c("-o","--outdir"), default = "./",
			  			 help = "directory in which to dump output" )
)
args <- parse_args(OptionParser(option_list = opts))
print(args)

## simulate HRxB6 AILs
ail.rez <- ldply(1:(args$nruns), function(i) {
	
	ped <- sim_ail_pedigree(ngen = 3, npairs = 20, design = "nosib")
	geno <- sim.until.fix(ped, npairs = 20, founders = 1)
	summ <- summarise.run(geno, founders = c(1,5))
	summ$run <- i
	summ$batch <- args$batch
	return(summ)
	
}, .progress = "text")

save(ail.rez, file = file.path(args$outdir, paste0("ail.batch", args$batch, ".Rdata")))