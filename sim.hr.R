#! /usr/bin/env Rscript

library(simcross)
library(plyr)
library(optparse)

source("simcross.R")

opts <- list( make_option( c("-n","--nruns"), type = "integer", default = 1,
						   help = "number of simulation runs" ),
			  make_option( c("-b","--batch"), type = "integer", default = 1,
			  			 help = "batch number to be added to output" ),
			  make_option( c("-o","--outdir"), default = "./",
			  			 help = "directory in which to dump output" )
			  )
args <- parse_args(OptionParser(option_list = opts))
print(args)

## simulate HR lines
hr.rez <- ldply(1:(args$nruns), function(i) {

	#message("--- run ", i, " ---")
	ped <- outbred.from.founders(ngen = 1, npairs = 10, nfounders = 4, design = "nosib")
	geno <- sim.until.fix(ped, npairs = 10, founders = c(1:3,5:7))
	summ <- summarise.run(geno, founders = c(1:3,5:7))
	summ$run <- i
	summ$batch <- args$batch
	return(summ)

}, .progress = "text")

## save result
save(hr.rez, file = file.path(args$outdir, paste0("hr.batch", args$batch, ".Rdata")))
