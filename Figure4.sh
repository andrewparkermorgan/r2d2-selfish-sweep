#! /bin/bash

OUTDIR=./simulations

## simulate HR selection lines (cf. Swallow et al (1998) Behav Genet 28:227-237)
./sim.hr.R --outdir $OUTDIR --nruns 1000

## simulate the HRxC57BL/6J advanced intercross line (cf. Kelley et al (2010) Physiol Genomics 42:190-200)
./sim.ail.R --outdir $OUTDIR --nruns 1000

## make expanded version of Figure 4
./Figure4.R