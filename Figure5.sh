#! /bin/bash

## do simulations of meiotic drive with modifier alleles
OUTDIR=./simulations
NRUNS=100

./simtrd.py --config scenarios/m0.yaml --nruns $NRUNS --out $OUTDIR/m0
mv $OUTDIR/m0.summary.txt $OUTDIR/simulations/n100.t5.s80.p80.m0.txt

./simtrd.py --config scenarios/m0.yaml --nruns $NRUNS --out $OUTDIR/m1
mv $OUTDIR/m1.summary.txt $OUTDIR/simulations/n100.t5.s80.p80.m1.txt

./simtrd.py --config scenarios/m0.yaml --nruns $NRUNS --out $OUTDIR/m2
mv $OUTDIR/m2.summary.txt $OUTDIR/simulations/n100.t5.s80.p80.m2.txt

./simtrd.py --config scenarios/m0.yaml --nruns $NRUNS --out $OUTDIR/m3
mv $OUTDIR/m3.summary.txt $OUTDIR/simulations/n100.t5.s80.p80.m3.txt

./simtrd.py --config scenarios/m0.yaml --nruns $NRUNS --out $OUTDIR/m4
mv $OUTDIR/m4.summary.txt $OUTDIR/simulations/n100.t5.s80.p80.m4.txt

# make Figure 5
./Figure5.R
