#Copyright 1998-2004 Timor Kadir.
#Kadir/Brady Feature detector (Scale Saliency) Code 
#Makefile for *nix-like systems

MEX=mex
MCC=mcc
CFLAGS=-O -DDEBUG -DFASTPLOG #-DSHOWPROG 
SRC= CalcEntropyScalesOpt.c Hists.c
all: CalcEntropyScalesOpt.mexglx

CalcScaleSaliency:CalcEntropyScalesOpt.c Hists.c
	$(MEX) $(CFLAGS) $^ -o $@	
clean:
	rm -f *.mexlx *.o *.mexglx 
	rm -Rf ScaleSaliency_Public

Public: CalcEntropyScalesOpt.c Hists.c 
	$(MEX) -V5 -DFASTPLOG -DPUBLICRELEASE $^ -o CalcScaleSaliency
	mkdir -p ScaleSaliency_Public
	cp CalcScaleSaliency.* ScaleSaliency_Public
	cp GreedyCluster.m ScaleSaliency_Public
	cp image_0001.jpg image_0001_out.jpg README drawcircle.m GoSalScale.m ScaleSaliency_Public
