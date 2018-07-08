#!/usr/bin/python
"""
This script makes the MC histogram for comparison to data prediction
"""
import RA2b

## set the luminosity of the MC samples
doLumi = 36.3

## get Zinv and ttZ histograms
zinv = RA2b.getHist('zinv',applyPuWeight=True,doLumi=doLumi,removeZkfactor=False)
ttzvv = RA2b.getHist('ttzvv',applyPuWeight=True,doLumi=doLumi)

## add the two samples
zinv.Add(ttzvv)

## store the histogrma 
zinv.SaveAs('../plots/histograms/ZinvMCttzMC174bin.root')
