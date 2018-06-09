#!/usr/bin/python
"""
This script makes the MC histogram for comparison to data prediction
"""
import RA2b
# import RA2b_rev as RA2b

import ROOT
outFile = ROOT.TFile('../plots/histograms/ZinvMCttzMC174bin_35.9ifb_test.root', "recreate")

## set the luminosity of the MC samples
doLumi = 35.9

## get Zinv and ttZ histograms
zvv = RA2b.getHist('zinv',applyPuWeight=False,doLumi=doLumi,removeZkfactor=False)
zvv.SetName(str(zvv.GetName()).replace("zinv", "zvv"))
print zvv.GetName()
zvv.SetDirectory(outFile.GetDirectory(""))
zvv.Draw()

ttzvv = RA2b.getHist('ttzvv',applyPuWeight=False,doLumi=doLumi)
ttzvv.SetDirectory(outFile.GetDirectory(""))
ttzvv.Draw()

# zinv = RA2b.getHist('zinv',treeLoc = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12",applyPuWeight="../plots/histograms/PileupHistograms_0121_69p2mb_pm4p6.root",doLumi=doLumi,removeZkfactor=False, doNoisy=True)
# ttzvv = RA2b.getHist('ttzvv',treeLoc = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12",applyPuWeight="../plots/histograms/PileupHistograms_0121_69p2mb_pm4p6.root",doLumi=doLumi)

## add the two samples
zinv = zvv.Clone()
zinv.SetName(str(zinv.GetName()).replace("zvv", "zinv"))
zinv.Add(ttzvv)
print zinv.GetName()
zinv.SetDirectory(outFile.GetDirectory(""))
zinv.Draw()

## store the histograms
# zinv.SaveAs('../plots/histograms/ZinvMCttzMC174bin_36.3ifb_KPPU.root')
# zvv.SaveAs('../plots/histograms/ZvvMC174bin_35.9ifb_noPU.root')
# ttzvv.SaveAs('../plots/histograms/ttzMC174bin_35.9ifb_noPU.root')
# zinv.SaveAs('../plots/histograms/ZinvMCttzMC174bin_35.9ifb_noPU.root')

outFile.Write()
