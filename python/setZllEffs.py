#!/usr/bin/python
"""Run this script to set the zll related efficiencies in
plots/histograms/effHists.root
Purity fits are ran using RA2b.getZmassFitPlot()
current binning for purity is 
(NJ_2: Nb0, Nb1, Nb2)
(NJ_3-4:  Nb0, Nb1, Nb2+) 
(NJ_5+:  Nb0, Nb1, Nb2+) 
Trigger efficiencies are hard coded below.
Lepton scale factors are two separate files
proveded by Frank."""

import RA2b
import ROOT

########## trigger effs from manuel ##############
trig_m = (0.988,0.004)
trig_e = (0.988,0.004)

########## run fits to get purity ################
fit_2j = []
fit_3to4j = []
fit_5jplus = []

# nb=12 is nb>=2
for nb in [0,1,12]:
    fit_2j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=1))
    fit_3to4j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=2))
    fit_5jplus.append(RA2b.getZmassFitPlot(bJetBin=nb,extraCuts='NJets>=5'))

fits = [fit_2j,fit_3to4j,fit_5jplus,fit_5jplus,fit_5jplus]

########## get the efficiency file ################
effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")


######### set the purities found above ############
h_pur_m = effFile.Get("h_pur_m")
h_pur_e = effFile.Get("h_pur_e")

Bin = 1
for nj in range(1,6):
    for nb in range(4):
        if(nb==3):
            if(nj==1):
                continue
            else:
                h_pur_m.SetBinContent(Bin,fits[nj-1][nb-1][0][0])
                h_pur_m.SetBinError(Bin,fits[nj-1][nb-1][0][1])
                h_pur_e.SetBinContent(Bin,fits[nj-1][nb-1][1][0])
                h_pur_e.SetBinError(Bin,fits[nj-1][nb-1][1][1])
        else:
            h_pur_m.SetBinContent(Bin,fits[nj-1][nb][0][0])
            h_pur_m.SetBinError(Bin,fits[nj-1][nb][0][1])
            h_pur_e.SetBinContent(Bin,fits[nj-1][nb][1][0])
            h_pur_e.SetBinError(Bin,fits[nj-1][nb][1][1])
        Bin+=1
h_pur_m.Write(h_pur_m.GetName(),2)
h_pur_e.Write(h_pur_e.GetName(),2)

######### set the trig effs from manuel ############
h_trig_m1 = effFile.Get("h_trig_m1")
h_trig_e1 = effFile.Get("h_trig_e1")

h_trig_m1.SetBinContent(1,trig_m[0])
h_trig_m1.SetBinError(1,trig_m[1])
h_trig_e1.SetBinContent(1,trig_e[0])
h_trig_e1.SetBinError(1,trig_e[1])

h_trig_m1.Write(h_trig_m1.GetName(),2)
h_trig_m1.Write(h_trig_e1.GetName(),2)

effFile.Close()

