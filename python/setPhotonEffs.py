#!/usr/bin/python
"""Run this script to set the photon related efficiencies in
plots/histograms/effHists.root
All efficiency factors are hard coded below"""

import ROOT

########## photon inputs: (value, absolute error) #####
## photon efficiency scale factor
## apply to MC
######## obsolete, use SFcorrections.Photons.root #####
SF = (0.99, 0.01)

## photon fragmentation
## apply to data
######## obsolete, use fragmentation.root #############
frag = (0.92, 0.07)

## photon trigger efficiency 
## Andrew provided these (see Nov. 22 email)
## apply as a SF to MC (don't apply any simulated triggers to MC)
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<300, 300-350, 350-500, 500-750, 750+]
trig_eb = [(0.969, 0.002),
           (0.983, 0.001),
           (0.985, 0.001),
           (0.984, 0.001),
           (0.979, 0.004),]

trig_ec = [(0.953, 0.004),
           (0.974, 0.003),
           (0.984, 0.001),
           (0.989, 0.003),
           (0.980, 0.019),]

## photon purity (Andrew provides these)
## apply to data
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<225, 225-250, 250-300, 300-350, 350-500, 500+]
pur_eb = [(0.9580, 0.0267),
          (0.9623, 0.0028),
          (0.9686, 0.0026),
          (0.9571, 0.0225),
          (0.9491, 0.0493),
          (0.9524, 0.1010),]

pur_ec = [(0.8879, 0.0108),
          (0.8570, 0.0114),
          (0.8761, 0.0098),
          (0.8768, 0.0259),
          (0.9039, 0.0250),
          (0.9261, 0.0740),]
###################################################

########## get the efficiency file ################
effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")

########## update the efficiency file #############
h_SF_g1 = effFile.Get("h_SF_g1")
h_SF_g1.SetBinContent(1,SF[0])
h_SF_g1.SetBinError(1,SF[1])
h_SF_g1.Write(h_SF_g1.GetName(),2)

h_frag1 = effFile.Get("h_frag1")
h_frag1.SetBinContent(1,frag[0])
h_frag1.SetBinError(1,frag[1])
h_frag1.Write(h_frag1.GetName(),2)

h_trig_eb = effFile.Get("h_trig_eb")
for i in range(len(trig_eb)):
    h_trig_eb.SetBinContent(i+1,trig_eb[i][0])
    h_trig_eb.SetBinError(i+1,trig_eb[i][1])
h_trig_eb.Write(h_trig_eb.GetName(),2)

h_trig_ec = effFile.Get("h_trig_ec")
for i in range(len(trig_ec)):
    h_trig_ec.SetBinContent(i+1,trig_ec[i][0])
    h_trig_ec.SetBinError(i+1,trig_ec[i][1])
h_trig_ec.Write(h_trig_ec.GetName(),2)

h_pur_eb = effFile.Get("h_pur_eb")
for i in range(len(pur_eb)):
    h_pur_eb.SetBinContent(i+1,pur_eb[i][0])
    h_pur_eb.SetBinError(i+1,pur_eb[i][1])
h_pur_eb.Write(h_pur_eb.GetName(),2)

h_pur_ec = effFile.Get("h_pur_ec")
for i in range(len(pur_ec)):
    h_pur_ec.SetBinContent(i+1,pur_ec[i][0])
    h_pur_ec.SetBinError(i+1,pur_ec[i][1])
h_pur_ec.Write(h_pur_ec.GetName(),2)


effFile.Close()
