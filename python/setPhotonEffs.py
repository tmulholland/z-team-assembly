#!/usr/bin/python
"""Run this script to set the photon related efficiencies in
plots/histograms/effHists.root
All efficiency factors are hard coded below"""

import ROOT

########## photon inputs: (value, absolute error) #####
## photon efficiency scale factor
## apply to MC
SF = (0.99, 0.01)

## photon fragmentation
## apply to data
frag = (0.92, 0.07)

## photon trigger efficiency (Kevin S. provided these)
## apply as a SF to MC (don't apply any simulated triggers to MC)
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<300, 300-350, 350-500, 500-750, 750+]
trig_eb = [(0.98821, 0.00230),
           (0.98692, 0.00293),
           (0.98347, 0.00246),
           (0.97952, 0.00562),
           (0.96992, 0.02312),]

trig_ec = [(0.96296, 0.00438),
           (0.94843, 0.00567),
           (0.93808, 0.00508),
           (0.96579, 0.01191),
           (1.00000, 0.06836),]

## photon purity (Andrew provides these)
## apply to data
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<225, 225-250, 250-300, 300-350, 350-500, 500+]
pur_eb = [(0.9580, 0.0267),
          (0.9623, 0.0028),
          (0.9686, 0.0026),
          (0.9557, 0.0343),
          (0.9709, 0.0230),
          (0.9873, 0.0116),]

pur_ec = [(0.8879, 0.0108),
          (0.8570, 0.0114),
          (0.8761, 0.0098),
          (0.8748, 0.0246),
          (0.9044, 0.0281),
          (0.9467, 0.0330),]
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
