#!/usr/bin/python
"""
getExtrapDatFiles.py generates the extrapolation 
dat files for integration. Can be run from
command line, interactively or with the driver script:
scripts/extrapDatFileDriver.csh
"""
import sys
import ROOT
import RA2b

## doSample is a list of dat files to generate ##
doSample = []
## Check if samples are given at runtime (e.g. sig hdp ldp)
## can give more than one and will run over each sequentially
if len(sys.argv) > 1:
    for i in range(1,len(sys.argv)):
        doSample.append(sys.argv[i])
else:
    ## default runs over all 
    doSample = ['sig','hdp','ldp']
#################################################

for sample in doSample:

    ## define kinematic range ##
    kinRange = [] 
    ## if qcd binning add 11-13 to beginning of kinRange
    if(sample=='ldp' or sample=='hdp'):
        kinRange+=range(11,14)

    ## all samples use 10 nominal kinematic bins
    kinRange+=range(1,11)
    ############################

    ## the individual dilepton yields per bin
    zee = RA2b.getHist('zee',dphiCut=sample,kinRange=kinRange)
    zmm = RA2b.getHist('zmm',dphiCut=sample,kinRange=kinRange)
    
    ## the stat uncertainty on the extrapolation factors
    zll_Extrap = RA2b.getExtrapolation(['zee','zmm'],njSplit=False,doFactorization=True,kinRange=-1,dphiCut=sample)
    statErr = []
    for i in range(1,20):
        if(zll_Extrap.GetBinContent(i)==0):
            statErr.append(1.)
        else:
            statErr.append(zll_Extrap.GetBinError(i)/zll_Extrap.GetBinContent(i))
    
    ####################################################################
    ## get the purity uncertainties from the effHists.root file ########
    ####################################################################
    effFile = ROOT.TFile("../plots/histograms/effHists.root")
    
    h_pur_m = effFile.Get("h_pur_m")
    h_pur_e = effFile.Get("h_pur_e")
    
    ## NJets=2 purity as a function of nb
    pursyst_2j = [0]
    for i in range(2,4):
        avePur = (h_pur_e.GetBinError(i)/h_pur_e.GetBinContent(i)+h_pur_m.GetBinError(i)/h_pur_m.GetBinContent(i))/2
        pursyst_2j.append(avePur)

    ## NJets=3-4 purity as a function of nb
    pursyst_3to4j = [0]
    for i in range(5,8):
        avePur = (h_pur_e.GetBinError(i)/h_pur_e.GetBinContent(i)+h_pur_m.GetBinError(i)/h_pur_m.GetBinContent(i))/2
        pursyst_3to4j.append(avePur)
    
    ## NJets=5+ purity as a function of nb
    pursyst_5jplus = [0]
    for i in range(9,12):
        avePur = (h_pur_e.GetBinError(i)/h_pur_e.GetBinContent(i)+h_pur_m.GetBinError(i)/h_pur_m.GetBinContent(i))/2
        pursyst_5jplus.append(avePur)
    
    ## comine the uncertainties into one list to be called later
    purSyst = [pursyst_2j,pursyst_3to4j,pursyst_5jplus,pursyst_5jplus,pursyst_5jplus]
  
    ## This is needed else applying efficiencies later does weird things
    effFile.Close()
    
    ####################################################################
    ## end get purity uncertainties ####################################
    ####################################################################

    ## get the 0b photon data
    ## this is needed to do the split and merge type extrapolation
    ## using hdp even for ldp because it gets divided out
    ## and we gain from the increased stats 
    h_pho0b_split = RA2b.get0bPrediction('photon',kinRange=kinRange)
    h_pho0b_merged = RA2b.mergeHist(h_pho0b_split,kinRange=kinRange)
    
    ## get the extrapolation factors from data with purities applied
    h_zeemmExtrp = RA2b.getExtrapolation(['zee','zmm'],doFactorization=True,applyEffs=True,dphiCut=sample,kinRange=kinRange)

    ## clone the data extrapolation factors before MC correction for syst calc
    h_zeemmNoCorr = h_zeemmExtrp.Clone()
    
    ## get the correction factors from MC
    dy_mc_corr = RA2b.getMonteCarloCorrection(['dyee','dymm'],applySF=False,applyPuWeight=True,dphiCut=sample,kinRange=kinRange)
    dyll_mc_corr = RA2b.getMonteCarloCorrection(['dyee','dymm'],njSplit=False,kinRange=-1,applySF=False,applyPuWeight=True)
    
    ## get the binomial correction for syst calc
    h_binom=RA2b.getBinomialCorrection(kinRange=kinRange)
    
    ## apply the mc correction factors 
    h_zeemmExtrp.Multiply(dy_mc_corr)
    
    # to compute up and low syst
    h_zeemmExtrpNoTTZ = h_zeemmExtrp.Clone() 
    
    ## apply the binomial correction factors for syst calc
    h_binom.Multiply(h_zeemmNoCorr)
    
    ## clone the split photon sample as place holders for syst calc 
    h_pho0b_binom = h_pho0b_split.Clone()
    h_pho0b_nocorr = h_pho0b_split.Clone()
    h_pho0b_mid = h_pho0b_split.Clone()

    ## multiply by the various extrapolation factors
    h_pho0b_split.Multiply(h_zeemmExtrp)
    h_pho0b_mid.Multiply(h_zeemmExtrpNoTTZ)
    h_pho0b_binom.Multiply(h_binom)
    h_pho0b_nocorr.Multiply(h_zeemmNoCorr)
    
    ## merge the hists back to the 174 binning scheme
    h_pho_merged = RA2b.mergeHist(h_pho0b_split,kinRange=kinRange)
    h_pho_mid = RA2b.mergeHist(h_pho0b_mid,kinRange=kinRange)
    h_pho_upsyst = RA2b.mergeHist(h_pho0b_binom,kinRange=kinRange)
    h_pho_nocorr = RA2b.mergeHist(h_pho0b_nocorr,kinRange=kinRange)
    
    ## compute the upper and lower syst hists
    h_pho_dnsyst = h_pho_mid.Clone()
    h_pho_upsyst.Add(h_pho_mid,-1)
    h_pho_dnsyst.Add(h_pho_nocorr,-1)
    
    ## divide mid to get fractional error
    h_pho_upsyst.Divide(h_pho_mid)
    h_pho_dnsyst.Divide(h_pho_mid)
    
    ## scale by 1/sqrt(3) : full width at half max
    h_pho_upsyst.Scale(1/3.**0.5)
    h_pho_dnsyst.Scale(1/3.**0.5)
    
    ## divide out the photon normalization
    ## this gives the extrapolation factors themselves
    h_extrap = h_pho_merged.Clone()
    h_extrap.Divide(h_pho0b_merged)
    
    ## kin systematics derived using
    ## analysis of the closure plot
    systKin = [0,0.07,0.10,0.20]
    
    ## formatting dictionaries 
    njDict = {1: ' 2 ',
              2: '3-4',
              3: '5-6',
              4: '7-8',
              5: '9+ '}
    nbDict = {0: '  0   ',
              1: '  1   ',
              2: '  2   ',
              3: ' >=3  '}
    
    ## this block is here to ensure that even when the corresponding
    ## photon 0b bin is zero, the extrapolation factors are non zero

    ## remove  low HT bins
    kinRangeCut = list(kinRange)
    for i in [11,1,4]:
        if i in kinRangeCut:
            kinRangeCut.remove(i)

    ## first NJets=9 bin
    Bin = 3*4*len(kinRange)+4*len(kinRangeCut)-len(kinRange)+1

    nonZeroBin = Bin-len(kinRangeCut)
    for nj in [5]:
        for nb in range(4):
            nonZeroBin += len(kinRangeCut)
            for kin in kinRangeCut:
                h_extrap_9p = h_extrap.GetBinContent(nonZeroBin)
                h_upsyst_9p = h_pho_upsyst.GetBinContent(nonZeroBin)
                h_dnsyst_9p = h_pho_dnsyst.GetBinContent(nonZeroBin)
                h_extrap.SetBinContent(Bin,h_extrap_9p)
                h_pho_upsyst.SetBinContent(Bin,h_upsyst_9p)
                h_pho_dnsyst.SetBinContent(Bin,h_dnsyst_9p)
                Bin+=1

    
    ## had to add this funny business to 
    ## remove the bins that don't have
    ## any events
    removedBins = RA2b.getRemovedBins(kinRange=kinRange)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]
    
    uncutBin = 0
    Bin = 1
    nbnjSubtract = 0
    print "Njbin | Nbbin | Nmumu |  Nee  | Nb/0b   | stat | MC stat | syst+ | syst- | sysKin | sysPur"
    for nj in range(1,6):
        for nb in range(4):
            if(nj==1 and nb==3):
                nbnjSubtract+=1
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue
                line = " "+njDict[nj]+"  |  "+nbDict[nb]+"  |   "
                line+= str(zmm.GetBinContent(Bin))+"  |"+str(zee.GetBinContent(Bin))
                line+= " | "+str(h_extrap.GetBinContent(Bin))
                line+= " | "+str(statErr[(nj-1)*4+nb-nbnjSubtract])
                line+= " | "+str(dyll_mc_corr.GetBinError((nj-1)*4+nb+1-nbnjSubtract)/max(dyll_mc_corr.GetBinContent((nj-1)*4+nb+1-nbnjSubtract),1))
                line+= " | "+str(abs(h_pho_upsyst.GetBinContent(Bin)))
                line+= " | "+str(abs(h_pho_dnsyst.GetBinContent(Bin)))
                line+= " | "+str(systKin[nb])
                line+= " | "+str(purSyst[nj-1][nb])
                print line
                Bin+=1
