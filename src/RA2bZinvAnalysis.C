//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#include "RA2bZinvAnalysis.h"
#include "fillFileMap.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
// #include <TLegend.h>
// #include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TCut.h>
#include <TTreeCache.h>

#include <TMath.h>
using TMath::Sqrt; using TMath::Power;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>

ClassImp(RA2bZinvAnalysis)

// ======================================================================================

RA2bZinvAnalysis::RA2bZinvAnalysis() :
  era_("2016"),
  ntupleVersion_("V12"),
  isSkim_(true),
  isMC_(false),
  intLumi_(1),
  treeLoc_(""),
  treeName_("tree"),
  deltaPhi_("nominal"),
  // deltaPhi_("hdp"),
  applyMassCut_(true),
  applyPtCut_(true),
  applyMinDeltaRCut_(true),
  // applySF_(false),
  // njSplit_(false),
  useTreeCCbin_(true),  // only in skims
  applyBTagSF_(false),  // overridden false if !isMC_
  applyPuWeight_(true),  // overridden false if !isMC_
  customPuWeight_(true),  // Substitute Kevin P recipe for the PuWeight in the tree
  puWeight(1),  // overridden from tree if isMC_
  Weight(1),  // overridden from tree if isMC_
  TrueNumInteractions(20),  // overridden from tree if isMC_
  RA2bin(0),  // overridden from tree if isSkim
  NElectrons(0),  // overriden from tree if >=V15
  NMuons(0)  // overriden from tree if >=V15
{
  if (!isMC_) {
    applyBTagSF_ = false;
    applyPuWeight_ = false;
  }
  if (!isSkim_) {
    useTreeCCbin_ = false;
    treeName_ = "TreeMaker2/PreSelection";  // For ntuple
  } else {
    treeName_ = "tree";  // For skims
  }
  if (ntupleVersion_ == "V12") {
    // treeLoc_ = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12";  // Colorado, owned by wtford (Zjets only)
    treeLoc_ = "/nfs/data38/cms/mulholland/lpcTrees/Skims/Run2ProductionV12";  // Colorado, owned by mulholland
    // treeLoc_ = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV12";  // xrootd
    // treeLoc_ = "/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV12";  // from cmslpc
  } else if (ntupleVersion_ == "V15") {
    treeLoc_ = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV15";  // ntuples, xrootd
    // treeLoc_ = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV15";  // xrootd
    // treeLoc_ = "/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV15";  // from cmslpc
    // tmt_ = new TreeMkr_unskimmed_data_V15;
  }
  if (era_ == TString("2016")) {
    intLumi_ = 35.9;

    if (applyPuWeight_ && customPuWeight_) {
      TFile* pufile = TFile::Open("../../Analysis/corrections/PileupHistograms_0121_69p2mb_pm4p6.root","READ");
      puHist_ = (TH1*) pufile->Get("pu_weights_down");
    }

    if (isMC_) BTagSFfile_ = "../../Analysis/btag/CSVv2_Moriond17_B_H_mod.csv";

    // Needed branches
    activeBranches_.push_back("NJets");
    activeBranches_.push_back("BTags");
    activeBranches_.push_back("HT");
    activeBranches_.push_back("MHT");
    activeBranches_.push_back("JetID");
    activeBranches_.push_back("Jets");
    activeBranches_.push_back("Jets_hadronFlavor");
    activeBranches_.push_back("Jets_HTMask");
    activeBranches_.push_back("isoElectronTracks");
    activeBranches_.push_back("isoMuonTracks");
    activeBranches_.push_back("isoPionTracks");
    activeBranches_.push_back("DeltaPhi1");
    activeBranches_.push_back("DeltaPhi2");
    activeBranches_.push_back("DeltaPhi3");
    activeBranches_.push_back("DeltaPhi4");
    if (isSkim_) {
      activeBranches_.push_back("RA2bin");
    } else {
      activeBranches_.push_back("NJetsclean");
      activeBranches_.push_back("BTagsclean");
      activeBranches_.push_back("HTclean");
      activeBranches_.push_back("MHTclean");
      activeBranches_.push_back("JetIDclean");
      activeBranches_.push_back("Jetsclean");
      activeBranches_.push_back("Jetsclean_hadronFlavor");
      activeBranches_.push_back("Jetsclean_HTMask");
      activeBranches_.push_back("isoElectronTracksclean");
      activeBranches_.push_back("isoMuonTracksclean");
      activeBranches_.push_back("isoPionTracksclean");
      activeBranches_.push_back("DeltaPhi1clean");
      activeBranches_.push_back("DeltaPhi2clean");
      activeBranches_.push_back("DeltaPhi3clean");
      activeBranches_.push_back("DeltaPhi4clean");
    }
    if (ntupleVersion_ != "V12") {
      activeBranches_.push_back("NMuons");
      activeBranches_.push_back("NElectrons");
    }
    activeBranches_.push_back("Muons");
    activeBranches_.push_back("Electrons");
    activeBranches_.push_back("ZCandidates");
    activeBranches_.push_back("Photons");
    activeBranches_.push_back("Photons_nonPrompt");
    activeBranches_.push_back("Photons_fullID");
    activeBranches_.push_back("NVtx");
    activeBranches_.push_back("TriggerPass");
    activeBranches_.push_back("TriggerPrescales");
    activeBranches_.push_back("HBHENoiseFilter");
    activeBranches_.push_back("HBHEIsoNoiseFilter");
    activeBranches_.push_back("eeBadScFilter");
    activeBranches_.push_back("EcalDeadCellTriggerPrimitiveFilter");
    activeBranches_.push_back("globalTightHalo2016Filter");
    activeBranches_.push_back("BadChargedCandidateFilter");
    activeBranches_.push_back("BadPFMuonFilter");
    if (isMC_) {
      activeBranches_.push_back("puWeight");
      activeBranches_.push_back("Weight");
      activeBranches_.push_back("madMinPhotonDeltaR");
    }

    kinThresholds_.push_back({300, 300, 500, 1000});  // mht threshold, {ht thresholds}
    kinThresholds_.push_back({350, 350, 500, 1000});
    kinThresholds_.push_back({500, 500, 1000});
    kinThresholds_.push_back({750, 750, 1500});
    kinThresholds_.push_back({250, 300, 500, 1000}); // QCD control bins
    nJetThresholds_ = {2, 3, 5, 7, 9};
    nbThresholds_ = {0, 1, 2, 3};

    Int_t bin = 0;
    for (unsigned j = 0; j < nJetThresholds_.size(); ++j) {
      for (unsigned b = 0; b < nbThresholds_.size(); ++b) {
	if (nbThresholds_[b] > nJetThresholds_[j]) continue;  // Exclude (Njets0, Nb3)
	unsigned mmax = deltaPhi_ == TString("nominal") ? kinThresholds_.size()-1 : kinThresholds_.size();
	int k = -1;
	for (unsigned m = 0; m < mmax; ++m) {
	  for (unsigned h = 1; h < kinThresholds_[m].size(); ++h) {
	    k++;
	    if (j > 2 && (m < 2 || m == kinThresholds_.size()-1) && h == 1) continue;   // Exclude (Njets3,4; HT0,3,(6))
	    std::vector<int> jbk = {int(j), int(b), k};
	    bin++;
	    toCCbin_[jbk] = bin;
	    // cout << "Filling toCCbin; j = " << j << ", b = "  << b << ", k = " << k << ", bin = " << bin << endl;
	  }
	}
      }
    }
  } // era 2016

  kinSize_ = 0;
  int mmax = (deltaPhi_ == TString("nominal")) ? kinThresholds_.size() -1 : kinThresholds_.size();
  for (int i = 0; i < mmax; ++i)
    kinSize_ += kinThresholds_[i].size();

  fillFileMap();
  fillCutMaps();

}  // ======================================================================================

TChain*
RA2bZinvAnalysis::getChain(const char* sample, Int_t* fCurrent, bool setBrAddr) {
  TString theSample(sample);
  TString key;
  if (theSample.Contains("zinv")) key = TString("zinv");
  else if (theSample.Contains("ttzvv")) key = TString("ttzvv");
  else if (theSample.Contains("dymm")) key = TString("dymm");
  else if (theSample.Contains("dyee")) key = TString("dyee");
  else if (theSample.Contains("ttzmm")) key = TString("ttzmm");
  else if (theSample.Contains("ttzee")) key = TString("ttzee");
  else if (theSample.Contains("VVmm")) key = TString("VVmm");
  else if (theSample.Contains("VVee")) key = TString("VVee");
  else if (theSample.Contains("ttmm")) key = TString("ttmm");
  else if (theSample.Contains("ttee")) key = TString("ttee");
  else if (theSample.Contains("zmm")) key = TString("zmm");
  else if (theSample.Contains("zee")) key = TString("zee");

  TChain* chain = new TChain(treeName_);
  std::vector<TString> files;
  try {files = fileMap_.at(key);}
  catch (const std::out_of_range& oor) {
    std::cerr << oor.what() << " getChain: sample key (" << key << ") not found in fileMap" << endl;
    exit(2);
  }
  for (auto file : files) {
    cout << file << endl;
    chain->Add(file);
  }
  cout << "Initial size of cache for chain = " << chain->GetCacheSize() << endl;
  TTreeCache::SetLearnEntries(1);
  chain->SetCacheSize(200*1024*1024);
  chain->SetCacheEntryRange(0, chain->GetEntries());
  for (auto theBranch : activeBranches_) chain->AddBranchToCache(theBranch, true);
  // chain->AddBranchToCache("*", true);
  chain->StopCacheLearningPhase();
  cout << "Reset size of cache for chain = " << chain->GetCacheSize() << endl;

  if (fCurrent != nullptr) *fCurrent = -1;
  // if (setBrAddr) tmt_->Init(chain);
  if (setBrAddr) setBranchAddress(chain);

  chain->SetBranchStatus("*", 0);  // disable all branches
  for (auto theBranch : activeBranches_) chain->SetBranchStatus(theBranch, 1);

  return chain;
}  // ======================================================================================

TCut
RA2bZinvAnalysis::getCuts(const TString sample) {

  TCut cuts;
  TString sampleKey;
  try {sampleKey = sampleKeyMap_.at(sample);}
  catch (const std::out_of_range& oor) {
    std::cerr << oor.what() << " getCuts called with invalid sample = " << sample << endl;
    return cuts;
  }

  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut_ = "ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut_ = "ZCandidates.Pt()>=200.";
    // ptCut_ = "ZCandidates.Pt()>=100.";  // Troy revision

  TString photonDeltaRcut;
  std::vector<TString> trigger;
  try {trigger = triggerMap_.at(sample);}
  catch (const std::out_of_range& oor) {trigger.clear();}

  // commonCuts_ = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  if (trigger.empty()) {
    commonCuts_ = "JetID==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-
  } else {
    commonCuts_ = "JetID==1 && globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-
    int Ntrig = trigger.size();
    if (Ntrig > 1) trigCuts_ += TString("(");
    for (auto theTrigger : trigger)
      trigCuts_ += TString("(TriggerPass[")+theTrigger+TString("]==1) + ");
    if (Ntrig > 1) trigCuts_.Replace(trigCuts_.Length()-3, 3, ")");
  }
  // cout << "trigCuts_ = " << trigCuts_ << endl;

  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "Photons[0].Pt()>=200.";
    if (trigger.empty() && applyMinDeltaRCut_) photonDeltaRcut = "madMinPhotonDeltaR>=0.4";
  }
				       
    // 	if(extraCuts!=None):
    //     cuts+=extraCuts

  HTcut_ = std::string("HT>=") + std::to_string(kinThresholds_[0][1]);
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  NJetscut_ = std::string("NJets>=") + std::to_string(nJetThresholds_[0]);
  objcut_ = objCutMap_.at(sampleKey);
  minDphicut_ = minDphiCutMap_.at(deltaPhi_);

  cuts += objcut_;
  cuts += HTcut_;
  cuts += NJetscut_;
  cuts += MHTcut_;
  cuts += minDphicut_;
  cuts += massCut_;
  cuts += ptCut_;
  cuts += photonDeltaRcut;
  cuts += commonCuts_;
  cuts += trigCuts_;
    // 	  if(applySF):
    //     cuts*=bJetCutsSF[bJetBin]
    // 	else:
    //     cuts+=bJetCuts[bJetBin]
  if(applyPuWeight_ && !customPuWeight_) cuts *= "puWeight*(1)";
    // 	  if(type(extraWeight) is str):
    //     extraWeight+="*(1)"
    //     cuts*=extraWeight

  return cuts;
 
}  // ======================================================================================

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //
  TCut baselineCuts = getCuts(sample);
  cout << "baseline = " << endl << baselineCuts << endl;
  Int_t fCurrent; //!current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);
  TObjArray* forNotify = new TObjArray;

  cutHistos cutHistFiller(chain, forNotify);
  for (auto & hg : histograms) {
    hg->hist = new TH1F(hg->name, hg->title, hg->Nbins, hg->lowEdge, hg->highEdge);
    hg->hist->GetXaxis()->SetTitle(hg->axisTitles.first);
    hg->hist->GetYaxis()->SetTitle(hg->axisTitles.second);
    hg->hist->SetOption("HIST");
    hg->hist->SetMarkerSize(0);
    hg->NminusOneCuts = baselineCuts;
    for (auto cutToOmit : hg->omitCut) hg->NminusOneCuts(*cutToOmit) = "1";
    cout << "For sample " << sample << ", histo " << hg->name  << ", hg->omitCut = ";
    for (auto cutToOmit : hg->omitCut) cout << *cutToOmit << " ";
    cout << ", cuts = " << endl << hg->NminusOneCuts << endl;
    hg->NminusOneFormula = new TTreeFormula(hg->name, hg->NminusOneCuts, chain);
    forNotify->Add(hg->NminusOneFormula);
  }
  chain->SetNotify(forNotify);

  Long64_t Nentries = chain->GetEntries();
  cout << "Nentries in tree = " << Nentries << endl;
  int count = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (count % 100000 == 0) cout << "Entry number " << count << endl;

    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) cout << "Current file in chain: " << thisFile->GetName() << endl;
    }
    chain->GetEntry(entry);
    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    // cout << "NJets, Jets.size = " << NJets << ", " << Jets->size() << endl;
    // cout << "          " << HT << ",            " << MHT << endl;
    // cout << "ev" << HT << "  " << MHT << "  " << NJets << "  " << BTags << endl;

    Double_t eventWt = 1;
    Double_t PUweight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    if (isMC_) {
      eventWt = 1000*intLumi_*Weight*PUweight;
      if (eventWt < 0) eventWt *= -1;
    }


    for (auto & hg : histograms) {
      if (hg->name.Contains(TString("hCut"))) {
	cutHistFiller.fill(hg->hist, eventWt);
	continue;
      }
      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt != 0) {
	if (hg->dvalue != nullptr) {
	  // cout << "  " << *(hg->dvalue);
	  hg->hist->Fill(*(hg->dvalue), selWt*eventWt);
	}
	else if (hg->ivalue != nullptr) {
	  // cout << "  " << *(hg->ivalue);
	  hg->hist->Fill(Double_t(*(hg->ivalue)), selWt*eventWt);
	}
	else if (hg->filler != nullptr) (this->*(hg->filler))(hg->hist, selWt*eventWt);
	else cerr << "No method to fill histogram provided for " << hg->name << endl;
      }
    }  // loop over histograms
    // cout << endl;
  }  // loop over entries

  delete forNotify;
  delete chain->GetCurrentFile();

}  // ======================================================================================

std::vector<TH1F*>
RA2bZinvAnalysis::makeHistograms(const char* sample) {
  //
  // Define histograms, variables to fill, and cuts to be modified.
  // Return a vector of histograms.
  // For an event variable of type double (int), set member dvalue (ivalue)
  // to point to the tree variable.  For other types of tree branches, supply
  // a function to fill the histogram, and set member filler to point to it.
  //
  // When drawing, control number of digits in axis labels by
  // TGaxis myTGaxis;
  // myTGaxis.SetMaxDigits(4);

  std::vector<hist1D*> histograms;

  hist1D hHT;
  hHT.name = TString("hHT_") + TString(sample);  hHT.title = "HT";
  hHT.Nbins = 60;  hHT.lowEdge = 0;  hHT.highEdge = 3000;
  hHT.axisTitles.first = "HT [GeV]";  hHT.axisTitles.second = "Events/50 GeV";
  hHT.dvalue = &HT;  hHT.omitCut.push_back(&HTcut_);
  histograms.push_back(&hHT);

  hist1D hMHT;
  hMHT.name = TString("hMHT_") + TString(sample);  hMHT.title = "MHT";
  hMHT.Nbins = 60;  hMHT.lowEdge = 0;  hMHT.highEdge = 3000;
  hMHT.axisTitles.first = "MHT [GeV]";  hMHT.axisTitles.second = "Events/50 GeV";
  hMHT.dvalue = &MHT;  hMHT.omitCut.push_back(&MHTcut_);  hMHT.omitCut.push_back(&ptCut_);
  histograms.push_back(&hMHT);

  hist1D hNJets;
  hNJets.name = TString("hNJets_") + TString(sample);  hNJets.title = "NJets";
  hNJets.Nbins = 20;  hNJets.lowEdge = 0;  hNJets.highEdge = 20;
  hNJets.axisTitles.first = "N (jets)";  hNJets.axisTitles.second = "Events/bin";
  hNJets.ivalue = &NJets;  hNJets.omitCut.push_back(&NJetscut_);
  histograms.push_back(&hNJets);

  hist1D hBTags;
  hBTags.name = TString("hBTags_") + TString(sample);  hBTags.title = "BTags";
  hBTags.Nbins = 20;  hBTags.lowEdge = 0;  hBTags.highEdge = 20;
  hBTags.axisTitles.first = "N (b jets)";  hBTags.axisTitles.second = "Events/bin";
  hBTags.ivalue = &BTags;
  histograms.push_back(&hBTags);

  hist1D hZmass;
  hZmass.name = TString("hZmass_") + TString(sample);  hZmass.title = "Z mass";
  hZmass.Nbins = 30;  hZmass.lowEdge = 60;  hZmass.highEdge = 120;
  hZmass.axisTitles.first = "M(Z) [GeV]";  hZmass.axisTitles.second = "Events/2 GeV";
  hZmass.filler = &RA2bZinvAnalysis::fillZmass;  hZmass.omitCut.push_back(&massCut_);
  histograms.push_back(&hZmass);

  hist1D hZpt;
  hZpt.name = TString("hZpt_") + TString(sample);  hZpt.title = "Z Pt";
  hZpt.Nbins = 60;  hZpt.lowEdge = 0;  hZpt.highEdge = 3000;
  hZpt.axisTitles.first = "Pt(Z) [GeV]";  hZpt.axisTitles.second = "Events/50 GeV";
  hZpt.filler = &RA2bZinvAnalysis::fillZpt;  hZpt.omitCut.push_back(&ptCut_);  hZpt.omitCut.push_back(&MHTcut_);
  histograms.push_back(&hZpt);

  hist1D hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + TString(sample);  hCutFlow.title = "Cut flow";
  hCutFlow.Nbins = 10;  hCutFlow.lowEdge = 0;  hCutFlow.highEdge = 10;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  histograms.push_back(&hCutFlow);

  hist1D hCuts;
  hCuts.name = TString("hCuts_") + TString(sample);  hCuts.title = "Cuts passed";
  hCuts.Nbins = 10;  hCuts.lowEdge = 0;  hCuts.highEdge = 10;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  histograms.push_back(&hCuts);

  bookAndFillHistograms(sample, histograms);

  std::vector<TH1F*> theHists;
  theHists.push_back(hHT.hist);
  theHists.push_back(hMHT.hist);
  theHists.push_back(hNJets.hist);
  theHists.push_back(hBTags.hist);
  theHists.push_back(hZmass.hist);
  theHists.push_back(hZpt.hist);

  return theHists;

}  // ======================================================================================

TH1F*
RA2bZinvAnalysis::makeCChist(const char* sample) {
  //
  // Book the analysis-binned histogram, fill, and return it.
  //
  Int_t MaxBins = toCCbin_.size();
  cout << "MaxBins = " << MaxBins << endl;

  TH1F* hCCbins = new TH1F("hCCbins", "Zinv background estimate", MaxBins, 0.5, MaxBins+0.5);
  hCCbins->SetOption("HIST");
  hCCbins->SetMarkerSize(0);

  TObjArray* forNotify = new TObjArray;  // Allow for more than one TObject to notify of a new file
  
  Int_t fCurrent; //!current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);

  // chain->Print();

  BTagCorrector* btagcorr = nullptr;
  if (applyBTagSF_) {
    btagcorr = new BTagCorrector;
    btagcorr->SetCalib(BTagSFfile_);
  }

  // Get the baseline cuts, make a TreeFormula, and add it to the list
  // of objects to be notified as new files in the chain are encountered
  TCut baselineCuts = getCuts(sample);
  cout << "baseline = " << baselineCuts << endl;
  // See https://root-forum.cern.ch/t/how-to-evaluate-a-formula-for-a-single-tree-event/12283
  TTreeFormula* baselineTF = new TTreeFormula("baselineTF", baselineCuts, chain);
  forNotify->Add(baselineTF);
  chain->SetNotify(forNotify);

  // Traverse the events in the chain
  Long64_t Nentries = chain->GetEntries();
  cout << "Nentries in tree = " << Nentries << endl;
  int count = 0, outCount = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (count < 20 || count % 100000 == 0) cout << "Entry number " << count << endl;
    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) {
    	cout << "Current file in chain: " << thisFile->GetName() << endl;
    	if (btagcorr) btagcorr->SetEffs(thisFile);
      }
    }
    chain->GetEntry(entry);
    cleanVars();  // If unskimmed input, copy <var>clean to <var>

    UInt_t binCC = 0;

    // Apply baseline selection
    baselineTF->GetNdata();
    double selWt = baselineTF->EvalInstance(0);
    if (selWt == 0) continue;
    outCount++;

    // Compute event weight factors
    Double_t eventWt = 1;
    Double_t PUweight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    if (count < 20 || count % 10000 == 0) cout << "PUweight = " << PUweight << endl;
    if (isMC_) {
      eventWt = 1000*intLumi_*Weight*selWt*PUweight;
      if (eventWt < 0) eventWt *= -1;
    }
    if (count < 20 || count % 10000 == 0) cout << "eventWt = " << eventWt << endl;

    if (useTreeCCbin_ && !applyBTagSF_) {
      binCC = RA2bin;
      hCCbins->Fill(Double_t(binCC), eventWt);
    } else {
      // Calculate binCC
      std::vector<int> jbk;
      int binKin = kinBin(HT, MHT);
      if (binKin < 0) continue;
      int binNjets = nJetThresholds_.size()-1;
      while (NJets < nJetThresholds_[binNjets]) binNjets--;
      if (!applyBTagSF_) {
	int binNb = nbThresholds_.size()-1;
	while (BTags < nbThresholds_[binNb]) binNb--;
	jbk = {binNjets, binNb, binKin};
	try {
	  binCC = toCCbin_.at(jbk);
	}
	catch (const std::out_of_range& oor) {
          // Omitted bins j = 3,4, k = 0,3
	  // std::cerr << "jpk out of range: " << oor.what() << ": j = " << binNjets
	  // 	    << ", b = " << binNb << ", k = " << binKin << ", RA2bin = " << RA2bin << '\n';
	  continue;
	}
	// if (outCount < 100) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", binCC = " << binCC << ", RA2bin = " << RA2bin << endl;
	hCCbins->Fill(Double_t(binCC), eventWt);
      } else {
        // apply BTagSF to all Nb bins
	// if (count < 20 || count % 10000 == 0) cout << "Size of input Jets = " << Jets->size() << ", Jets_hadronFlavor = " << Jets_hadronFlavor->size() << " Jets_HTMask = " << Jets_HTMask->size() << endl;
	vector<double> probNb = btagcorr->GetCorrections(Jets, Jets_hadronFlavor, Jets_HTMask);
	for (int binNb = 0; binNb < (int) nbThresholds_.size(); ++binNb) {
	  jbk = {binNjets, binNb, binKin};
	  try {binCC = toCCbin_.at(jbk);}
	  catch (const std::out_of_range& oor) {continue;}
	  if (count % 100000 == 0) cout << "j = " << binNjets << ", NbTags = " << BTags << ", b = " << binNb << ", b wt = " << probNb[binNb] << ", k = " << binKin << ", binCC = " << binCC << endl;
	  hCCbins->Fill(Double_t(binCC), eventWt*probNb[binNb]);
	}
      }  // if apply BTagSF
    }  // if useTreeCCbin
  }  // End loop over entries

  delete forNotify;
  delete baselineTF;
  delete chain->GetCurrentFile();
  if (btagcorr) delete btagcorr;

  return hCCbins;

  }  // ======================================================================================

int
RA2bZinvAnalysis::kinBin(double& ht, double& mht) {
  int theBin = -1;
  int NmhtBins = kinThresholds_.size() - 1;
  if (deltaPhi_ != TString("nominal")) {
    // ldp or hdp
    if (mht < kinThresholds_[NmhtBins][0] || ht < kinThresholds_[NmhtBins][1]) return theBin;
    if (mht < kinThresholds_[0][0]) {
      theBin += 10;
      for (unsigned j = 1; j < kinThresholds_[NmhtBins].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[NmhtBins][j] && (j == kinThresholds_[NmhtBins].size()-1 || ht < kinThresholds_[NmhtBins][j+1])) return theBin;
      }
    }
  }
  if (mht < kinThresholds_[0][0] || ht < kinThresholds_[0][1]) return theBin;
  for (int i = 0; i < NmhtBins; ++i) {
    if (mht >= kinThresholds_[i][0] && (i == NmhtBins-1 || mht < kinThresholds_[i+1][0])) {
      for (unsigned j = 1; j < kinThresholds_[i].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[i][j] && (j == kinThresholds_[i].size()-1 || ht < kinThresholds_[i][j+1])) return theBin;
      }
    } else {
      theBin += kinThresholds_[i].size() - 1;
    }
  }
  return -2;  // Outside binned area (e.g., MHT > HT), though within preselection
}  // ======================================================================================

void
RA2bZinvAnalysis::fillCutMaps() {

  sampleKeyMap_["sig"] = "sig";
  sampleKeyMap_["T1bbbb"] = "sig";
  sampleKeyMap_["T1qqqq"] = "sig";
  sampleKeyMap_["T1tttt"] = "sig";
  sampleKeyMap_["zinv"] = "sig";
  sampleKeyMap_["topW"] = "sig";
  sampleKeyMap_["topWsle"] = "sle";
  sampleKeyMap_["topWslm"] = "slm";
  sampleKeyMap_["topsle"] = "sle";
  sampleKeyMap_["topslm"] = "slm";
  sampleKeyMap_["wjetssle"] = "sle";
  sampleKeyMap_["wjetsslm"] = "slm";
  sampleKeyMap_["qcd"] = "sig";
  sampleKeyMap_["qcdIDP"] = "sig";
  sampleKeyMap_["zinvIDP"] = "sig";
  sampleKeyMap_["wjets"] = "sig";
  sampleKeyMap_["other"] = "sig";
  sampleKeyMap_["ttzvv"] = "ttz";
  sampleKeyMap_["zmm"] = "zmm";
  sampleKeyMap_["zmm20"] = "zmm";
  sampleKeyMap_["dymm"] = "zmm";
  sampleKeyMap_["ttmm"] = "zmm";
  sampleKeyMap_["ttzmm"] = "zmm";
  sampleKeyMap_["VVmm"] = "zmm";
  sampleKeyMap_["tribosonmm"] = "zmm";
  sampleKeyMap_["zee"] = "zee";
  sampleKeyMap_["zee20"] = "zee";
  sampleKeyMap_["dyee"] = "zee";
  sampleKeyMap_["ttzee"] = "zee";
  sampleKeyMap_["ttee"] = "zee";
  sampleKeyMap_["VVee"] = "zee";
  sampleKeyMap_["tribosonee"] = "zee";
  sampleKeyMap_["zll"] = "zll";
  sampleKeyMap_["zll20"] = "zll";
  sampleKeyMap_["dyll"] = "zll";
  sampleKeyMap_["ttzll"] = "zll";
  sampleKeyMap_["ttll"] = "zll";
  sampleKeyMap_["VVll"] = "zll";
  sampleKeyMap_["photon"] = "photon";
  sampleKeyMap_["photon20"] = "photon";
  sampleKeyMap_["gjets"] = "photon";
  sampleKeyMap_["gjetsold"] = "photon";
  sampleKeyMap_["ttgjets"] = "photon";
  sampleKeyMap_["gjetsqcd"] = "photonqcd";

  if (ntupleVersion_ == "V12") {
    objCutMap_["sig"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
    objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
    objCutMap_["zll"] = "((@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0) || (@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0))";
    objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && Photons[0].Pt()>=200 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["ttz"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
    objCutMap_["slm"] = "@Muons.size()==1 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";
    objCutMap_["sle"] = "@Muons.size()==0 && @Electrons.size()==1 && isoMuonTracks==0 && isoPionTracks==0";
  } else if (ntupleVersion_ == "V15") {
    objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
    objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
    objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0) || (NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0))";
    objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && Photons[0].Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
    objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
    objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0";
    objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracks==0 && isoPionTracks==0";
  }

  minDphiCutMap_["nominal"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
  minDphiCutMap_["hdp"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
  minDphiCutMap_["ldp"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";

  MHTCutMap_["nominal"] = "MHT>=300";
  MHTCutMap_["hdp"] = "MHT>=250";
  MHTCutMap_["ldp"] = "MHT>=250";

  if (ntupleVersion_ == "V12") {
  /*
    V12, 2016:
    3: HLT_Ele105_CaloIdVT_GsfTrkIdT_v  *
    4: HLT_Ele115_CaloIdVT_GsfTrkIdT_v  *
    5: HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v
    6: HLT_Ele15_IsoVVVL_PFHT350_v  *
    7: HLT_Ele15_IsoVVVL_PFHT400_v  *
    8: HLT_Ele15_IsoVVVL_PFHT600_v
    9: HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    10: HLT_Ele25_eta2p1_WPTight_Gsf_v 
    11: HLT_Ele27_WPTight_Gsf_v  *
    12: HLT_Ele27_eta2p1_WPLoose_Gsf_v  *
    13: HLT_Ele45_WPLoose_Gsf_v
    14: HLT_Ele50_IsoVVVL_PFHT400_v

    18: HLT_IsoMu24_v  *
    19: HLT_IsoTkMu22_v
    20: HLT_IsoTkMu24_v  *
    21: HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v
    22: HLT_Mu15_IsoVVVL_PFHT350_v  *
    23: HLT_Mu15_IsoVVVL_PFHT400_v  *
    24: HLT_Mu15_IsoVVVL_PFHT600_v
    25: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    26: HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
    27: HLT_Mu45_eta2p1_v
    28: HLT_Mu50_IsoVVVL_PFHT400_v
    29: HLT_Mu50_v  *

    50: HLT_Photon135_PFMET100_v
    51: HLT_Photon165_HE10_v
    52: HLT_Photon175_v  *
    53: HLT_Photon90_CaloIdL_PFHT500_v
    54: HLT_Photon90_CaloIdL_PFHT600_v
    55: HLT_TkMu50_v

    30: HLT_PFHT200_v
    31: HLT_PFHT250_v
    32: HLT_PFHT300_PFMET100_v
    33: HLT_PFHT300_PFMET110_v
    34: HLT_PFHT300_v
    35: HLT_PFHT350_v
    36: HLT_PFHT400_v
    37: HLT_PFHT475_v
    38: HLT_PFHT600_v
    39: HLT_PFHT650_v
    40: HLT_PFHT800_v
    41: HLT_PFHT900_v
    42: HLT_PFMET100_PFMHT100_IDTight_v  *
    43: HLT_PFMET110_PFMHT110_IDTight_v  *
    44: HLT_PFMET120_PFMHT120_IDTight_v  *
    45: HLT_PFMET90_PFMHT90_IDTight_v
    46: HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v  *
    47: HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v  *
    48: HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v  *
    49: HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v
  */
    triggerMap_["zmm"] = {"18", "20", "22", "23", "29"};
    triggerMap_["zee"] = {"3", "4", "6", "7", "11", "12"};
    triggerMap_["photon"] = {"52"};  // re-miniAOD; 51 for ReReco/PromptReco
    triggerMap_["sig"] = {"42", "43", "44", "46", "47", "48"};
  } else if (ntupleVersion_ == "V15") {
  /*
    V15, 2018
    17: HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v
    18: HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v
    ...
    21: HLT_Ele15_IsoVVVL_PFHT350_v  *
    22: HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v
    23: HLT_Ele15_IsoVVVL_PFHT400_v  *
    24: HLT_Ele15_IsoVVVL_PFHT450_v
    25: HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v
    26: HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v
    27: HLT_Ele15_IsoVVVL_PFHT600_v
    28: HLT_Ele20_eta2p1_WPLoose_Gsf_v  *
    29: HLT_Ele20_WPLoose_Gsf_v
    30: HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    31: HLT_Ele25_eta2p1_WPTight_Gsf_v
    32: HLT_Ele27_WPTight_Gsf_v  x
    33: HLT_Ele27_eta2p1_WPLoose_Gsf_v  x
    34: HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v
    35: HLT_Ele32_WPTight_Gsf_v  *
    36: HLT_Ele35_WPTight_Gsf_v
    37: HLT_Ele45_WPLoose_Gsf_v
    38: HLT_Ele50_IsoVVVL_PFHT400_v
    39: HLT_Ele50_IsoVVVL_PFHT450_v
    40: HLT_Ele105_CaloIdVT_GsfTrkIdT_v  *
    41: HLT_Ele115_CaloIdVT_GsfTrkIdT_v  *
    42: HLT_Ele135_CaloIdVT_GsfTrkIdT_v
    43: HLT_Ele145_CaloIdVT_GsfTrkIdT_v

    19: HLT_DoubleMu8_Mass8_PFHT300_v
    20: HLT_DoubleMu8_Mass8_PFHT350_v
    ...
    44: HLT_IsoMu16_eta2p1_MET30_v
    45: HLT_IsoMu20_v
    46: HLT_IsoMu22_v
    47: HLT_IsoMu22_eta2p1_v
    48: HLT_IsoMu24_v  *  prescaled in late 2017 --Owen
    49: HLT_IsoMu24_eta2p1_v
    50: HLT_IsoMu27_v
    51: HLT_IsoTkMu22_v
    52: HLT_IsoTkMu24_v  *
    53: HLT_Mu15_IsoVVVL_PFHT350_v  *
    54: HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v
    55: HLT_Mu15_IsoVVVL_PFHT400_v  *
    56: HLT_Mu15_IsoVVVL_PFHT450_v
    57: HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v
    58: HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v
    59: HLT_Mu15_IsoVVVL_PFHT600_v
    60: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    61: HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
    62: HLT_Mu45_eta2p1_v
    63: HLT_Mu50_v  *
    64: HLT_Mu50_IsoVVVL_PFHT400_v
    65: HLT_Mu50_IsoVVVL_PFHT450_v
    66: HLT_Mu55_v

    135: HLT_Photon90_CaloIdL_PFHT500_v
    136: HLT_Photon90_CaloIdL_PFHT600_v
    137: HLT_Photon90_CaloIdL_PFHT700_v
    138: HLT_Photon135_PFMET100_v
    139: HLT_Photon165_HE10_v
    140: HLT_Photon165_R9Id90_HE10_IsoM_v
    141: HLT_Photon175_v  *
    142: HLT_Photon200_v
    143: HLT_Photon300_NoHE_v

    67: HLT_PFHT200_v
    68: HLT_PFHT250_v
    69: HLT_PFHT300_v
    70: HLT_PFHT300_PFMET100_v
    71: HLT_PFHT300_PFMET110_v
    72: HLT_PFHT350_v
    73: HLT_PFHT370_v
    74: HLT_PFHT380_SixPFJet32_v
    75: HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v
    76: HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v
    77: HLT_PFHT400_v
    78: HLT_PFHT400_SixJet30_v
    79: HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v
    80: HLT_PFHT430_v
    81: HLT_PFHT430_SixJet40_BTagCSV_p056_v
    82: HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v
    83: HLT_PFHT430_SixPFJet40_v
    84: HLT_PFHT450_SixJet40_v
    85: HLT_PFHT450_SixJet40_BTagCSV_p056_v
    86: HLT_PFHT450_SixPFJet40_PFBTagCSV_1p5_v
    87: HLT_PFHT475_v
    88: HLT_PFHT500_PFMET100_PFMHT100_IDTight_v
    89: HLT_PFHT500_PFMET110_PFMHT110_IDTight_v
    90: HLT_PFHT510_v
    91: HLT_PFHT590_v
    92: HLT_PFHT600_v
    93: HLT_PFHT650_v
    94: HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v
    95: HLT_PFHT680_v
    96: HLT_PFHT700_PFMET85_PFMHT85_IDTight_v
    97: HLT_PFHT700_PFMET95_PFMHT95_IDTight_v
    98: HLT_PFHT780_v
    99: HLT_PFHT800_v
    100: HLT_PFHT800_PFMET75_PFMHT75_IDTight_v
    101: HLT_PFHT800_PFMET85_PFMHT85_IDTight_v
    102: HLT_PFHT890_v
    103: HLT_PFHT900_v
    104: HLT_PFHT1050_v
    105: HLT_PFJet500_v
    106: HLT_PFJet550_v
    107: HLT_PFMET90_PFMHT90_IDTight_v
    108: HLT_PFMET100_PFMHT100_IDTight_v  *
    109: HLT_PFMET100_PFMHT100_IDTight_PFHT60_v
    110: HLT_PFMET110_PFMHT110_IDTight_v  x
    111: HLT_PFMET110_PFMHT110_IDTight_PFHT60_v
    112: HLT_PFMET120_PFMHT120_IDTight_v  *
    113: HLT_PFMET120_PFMHT120_IDTight_PFHT60_v
    114: HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v
    115: HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v
    116: HLT_PFMET130_PFMHT130_IDTight_v
    117: HLT_PFMET130_PFMHT130_IDTight_PFHT60_v
    118: HLT_PFMET140_PFMHT140_IDTight_v
    119: HLT_PFMET140_PFMHT140_IDTight_PFHT60_v
    120: HLT_PFMET500_PFMHT500_IDTight_CalBTagCSV_3p1_v
    121: HLT_PFMET700_PFMHT700_IDTight_CalBTagCSV_3p1_v
    122: HLT_PFMET800_PFMHT800_IDTight_CalBTagCSV_3p1_v
    123: HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v
    124: HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v  *
    125: HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v
    126: HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v  x
    127: HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v
    128: HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v  *
    129: HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v
    130: HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v
    131: HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v
    132: HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v
    133: HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v
    134: HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60_v
   */
    triggerMap_["zmm"] = {"48", "52", "53", "55", "63"};
    triggerMap_["zee"] = {"21", "23", "28", "35", "40", "41"};
    triggerMap_["photon"] = {"141"};
    triggerMap_["sig"] = {"108", "112", "124", "128"};
  }

  triggerMap_["zll"].reserve(triggerMap_["zmm"].size() + triggerMap_["zee"].size());
  triggerMap_["zll"] = triggerMap_["zmm"];
  triggerMap_["zll"].insert(triggerMap_["zll"].end(), triggerMap_["zee"].begin(), triggerMap_["zee"].end());
  triggerMap_["sle"] = triggerMap_["sig"];
  triggerMap_["slm"] = triggerMap_["sig"];

}  // ======================================================================================

// void
// RA2bZinvAnalysis::fillCutFlow(TH1F* hcf, Double_t wt) {
//   hcf->Fill(0.0, wt);
//   if (HTcut_) {
//     hcf->Fill(1.0, wt);
//     if (MHTcut_) {
//       hcf->Fill(2.0, wt);
//     }
//   }
// }

RA2bZinvAnalysis::cutHistos::cutHistos(TChain* chain, TObjArray* forNotify) : forNotify_(forNotify) {
  HTcutf_ = new TTreeFormula("HTcut", HTcut_, chain);  forNotify->Add(HTcutf_);
  MHTcutf_ = new TTreeFormula("MHTcut", MHTcut_, chain);  forNotify->Add(MHTcutf_);
  NJetscutf_ = new TTreeFormula("NJetscut", NJetscut_, chain);  forNotify->Add(NJetscutf_);
  minDphicutf_ = new TTreeFormula("minDphicut", minDphicut_, chain);  forNotify->Add(minDphicutf_);
  objcutf_ = new TTreeFormula("objcut", objcut_, chain);  forNotify->Add(objcutf_);
  ptcutf_ = new TTreeFormula("ptcut", ptCut_, chain);  forNotify->Add(ptcutf_);
  masscutf_ = new TTreeFormula("masscut", massCut_, chain);  forNotify->Add(masscutf_);
  trigcutf_ = new TTreeFormula("trigcut", trigCuts_, chain);  forNotify->Add(trigcutf_);
  commoncutf_ = new TTreeFormula("commoncut", commonCuts_, chain);  forNotify->Add(commoncutf_);
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1F* hcf, Double_t wt) {
  hcf->GetXaxis()->SetBinLabel(1, "None");
  hcf->GetXaxis()->SetBinLabel(2, "HT");
  hcf->GetXaxis()->SetBinLabel(3, "MHT");
  hcf->GetXaxis()->SetBinLabel(4, "NJets");
  hcf->GetXaxis()->SetBinLabel(5, "mnDphi");
  hcf->GetXaxis()->SetBinLabel(6, "objects");
  hcf->GetXaxis()->SetBinLabel(7, "Zpt");
  hcf->GetXaxis()->SetBinLabel(8, "Zmass");
  hcf->GetXaxis()->SetBinLabel(9, "Trigger");
  hcf->GetXaxis()->SetBinLabel(10, "Filters");
  hcf->GetXaxis()->LabelsOption("vu");

  HTcutf_->GetNdata();
  MHTcutf_->GetNdata();
  NJetscutf_->GetNdata();
  minDphicutf_->GetNdata();
  objcutf_->GetNdata();
  ptcutf_->GetNdata();
  masscutf_->GetNdata();
  trigcutf_->GetNdata();
  commoncutf_->GetNdata();

  hcf->Fill(0.5, wt);
  if (TString(hcf->GetName()).Contains(TString("Flow"))) {
    if (HTcutf_->EvalInstance(0)) {
      hcf->Fill(1.5, wt);
      if (MHTcutf_->EvalInstance(0)) {
	hcf->Fill(2.5, wt);
	if (NJetscutf_->EvalInstance(0)) {
	  hcf->Fill(3.5, wt);
	  if (minDphicutf_->EvalInstance(0)) {
	    hcf->Fill(4.5, wt);
	    if (objcutf_->EvalInstance(0)) {
	      hcf->Fill(5.5, wt);
	      if (ptcutf_->EvalInstance(0)) {
		hcf->Fill(6.5, wt);
		if (masscutf_->EvalInstance(0)) {
		  hcf->Fill(7.5, wt);
		  if (trigcutf_->EvalInstance(0)) {
		    hcf->Fill(8.5, wt);
		    if (commoncutf_->EvalInstance(0)) {
		      hcf->Fill(9.5, wt);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } else {
    if (HTcutf_->EvalInstance(0)) hcf->Fill(1.5, wt);
    if (MHTcutf_->EvalInstance(0)) hcf->Fill(2.5, wt);
    if (NJetscutf_->EvalInstance(0)) hcf->Fill(3.5, wt);
    if (minDphicutf_->EvalInstance(0)) hcf->Fill(4.5, wt);
    if (objcutf_->EvalInstance(0)) hcf->Fill(5.5, wt);
    if (ptcutf_->EvalInstance(0)) hcf->Fill(6.5, wt);
    if (masscutf_->EvalInstance(0)) hcf->Fill(7.5, wt);
    if (trigcutf_->EvalInstance(0)) hcf->Fill(8.5, wt);
    if (commoncutf_->EvalInstance(0)) hcf->Fill(9.5, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::runMakeClass(const char* sample, const char* ext) {
  TChain* chain = getChain(sample, nullptr, false);
  TString templateName("TreeMkrTemplate_");
  templateName += ext;
  chain->MakeClass(templateName.Data());

}  // ======================================================================================

void
RA2bZinvAnalysis::checkTrigPrescales(const char* sample) {
  Int_t fCurrent;  // current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);
  Long64_t Nentries = chain->GetEntries();
  Int_t countInFile = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) cout << "Current file in chain: " << thisFile->GetName() << endl;
      countInFile = 0;
    }
    chain->GetEntry(entry);
    countInFile++;
    if (countInFile == 1) {
      int trigNo = 0;
      for (auto & theTrigPrescale : *TriggerPrescales) {
	if (theTrigPrescale != 1) cout << trigNo << ":  " << theTrigPrescale << endl;
    	++trigNo;
      }
    }
  }

}  // ======================================================================================
