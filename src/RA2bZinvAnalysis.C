//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#ifndef RA2BZINVANALYSIS_H
#define RA2BZINVANALYSIS_H

#define ISMC
// #undef ISMC

#include <TChain.h>
#include <TTreeReaderValue.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTreeFormula.h>
#include <TChainElement.h>
#ifdef ISMC
#include "../../Analysis/btag/BTagCorrector.h"
#endif

class RA2bZinvAnalysis {

public:
  RA2bZinvAnalysis();

  virtual ~RA2bZinvAnalysis() {};

  TChain* getChain(const char* sample, Int_t* fCurrent = nullptr, bool setBrAddr = true);
  std::vector<TH1F*> makeHistograms(const char* sample);
  TH1F* makeCChist(const char* sample);
  TCut getCuts(const TString sampleKey);
  int kinBin(double& ht, double& mht);
  void runMakeClass(const char* sample, const char* ext);

  struct hist1D {
    TH1F* hist;
    const char* name;
    const char* title;
    std::pair<const char*, const char*> axisTitles;
    Int_t Nbins;
    Double_t lowEdge;
    Double_t highEdge;
    Double_t* dvalue;
    Int_t* ivalue;
    void (RA2bZinvAnalysis::*filler)(TH1F* h, double wt);
    std::vector<TString*> omitCut;
    TString NminusOneCuts;
    TTreeFormula* NminusOneFormula;
    hist1D() : dvalue(nullptr), ivalue(nullptr), filler(nullptr) {}
  };

private:
  const char* treeLoc_;
  const char* treeName_;
  TString era_;  // "2016", ...
  TString ntupleVersion_; // "V12", ...
  TString deltaPhi_;  // "nominal", "hdp", "ldp"
  bool applyMassCut_;
  bool applyPtCut_;
  bool applyMinDeltaRCut_;
  // bool applySF_;
  // bool njSplit_;
  bool useTreeCCbin_;
#ifdef ISMC
  bool applyBTagSF_;
  bool applyPuWeight_;
  bool customPuWeight_;
  TH1* puHist_;
  const char* BTagSFfile_;
#endif
  std::vector< std::vector<double> > kinThresholds_;
  std::vector<int> nJetThresholds_;
  std::vector<int> nbThresholds_;
  unsigned kinSize_;
  double intLumi_;

  typedef std::map<TString, std::vector<TString> > vstring_map;
  typedef std::map<TString, TString> string_map;
  typedef std::map<std::vector<int>, Int_t> ivector_map;
  vstring_map fileMap_;
  vstring_map triggerMap_;
  string_map objCutMap_;
  string_map minDphiCutMap_;
  string_map MHTCutMap_;
  string_map sampleKeyMap_;
  ivector_map toCCbin_;
  TString HTcut_;
  TString MHTcut_;
  TString NJetscut_;
  TString massCut_;
  TString ptCut_;

  void fillFileMap();
  void fillCutMaps();
  void bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms);

#ifdef ISMC
#include "LeafDeclaration_MC_V12.h"
#else
#include "LeafDeclaration_data_V12.h"
#endif

  // Functions to fill histograms with non-double, non-int types
  void fillZmass(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.M(), wt);}
  void fillZpt(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.Pt(), wt);}

  ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID
};

// ======================================================================================

#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
// #include <TLegend.h>
// #include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TString.h>
#include <TRegexp.h>
#include <TCut.h>
#include <TMath.h>
using TMath::Sqrt; using TMath::Power;

#define _CRT_SECURE_NO_WARNINGS

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
  useTreeCCbin_(true)
#ifdef ISMC
  , applyBTagSF_(false),
  applyPuWeight_(false),
  customPuWeight_(true)  // Substitute Kevin P recipe for the PuWeight in the tree
#endif
{
  if (era_ == TString("2016")) {
    if (ntupleVersion_ == "V12") {
      // treeLoc_ = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12";  // Colorado, owned by wtford (Zjets only)
      treeLoc_ = "/nfs/data38/cms/mulholland/lpcTrees/Skims/Run2ProductionV12";  // Colorado, owned by mulholland
      // treeLoc_ = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV12";  // xrootd
      // treeLoc_ = "/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV12";  // from cmslpc
    }
    intLumi_ = 35.9;

#ifdef ISMC
    if (applyPuWeight_ && customPuWeight_) {
      TFile* pufile = TFile::Open("../../Analysis/corrections/PileupHistograms_0121_69p2mb_pm4p6.root","READ");
      puHist_ = (TH1*) pufile->Get("pu_weights_down");
    }

    BTagSFfile_ = "../../Analysis/btag/CSVv2_Moriond17_B_H_mod.csv";
#endif

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

#include "fillFileMap_V12.h"

TChain*
RA2bZinvAnalysis::getChain(const char* sample, Int_t* fCurrent, bool setBrAddr) {
  TString theSample(sample);
  TString key;
  if (theSample.Contains("zinv")) key = TString("zinv");
  else if (theSample.Contains("ttzvv")) key = TString("ttzvv");
  else if (theSample.Contains("dymm")) key = TString("dymm");
  else if (theSample.Contains("dyee")) key = TString("dyee");
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

  if (fCurrent != nullptr) *fCurrent = -1;
  if (setBrAddr) setBranchAddress(chain);

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
  TString commonCuts;
  TString trigCuts;
  std::vector<TString> trigger;
  try {trigger = triggerMap_.at(sample);}
  catch (const std::out_of_range& oor) {trigger.clear();}

  // commonCuts = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  if (trigger.empty()) {
    commonCuts = "JetID==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-
  } else {
    commonCuts = "JetID==1 && globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-
    int Ntrig = trigger.size();
    if (Ntrig > 1) trigCuts += TString("(");
    for (auto theTrigger : trigger)
      trigCuts += TString("(TriggerPass[")+theTrigger+TString("]==1) + ");
    if (Ntrig > 1) trigCuts.Replace(trigCuts.Length()-3, 3, ")");
  }
  cout << "trigCuts = " << trigCuts << endl;

  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "Photons[0].Pt()>=200.";
    if (trigger.empty() && applyMinDeltaRCut_) photonDeltaRcut = "madMinPhotonDeltaR>=0.4";
  }
				       
    // 	if(extraCuts!=None):
    //     cuts+=extraCuts

  HTcut_ = std::string("HT>=") + std::to_string(kinThresholds_[0][1]);
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  NJetscut_ = std::string("NJets>=") + std::to_string(nJetThresholds_[0]);

  cuts += objCutMap_.at(sampleKey);
  cuts += HTcut_;
  cuts += NJetscut_;
  cuts += MHTcut_;
  cuts += minDphiCutMap_.at(deltaPhi_);
  cuts += massCut_;
  cuts += ptCut_;
  cuts += photonDeltaRcut;
  cuts += commonCuts;
  cuts += trigCuts;
    // 	  if(applySF):
    //     cuts*=bJetCutsSF[bJetBin]
    // 	else:
    //     cuts+=bJetCuts[bJetBin]
#ifdef ISMC
  if(applyPuWeight_ && !customPuWeight_) cuts *= "puWeight*(1)";
    // 	  if(type(extraWeight) is str):
    //     extraWeight+="*(1)"
    //     cuts*=extraWeight
#endif

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

    Double_t eventWt = 1;
    Double_t PUweight = 1;
#ifdef ISMC
    if (applyPuWeight_ && customPuWeight_) {
      // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    eventWt = 1000*intLumi_*Weight*PUweight;
    if (eventWt < 0) eventWt *= -1;
#endif

    for (auto & hg : histograms) {
      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt != 0) {
	if (hg->dvalue != nullptr) hg->hist->Fill(*(hg->dvalue), selWt*eventWt);
	else if (hg->ivalue != nullptr) hg->hist->Fill(Double_t(*(hg->ivalue)), selWt*eventWt);
	else if (hg->filler != nullptr) (this->*(hg->filler))(hg->hist, selWt*eventWt);
	else cerr << "No method to fill histogram provided for " << hg->name << endl;
      }
    }  // loop over histograms
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

#ifdef ISMC
  BTagCorrector* btagcorr = nullptr;
  if (applyBTagSF_) {
    btagcorr = new BTagCorrector;
    btagcorr->SetCalib(BTagSFfile_);
  }
#endif

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
#ifdef ISMC
    	if (btagcorr) btagcorr->SetEffs(thisFile);
#endif
      }
    }
    chain->GetEntry(entry);

    UInt_t binCC = 0;

    // Apply baseline selection
    baselineTF->GetNdata();
    double selWt = baselineTF->EvalInstance(0);
    if (selWt == 0) continue;
    outCount++;

    // Compute event weight factors
    Double_t eventWt = 1;
    Double_t PUweight = 1;
#ifdef ISMC
    if (applyPuWeight_ && customPuWeight_) {
      // This recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    if (count < 20 || count % 10000 == 0) cout << "PUweight = " << PUweight << endl;
    eventWt = 1000*intLumi_*Weight*selWt*PUweight;
    if (eventWt < 0) eventWt *= -1;
    if (count < 20 || count % 10000 == 0) cout << "eventWt = " << eventWt << endl;

    if (useTreeCCbin_ && !applyBTagSF_) {
#else
    if (useTreeCCbin_) {
#endif
      binCC = RA2bin;
      hCCbins->Fill(Double_t(binCC), eventWt);
    } else {
      std::vector<int> jbk;
      int binKin = kinBin(HT, MHT);
      if (binKin < 0) continue;
      int binNjets = nJetThresholds_.size()-1;
      while (NJets < nJetThresholds_[binNjets]) binNjets--;
#ifdef ISMC
      if (!applyBTagSF_) {
#endif
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
#ifdef ISMC
      } else {  // apply BTagSF to all Nb bins
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
#endif
    }  // if useTreeCCbin
  }  // End loop over entries

  delete forNotify;
  delete baselineTF;
  delete chain->GetCurrentFile();
#ifdef ISMC
  if (btagcorr) delete btagcorr;
#endif

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
  sampleKeyMap_["dibosonmm"] = "zmm";
  sampleKeyMap_["tribosonmm"] = "zmm";
  sampleKeyMap_["zee"] = "zee";
  sampleKeyMap_["zee20"] = "zee";
  sampleKeyMap_["dyee"] = "zee";
  sampleKeyMap_["ttzee"] = "zee";
  sampleKeyMap_["ttee"] = "zee";
  sampleKeyMap_["dibosonee"] = "zee";
  sampleKeyMap_["tribosonee"] = "zee";
  sampleKeyMap_["zll"] = "zll";
  sampleKeyMap_["zll20"] = "zll";
  sampleKeyMap_["dyll"] = "zll";
  sampleKeyMap_["ttzll"] = "zll";
  sampleKeyMap_["ttll"] = "zll";
  sampleKeyMap_["dibosonll"] = "zll";
  sampleKeyMap_["photon"] = "photon";
  sampleKeyMap_["photon20"] = "photon";
  sampleKeyMap_["gjets"] = "photon";
  sampleKeyMap_["gjetsold"] = "photon";
  sampleKeyMap_["ttgjets"] = "photon";
  sampleKeyMap_["gjetsqcd"] = "photonqcd";

  objCutMap_["sig"] = 
    "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
  objCutMap_["zmm"] = 
    "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
  objCutMap_["zee"] = 
    "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
  objCutMap_["zll"] = 
    "((@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0) || (@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0))";
  objCutMap_["photon"] = 
    "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
  objCutMap_["photonqcd"] = 
    "Sum$(Photons_nonPrompt)!=0 && Photons[0].Pt()>=200 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
  objCutMap_["ttz"] = 
    "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
  objCutMap_["slm"] = 
    "@Muons.size()==1 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";
  objCutMap_["sle"] = 
    "@Muons.size()==0 && @Electrons.size()==1 && isoMuonTracks==0 && isoPionTracks==0";

  minDphiCutMap_["nominal"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
  minDphiCutMap_["hdp"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
  minDphiCutMap_["ldp"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";

  MHTCutMap_["nominal"] = "MHT>=300";
  MHTCutMap_["hdp"] = "MHT>=250";
  MHTCutMap_["ldp"] = "MHT>=250";

  triggerMap_["zmm"] = {"22", "23", "29", "18", "20"};
  triggerMap_["zee"] = {"6", "7", "11", "12", "3", "4"};
  triggerMap_["zll"].reserve(triggerMap_["zmm"].size() + triggerMap_["zee"].size());
  triggerMap_["zll"] = triggerMap_["zmm"];
  triggerMap_["zll"].insert(triggerMap_["zll"].end(), triggerMap_["zee"].begin(), triggerMap_["zee"].end());
  triggerMap_["photon"] = {"52"};  // re-miniAOD; 51 for ReReco/PromptReco
  triggerMap_["sig"] = {"42", "43", "44", "46", "47", "48"};
  triggerMap_["sle"] = triggerMap_["sig"];
  triggerMap_["slm"] = triggerMap_["sig"];

}  // ======================================================================================

void
RA2bZinvAnalysis::runMakeClass(const char* sample, const char* ext) {
  TChain* chain = getChain(sample, nullptr, false);
  TString templateName("TreeMkrTemplate_");
  templateName += ext;
  chain->MakeClass(templateName.Data());

}  // ======================================================================================

#endif
