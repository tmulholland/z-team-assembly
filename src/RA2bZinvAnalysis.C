//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//
//  (see TreeMaker/TreeMaker/python/makeTreeFromMiniAOD_cff.py
//   for tree branch names)
//


#include "TChain.h"
#include <TTreeReaderValue.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTreeFormula.h>

class RA2bZinvAnalysis {

public:
  RA2bZinvAnalysis();

  virtual ~RA2bZinvAnalysis() {};

  TChain* getChain(const char* sample);
  std::vector<TH1F*> makeHistograms(const char* sample);
  TH1F* makeCChist(const char* sample);
  TCut getCuts(const TString sampleKey);
  int kinBin(double& ht, double& mht);
  struct hist1D {
    TH1F* hist;
    const char* name;
    const char* title;
    Int_t Nbins;
    Double_t lowEdge;
    Double_t highEdge;
    Double_t* dvalue;
    Int_t* ivalue;
    TString* omitCut;
    TString NminusOneCuts;
    TTreeFormula* NminusOneFormula;
  };

private:
  const char* treeLoc_;
  const char* treeName_;
  TString era_;  // "2016", ...
  TString ntupleVersion_; // "V12", ...
  TString deltaPhi_;  // "nominal", "hdp", "ldp"
  bool applyMassCut_;
  bool applyPtCut_;
  bool applySF_;
  bool njSplit_;
  bool useTreeCCbin_;
  bool applyPuWeight_;
  bool customPuWeight_;
  TH1* puHist_;
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
  string_map bTagSF_;
  ivector_map toCCbin_;
  TString HTcut_;
  TString MHTcut_;
  TString NJetscut_;

  void fillFileMap();
  void fillCutMaps();
  void bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms);

#include "LeafDeclaration_V12.h"

  ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID
};

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
  applySF_(false),
  njSplit_(false),
  useTreeCCbin_(true),
  applyPuWeight_(true),
  customPuWeight_(true)
{
  if (era_ == TString("2016")) {
    if (ntupleVersion_ == "V12") {
      // treeLoc_ = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12";
      treeLoc_ = "/nfs/data38/cms/mulholland/lpcTrees/Skims/Run2ProductionV12";
    }
    intLumi_ = 35.9;

    if (applyPuWeight_ && customPuWeight_) {
      TFile* pufile = TFile::Open("../plots/histograms/PileupHistograms_0121_69p2mb_pm4p6.root","READ");
      puHist_ = (TH1*) pufile->Get("pu_weights_down");
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

}

#include "fillFileMap_V12.h"

TChain*
RA2bZinvAnalysis::getChain(const char* sample) {
  TString theSample(sample);
  TString key;
  if (theSample.Contains("zinv")) key = TString("zinv");
  else if (theSample.Contains("ttzvv")) key = TString("ttzvv");

  TChain* chain = new TChain(treeName_);
  for (auto file : fileMap_.at(key)) {
    cout << file << endl;
    chain->Add(file);
  }

  setBranchAddress(chain);

  return chain;
}

TCut
RA2bZinvAnalysis::getCuts(const TString sample) {

  TCut cuts;
  TString sampleKey;
  try {sampleKey = sampleKeyMap_.at(sample);}
  catch (const std::out_of_range& oor) {
    std::cerr << oor.what() << " getCuts called with invalid sample = " << sample << endl;
    return cuts;
  }

  TString massCut;
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut = "ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  TString ptCut;
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut = "ZCandidates.Pt()>=200.";

  if ((sampleKey == "photon") && applyPtCut_)
    ptCut = "Photons[0].Pt()>=200.";

  TString photonDeltaRcut;
  TString commonCuts;
  TString trigCuts;
  std::vector<TString> trigger;
  try {trigger = triggerMap_.at(sample);}
  catch (const std::out_of_range& oor) {trigger.clear();}
  if (trigger.empty()) {
    commonCuts = "JetID==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";
    if (sampleKey=="photon") photonDeltaRcut = "madMinPhotonDeltaR>=0.4";
  } else {
    commonCuts = "JetID==1 && globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";
    int Ntrig = trigger.size();
    if (Ntrig > 1) trigCuts += TString("(");
    for (auto theTrigger : trigger)
      trigCuts += TString("(TriggerPass[")+theTrigger+TString("]==1) + ");
    if (Ntrig > 1) trigCuts.Replace(trigCuts.Length()-3, 3, ")");
  }
  cout << "trigCuts = " << trigCuts << endl;
				       
    // 	if(extraCuts!=None):
    //     cuts+=extraCuts

  HTcut_ = std::string("HT>=") + std::to_string(kinThresholds_[0][1]);
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  NJetscut_ = std::string("NJets>=") + std::to_string(nJetThresholds_[0]);
  cout << "HTcut_ = " << HTcut_ << ", NJetscut_ = " << NJetscut_ << endl;

  cuts += objCutMap_.at(sampleKey);
  cuts += HTcut_;
  cuts += NJetscut_;
  cuts += MHTcut_;
  cuts += minDphiCutMap_.at(deltaPhi_);
  cuts += massCut;
  cuts += ptCut;
  cuts += photonDeltaRcut;
  cuts += commonCuts;
  cuts += trigCuts;
    // 	  if(applySF):
    //     cuts*=bJetCutsSF[bJetBin]
    // 	else:
    //     cuts+=bJetCuts[bJetBin]
  if(applyPuWeight_ && !customPuWeight_) cuts *= "puWeight*(1)";
    // 	  if(type(extraWeight) is str):
    //     extraWeight+="*(1)"
    //     cuts*=extraWeight

  return cuts;
 
}

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms) {
  //
  // Define N - 1 cuts, book histograms.  Traverse the chain and fill.
  //
  TCut baselineCuts = getCuts(sample);
  cout << "baseline = " << endl << baselineCuts << endl;
  TChain* chain = getChain(sample);
  TObjArray* forNotify = new TObjArray;

  for (auto & hg : histograms) {
    hg->hist = new TH1F(hg->name, hg->title, hg->Nbins, hg->lowEdge, hg->highEdge);
    hg->hist->SetOption("HIST");
    hg->hist->SetMarkerSize(0);
    hg->NminusOneCuts = baselineCuts;
    if (hg->omitCut->Length() > 0) hg->NminusOneCuts(*(hg->omitCut)) = "1";
    cout << "For sample " << sample << ", histo " << hg->name << ", hg->omitCut = " << *(hg->omitCut) << ", cuts = " << endl
	 << hg->NminusOneCuts << endl;
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
    chain->GetEntry(entry);

    Double_t puWeight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      puWeight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    Double_t eventWt = 1000*intLumi_*Weight*puWeight;
    if (eventWt < 0) eventWt *= -1;

    for (auto & hg : histograms) {
      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt != 0) {
	Double_t x = 0;
	if (hg->dvalue != nullptr) x = *(hg->dvalue);
	else if (hg->ivalue != nullptr) x = *(hg->ivalue);
	hg->hist->Fill(x, selWt*eventWt);
      }
    }
  }
}

std::vector<TH1F*>
RA2bZinvAnalysis::makeHistograms(const char* sample) {
  //
  // Define histograms, variables to fill, and cuts to be modified.
  // Return a vector of histograms.
  //
  std::vector<hist1D*> histograms;

  TString nullstring("");

  hist1D hHT;
  hHT.name = TString("hHT_") + TString(sample);  hHT.title = "HT";
  hHT.Nbins = 60;  hHT.lowEdge = 0;  hHT.highEdge = 3000; // 
  hHT.dvalue = &HT;  hHT.ivalue = nullptr;  hHT.omitCut = &HTcut_;
  histograms.push_back(&hHT);

  hist1D hMHT;
  hMHT.name = TString("hMHT_") + TString(sample);  hMHT.title = "MHT";
  hMHT.Nbins = 60;  hMHT.lowEdge = 0;  hMHT.highEdge = 3000;
  hMHT.dvalue = &MHT;  hMHT.ivalue = nullptr;  hMHT.omitCut = &MHTcut_;
  histograms.push_back(&hMHT);

  hist1D hNJets;
  hNJets.name = TString("hNJets_") + TString(sample);  hNJets.title = "NJets";
  hNJets.Nbins = 20;  hNJets.lowEdge = 0;  hNJets.highEdge = 20;
  hNJets.ivalue = &NJets;  hNJets.dvalue = nullptr;  hNJets.omitCut = &NJetscut_;
  histograms.push_back(&hNJets);

  hist1D hBTags;
  hBTags.name = TString("hBTags_") + TString(sample);  hBTags.title = "BTags";
  hBTags.Nbins = 20;  hBTags.lowEdge = 0;  hBTags.highEdge = 20;
  hBTags.ivalue = &BTags;  hBTags.dvalue = nullptr;  hBTags.omitCut = &nullstring;
  histograms.push_back(&hBTags);

  bookAndFillHistograms(sample, histograms);

  std::vector<TH1F*> theHists;
  theHists.push_back(hHT.hist);
  theHists.push_back(hMHT.hist);
  theHists.push_back(hNJets.hist);
  theHists.push_back(hBTags.hist);

  return theHists;

}

TH1F*
RA2bZinvAnalysis::makeCChist(const char* sample) {
  //
  // Book the analysis-binned histogram; traverse the chain, fill, and return it.
  //
  Int_t MaxBins = toCCbin_.size();
  cout << "MaxBins = " << MaxBins << ", applyPuWeight = " << applyPuWeight_ << endl;

  TH1F* hCCbins = new TH1F("hCCbins", "Zinv background estimate", MaxBins, 0.5, MaxBins+0.5);
  hCCbins->SetOption("HIST");
  hCCbins->SetMarkerSize(0);
  
  TChain* chain = getChain(sample);
  // chain->Print();
  // chain->MakeCode();

  TCut baselineCuts = getCuts(sample);
  cout << "baseline = " << baselineCuts << endl;

  Long64_t Nentries = chain->GetEntries();
  cout << "Nentries in tree = " << Nentries << endl;
  int count = 0, outCount = 0;
  // See https://root-forum.cern.ch/t/how-to-evaluate-a-formula-for-a-single-tree-event/12283
  TTreeFormula* select = new TTreeFormula("select", baselineCuts, chain);
  chain->SetNotify(select);  // This is needed only for TChain.
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (count % 100000 == 0) cout << "Entry number " << count << endl;

    chain->LoadTree(entry);

    select->GetNdata();
    double selWt = select->EvalInstance(0);
    if (selWt == 0) continue;
    outCount++;
    // if (outCount < 100) cout << "selWt = " << selWt << endl;

    chain->GetEntry(entry);

    UInt_t binCC = 0;
    if (useTreeCCbin_) {
      binCC = RA2binBranch;
    } else {
      int binKin = kinBin(HT, MHT);
      if (binKin < 0) continue;
      int binNjets = nJetThresholds_.size()-1;
      while (NJets < nJetThresholds_[binNjets]) binNjets--;
      int binNb = nbThresholds_.size()-1;
      while (BTags < nbThresholds_[binNb]) binNb--;

      std::vector<int> jbk = {binNjets, binNb, binKin};
      try {
    	binCC = toCCbin_.at(jbk);
      }
      catch (const std::out_of_range& oor) {
    	// Omitted bins j = 3,4, k = 0,3
    	// std::cerr << "jpk out of range: " << oor.what()
    	// 	  << ": j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", RA2binBranch = " << RA2binBranch << '\n';
    	continue;
      }
      if (outCount < 100) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", binCC = " << binCC << ", RA2binBranch = " << RA2binBranch << ", Weight = " << Weight << endl;
    }

    Double_t puWeight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      puWeight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    Double_t eventWt = 1000*intLumi_*Weight*selWt*puWeight;
    if (eventWt < 0) eventWt *= -1;
    hCCbins->Fill(Double_t(binCC), eventWt);

  }  // End loop over entries

  return hCCbins;

}

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
}

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

  bTagSF_["0"] = "BTagsSF[0]*(1)";
  bTagSF_["1"] = "BTagsSF[1]*(1)";
  bTagSF_["2"] = "BTagsSF[2]*(1)";
  bTagSF_["3"] = "BTagsSF[3]*(1)";
  bTagSF_["none"] = "(1)";

  triggerMap_["zmm"] = {"22", "23", "29", "17"};
  triggerMap_["zee"] = {"6", "7", "12", "3"};
  triggerMap_["zll"].reserve(triggerMap_["zmm"].size() + triggerMap_["zee"].size());
  triggerMap_["zll"] = triggerMap_["zmm"];
  triggerMap_["zll"].insert(triggerMap_["zll"].end(), triggerMap_["zee"].begin(), triggerMap_["zee"].end());
  triggerMap_["photon"] = {"51"};
  triggerMap_["sig"] = {"29", "33"};

}
