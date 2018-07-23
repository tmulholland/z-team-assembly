//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#ifndef RA2BZINVANALYSIS_H
#define RA2BZINVANALYSIS_H

/* #define ISMC */

#include <TString.h>
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
    TString name;
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

#endif
