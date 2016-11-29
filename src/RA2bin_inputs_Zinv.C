//
//  Zinv inputs to RA2/b SUSY search fit
//

#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TMath.h"
using TMath::Sqrt; using TMath::Power;

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstring>

namespace Sample {
  enum sampleChoice {Signal, HDP, LDP};
}
using namespace Sample;

// void RA2bin_inputs_Zinv(Int_t doSample = Signal,
void RA2bin_inputs_Zinv(sampleChoice doSample = Signal,
			const TString gJetsFnRoot = TString("gJets"),
			const TString DRfnRoot = TString("DR"),
			const TString DYfnRoot = TString("DY"),
			const TString MCfileName = TString("ZinvMCttzMC174bin.root"),
			const Float_t MClumiScale = 1.0,
			const Bool_t expandDYkin = false) {
  //
  // Supply input flat files with names
  //  <gJetsFnRoot>_signal.dat for doSample = Signal
  //  <gJetsFnRoot>_htp.dat for doSample = HDP
  //  <gJetsFnRoot>_ltp.dat for doSample = LDP
  // If MC root file has lumi different from that of the data, supply 
  //  MClumiScale = dataLumi / MClumi.
  //

  gROOT->Reset();

  // use the 'plain' style for plots (white backgrounds, etc)
  cout << "...using style 'Plain'\n";

  gROOT->SetStyle("Plain");

  // use bold lines and markers
  gStyle->SetMarkerStyle(8);  // non-scalable dot
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(4);
  gStyle->SetFuncStyle(1);
  // gStyle->SetFuncWidth(2);

  //..Get rid of X error bars
  gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //  ===============================================================
  //
  // Prototypes for helper functions
  //
  Int_t getData_gJets(const char* fileName,
		      std::vector <std::vector<Float_t> >& Ngobs,
		      std::vector <std::vector<Float_t> >& NgobsEB,
		      std::vector <std::vector<Float_t> >& NgobsEE,
		      std::vector <std::vector<Float_t> >& ZgR,
		      std::vector <std::vector<Float_t> >& ZgRerr,
		      std::vector <std::vector<Float_t> >& ZgRerrUp,
		      std::vector <std::vector<Float_t> >& ZgRerrLow,
		      std::vector <std::vector<Float_t> >& gFrag,
		      std::vector <std::vector<Float_t> >& gFragErr,
		      std::vector <std::vector<Float_t> >& gPur,
		      std::vector <std::vector<Float_t> >& gPurErr,
		      std::vector <std::vector<Float_t> >& ZgDR,
		      std::vector <std::vector<Float_t> >& ZgDRerrUp,
		      std::vector <std::vector<Float_t> >& ZgDRerrLow);
  Int_t getData_DR(const char* fileName,
		   std::vector <std::vector<Float_t> >& ZgDR,
		   std::vector <std::vector<Float_t> >& ZgDRerrUp,
		   std::vector <std::vector<Float_t> >& ZgDRerrLow,
		   Float_t& DRscaleErr,
		   Float_t& DY0bPurErr,
		   Float_t& DYtrigEffErr,
		   Float_t& LeptonSFerr,
		   Float_t& btag_SFerr);
  Int_t getData_DR(const char* fileName,
		   std::vector <std::vector<Float_t> >& ZgDR,
		   std::vector <std::vector<Float_t> >& ZgDRerrUp,
		   std::vector <std::vector<Float_t> >& ZgDRerrLow);
  Int_t getData_DY(const char* fileName,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0stat,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0MCstat,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysUp,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysLow,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysKin,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysPur);
  void setCorrelationLabels(TH1* histo, const unsigned correlVar, const int minCorrelBin=9999);

  Float_t ignoreBin = -999.;
  TString dataFile_gJets;
  TString dataFile_DR;
  TString dataFile_DY;
  TString output_rootFile;
  TString output_plotFile;
  Float_t histoXlabelSize;
  Float_t canvasBottomMargin;
  Int_t MaxNjets, MaxNb, MaxKin, MaxKinDY;
  if (doSample == LDP || doSample == HDP) {
    MaxNjets = 5;
    MaxNb = 4;
    MaxKin = 13;
    histoXlabelSize = 0.015;
    canvasBottomMargin = 0.15;
    if (doSample == LDP) {
      dataFile_gJets = gJetsFnRoot+TString("_ldp.dat");
      dataFile_DR = DRfnRoot+TString("_ldp.dat");
      dataFile_DY = DYfnRoot+TString( "_ldp.dat");
      output_rootFile = "ZinvHistos_ldp.root";
      output_plotFile = "ZinvBGpred_ldp.png";
    } else if (doSample == HDP) {
      dataFile_gJets = gJetsFnRoot+TString("_hdp.dat");
      dataFile_DR = DRfnRoot+TString("_hdp.dat");
      dataFile_DY = DYfnRoot+TString("_hdp.dat");
      output_rootFile = "ZinvHistos_hdp.root";
      output_plotFile = "ZinvBGpred_hdp.png";
    }
  } else {
    MaxNjets = 5;
    MaxNb = 4;
    MaxKin = 10;
    histoXlabelSize = 0.04;
    canvasBottomMargin = 0.35;
    // dataFile_gJets = "gJets_signal.dat";
    dataFile_gJets = gJetsFnRoot+TString("_signal.dat");
    dataFile_DR = DRfnRoot+TString("_signal.dat");
    dataFile_DY = DYfnRoot+TString("_signal.dat");
    output_rootFile = "ZinvHistos.root";
    output_plotFile = "ZinvBGpred.png";
  }
  MaxKinDY = expandDYkin ? MaxKin : 1;
  Int_t MaxBins;
  if (doSample == Signal) {
    MaxBins = (MaxNjets*MaxNb-1)*MaxKin-2*2*MaxNb;  // Exclude (Njets0,Nb3) && (Njets3,4;HT0,3)
  } else {
    MaxBins = (MaxNjets*MaxNb-1)*MaxKin-2*3*MaxNb;  // Exclude (Njets0,Nb3) && (Njets3,4;HT0,3,6)
  }
  Int_t bin;  // Bin index for loops

  // N(b=0) from gamma+jets
  std::vector <std::vector<Float_t> > Ngobs(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > NgobsEB(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > NgobsEE(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgR(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerr(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrUp(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrLow(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrEffUp(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrEffLow(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gFrag(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gFragErr(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gFragErrEff(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gPur(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gPurErr(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gPurErrEff(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDR(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDRerrUp(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDRerrLow(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDR_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDRerrUp_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgDRerrLow_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  Float_t DRscaleErr=0, DY0bPurErr=0, DYtrigEffErr=0, LeptonSFerr=0, btagSFerr=0;

  // Ratio Nb/N0 from DY
  std::vector < std::vector <std::vector<Float_t> > >
    DYvalues(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYstat(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYMCstat(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYsysNjUp(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYsysNjLow(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYsysKin(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  std::vector < std::vector <std::vector<Float_t> > >
    DYsysPur(MaxNjets, std::vector< std::vector<Float_t> >(MaxNb, std::vector<Float_t>(MaxKinDY, 0)));
  // Read gJets, DY data from files
  if (0 != getData_gJets(dataFile_gJets, Ngobs, NgobsEB, NgobsEE, ZgR, ZgRerr, ZgRerrUp, ZgRerrLow, gFrag, gFragErr,
			 gPur, gPurErr, ZgDR_from_gJets, ZgDRerrUp_from_gJets, ZgDRerrLow_from_gJets)) {
    cout << "Failed to get data." << endl;
    return;
  }
  if (0 != getData_DR(dataFile_DR, ZgDR, ZgDRerrUp, ZgDRerrLow, DRscaleErr, DY0bPurErr, DYtrigEffErr, LeptonSFerr, btagSFerr)) {
    cout << "Failed to get data." << endl;
    return;
  }
  if (0 != getData_DY(dataFile_DY, DYvalues, DYstat, DYMCstat, DYsysNjUp, DYsysNjLow, DYsysKin, DYsysPur)) {
    cout << "Failed to get data." << endl;
    return;
  }

  Float_t gPurAv = 0, gPurErrAv = 0, NobsSum = 0, ZgRAv = 0, ZgRerrUpAv = 0, ZgRerrLowAv = 0;
  cout << "DRscaleErr = " << DRscaleErr << ", DY0bPurErr = " << DY0bPurErr << ", DYtrigEffErr = " << DYtrigEffErr << ", LeptonSFerr = " <<LeptonSFerr  << ", btagSFerr = " << btagSFerr << endl;
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      bin++;
      cout << bin << "  Ngobs = " << Ngobs[ijet][ikin] << "  ZgR = " << ZgR[ijet][ikin] 
           << "  ZgRerr = " << ZgRerr[ijet][ikin]
           << "  ZgRerrUp = " << ZgRerrUp[ijet][ikin]
           << "  ZgRerrLow = " << ZgRerrLow[ijet][ikin]
           << "  gPur = " << gPur[ijet][ikin]
           << "  gPurErr = " << gPurErr[ijet][ikin]
           << "  ZgDR = " << ZgDR[ijet][ikin]
           << "  ZgDRerrUp = " << ZgDRerrUp[ijet][ikin]
           << "  ZgDRerrLow = " << ZgDRerrLow[ijet][ikin] << endl;
      NobsSum += Ngobs[ijet][ikin];
      gPurAv += Ngobs[ijet][ikin]*gPur[ijet][ikin];
      gPurErrAv += Ngobs[ijet][ikin]*gPurErr[ijet][ikin]*gPur[ijet][ikin];  // Absolute error on gPur
      ZgRAv += Ngobs[ijet][ikin]*ZgR[ijet][ikin];
      ZgRerrUpAv += Ngobs[ijet][ikin]*ZgRerrUp[ijet][ikin]*ZgR[ijet][ikin];  // Absolute upper error on ZgR
      ZgRerrLowAv += Ngobs[ijet][ikin]*ZgRerrLow[ijet][ikin]*ZgR[ijet][ikin];  // Absolute lower error on ZgR
    }
  }
  cout << endl;
  gPurAv /= NobsSum;  // Average gPur
  gPurErrAv /= NobsSum;  // Average gPur absolute uncertainty
  ZgRAv /= NobsSum;  // Average ZgR
  ZgRerrUpAv /= NobsSum;  // Average ZgR absolute upper uncertainty
  ZgRerrLowAv /= NobsSum;  // Average ZgR absolute lower uncertainty
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<MaxKinDY; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	bin++;
	cout << bin << "  DYvalue = " << DYvalues[ijet][ib][ikin]
	     << "  DYstat = " << DYstat[ijet][ib][ikin]
	     << "  DYMCstat = " << DYMCstat[ijet][ib][ikin]
             << "  DYsysNjUp = " << DYsysNjUp[ijet][ib][ikin]
	     << "  DYsysNjLow = " << DYsysNjLow[ijet][ib][ikin]
             << "  DYsysKin = " << DYsysKin[ijet][ib][ikin]
	     << "  DYsysPur = " << DYsysPur[ijet][ib][ikin] << endl;
      }
    }
  }
  cout << endl << endl;

  TFile *MCfile = 0;
  if (doSample == Signal) {
    // Input file for MC histogram
    cout << "Initializing input file " << MCfileName << endl;
    MCfile = TFile::Open(MCfileName, "READ");
    if (!MCfile) {
      cout << "Didn't find MC file " << MCfileName << endl;
    }
  }

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:"); // Without these 2 lines
  rootdir->cd();                     // the histo disappears when the macro exits.
  
  TH1F* hTemplate = new TH1F("hTemplate", "template for Zinv histograms", MaxBins, 0.5, MaxBins+0.5);
  // Tune up the template histogram
  // hTemplate->Reset();
  // hTemplate->GetXaxis()->SetLabelSize(histoXlabelSize);
  // hTemplate->GetYaxis()->SetLabelSize(0.05);
  // hTemplate->GetYaxis()->SetTitleSize(0.05);
  hTemplate->SetOption("HIST");
  hTemplate->SetMarkerSize(0);

  // hCorrelTemplate is the starting point for x-axis labels that will encode the
  // correlation patterns for the integration fitter.
  TH1F* hCorrelTemplate = (TH1F*)hTemplate->Clone("hCorrelTemplate");
  if (doSample == Signal) {  // Leave x-axis labels as numbers for hdp, ldp samples
    // Set the x-axis labels for default (no) correlations
    bin = 0;
    for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
      for (Int_t ib=0; ib<MaxNb; ++ib) {
	if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
	for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
	  if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	  if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	  bin++;
	  char label[25];
	  Int_t iMHT = ikin/3;
	  if(ikin>=8) {
	    iMHT=3;
	  }
	  Int_t iHT = ikin;
	  sprintf(label, "%5s%1d%6s%1d%4s%1d%3s%1d", "_NJets", ijet, "_BTags", ib, "_MHT", iMHT, "_HT", iHT);
	  // cout << label << endl;
	  hCorrelTemplate->GetXaxis()->SetBinLabel(bin, label);
	}
      }
    }
  }

  // Output file
  TFile *ZinvHistos = TFile::Open(output_rootFile, "RECREATE");

  gStyle->SetOptStat("n");

  // Book histograms for N(b=0) yield factors from gamma+jets
  //  NZvv_pred = Ngobs * ZgR * gPur * ZgDR
  TH1F *hzvvgJNobs = (TH1F*)hCorrelTemplate->Clone("hzvvgJNobs");
  if (doSample == Signal) setCorrelationLabels(hzvvgJNobs, 4);  
  hzvvgJNobs->GetYaxis()->SetTitle("gamma+jets yield");

  TH1F *hgJstat = (TH1F*)hTemplate->Clone("hgJstat");
  hgJstat->GetYaxis()->SetTitle("gamma+jets yield stat error");

  TH1F *hzvvTF = (TH1F*)hCorrelTemplate->Clone("hzvvTF");
  hzvvTF->GetYaxis()->SetTitle("Zinv transfer factor");

  TH1F *hgJZgR = (TH1F*)hTemplate->Clone("hgJZgR");
  hgJZgR->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio");

  // TH1F *hzvvgJZgRerr = (TH1F*)hCorrelTemplate->Clone("hzvvgJZgRerr");
  // if (doSample == Signal) setCorrelationLabels(hzvvgJZgRerr, 4);  
  // hzvvgJZgRerr->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  TH1F *hgJZgRerr = (TH1F*)hTemplate->Clone("hgJZgRerr");
  hgJZgRerr->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  // The uncertainty gJZgRerrUp, Low cancels in between RZg and DR
  TH1F *hgJZgRerrUp = (TH1F*)hCorrelTemplate->Clone("hgJZgRerrUp");
  if (doSample == Signal) setCorrelationLabels(hgJZgRerrUp, 15);  
  hgJZgRerrUp->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  TH1F *hgJZgRerrLow = (TH1F*)hCorrelTemplate->Clone("hgJZgRerrLow");
  if (doSample == Signal) setCorrelationLabels(hgJZgRerrLow, 15);  
  hgJZgRerrLow->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  TH1F *hgJPur = (TH1F*)hTemplate->Clone("hgJPur");
  hgJPur->GetYaxis()->SetTitle("gamma+jets purity");

  TH1F *hzvvgJPurErr = (TH1F*)hCorrelTemplate->Clone("hzvvgJPurErr");
  if (doSample == Signal) setCorrelationLabels(hzvvgJPurErr, 15);  
  hzvvgJPurErr->GetYaxis()->SetTitle("gamma+jets purity error");

  TH1F *hZgDR = (TH1F*)hTemplate->Clone("hZgDR");
  hZgDR->GetYaxis()->SetTitle("gamma+jets double ratio");

  TH1F *hzvvZgDRerrUp = (TH1F*)hCorrelTemplate->Clone("hzvvZgDRerrUp");
  if (doSample == Signal) setCorrelationLabels(hzvvZgDRerrUp, 4);  
  hzvvZgDRerrUp->GetYaxis()->SetTitle("gamma+jets double ratio error");

  TH1F *hzvvZgDRerrLow = (TH1F*)hCorrelTemplate->Clone("hzvvZgDRerrLow");
  if (doSample == Signal) setCorrelationLabels(hzvvZgDRerrLow, 4);  
  hzvvZgDRerrLow->GetYaxis()->SetTitle("gamma+jets double ratio error");

  TH1F *hzvvScaleErr = (TH1F*)hCorrelTemplate->Clone("hzvvScaleErr");
  if (doSample == Signal) setCorrelationLabels(hzvvScaleErr, 15);  
  hzvvScaleErr->GetYaxis()->SetTitle("Zinv global scale error");

  // Book histograms for Nb/N0 ratio values and error components
  //  Nb/N0 pred = DYvalues
  TH1F *hDYvalue = (TH1F*)hTemplate->Clone("hDYvalue");
  hDYvalue->GetYaxis()->SetTitle("DY ratio to 0b value");

  TH1F *hzvvDYstat = (TH1F*)hCorrelTemplate->Clone("hzvvDYstat");
  if (doSample == Signal) setCorrelationLabels(hzvvDYstat, 11, 2999);  
  hzvvDYstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");

  TH1F *hzvvDYMCstat = (TH1F*)hCorrelTemplate->Clone("hzvvDYMCstat");
  if (doSample == Signal) setCorrelationLabels(hzvvDYMCstat, 3);  
  hzvvDYMCstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");

  TH1F *hzvvDYsysNjUp = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysNjUp");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysNjUp, 3);  
  hzvvDYsysNjUp->GetYaxis()->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");

  TH1F *hzvvDYsysNjLow = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysNjLow");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysNjLow, 3);  
  hzvvDYsysNjLow->GetYaxis()->SetTitle("DY ratio to 0b syst- error Nj extrapolation");

  TH1F *hzvvDYsysKin = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysKin");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysKin, 0);  
  hzvvDYsysKin->GetYaxis()->SetTitle("DY ratio to 0b syst error kinematics dependence");

  TH1F *hzvvDYsysPur = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysPur");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysPur, 15, 2299);  
  hzvvDYsysPur->GetYaxis()->SetTitle("DY ratio to 0b syst error purity");

  // For combined results
  Float_t statErr[MaxBins], sysUp[MaxBins], sysLow[MaxBins];
  TH1F *ZinvBGpred = (TH1F*)hTemplate->Clone("ZinvBGpred");
  ZinvBGpred->GetYaxis()->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background");
  ZinvBGpred->SetLineColor(kBlue);
  ZinvBGpred->SetMarkerColor(kBlue);
  ZinvBGpred->SetMarkerStyle(8);
  ZinvBGpred->SetMarkerSize(0.5);
  // ZinvBGpred->SetMinimum(-1);
  TH1F *ZinvBGsysUp = (TH1F*)hTemplate->Clone("ZinvBGsysUp");
  ZinvBGsysUp->GetYaxis()->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background");
  TH1F *ZinvBGsysLow = (TH1F*)hTemplate->Clone("ZinvBGsysLow");
  ZinvBGsysLow->GetYaxis()->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background");

  // For combined results, for bins with zero control-sample events
  TH1F *ZinvBG0EVpred = (TH1F*)hTemplate->Clone("ZinvBG0EVpred");
  ZinvBG0EVpred->GetYaxis()->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background, 0 events");
  ZinvBG0EVpred->SetLineColor(kBlue);
  ZinvBG0EVpred->SetMarkerColor(kBlue);
  ZinvBG0EVpred->SetMarkerStyle(8);
  ZinvBG0EVpred->SetMarkerSize(0.5);
  TH1F *ZinvBG0EVsysUp = (TH1F*)hTemplate->Clone("ZinvBG0EVsysUp");
  ZinvBG0EVsysUp->GetYaxis()->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");
  TH1F *ZinvBG0EVsysLow = (TH1F*)hTemplate->Clone("ZinvBG0EVsysLow");
  ZinvBG0EVsysLow->GetYaxis()->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");

  ofstream tableFile;  tableFile.open("table_for_AN.txt");
  char buf[256];

  // Fill histograms
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	bin++;
	Int_t ikinDY = MaxKinDY == MaxKin ? ikin : 0;

        //
        // Take into account that photon purity and ZgR syst errors partly cancel in the prediction
        //
	gPurErrEff[ijet][ikin] = fabs( gPurErr[ijet][ikin] - gPurErrAv / gPurAv );  // Frac. error on gPur/<gPur>
	ZgRerrEffUp[ijet][ikin] = fabs( ZgRerrUp[ijet][ikin] - ZgRerrUpAv / ZgRAv );  // Frac. upper error on ZgR/<ZgR>
	ZgRerrEffLow[ijet][ikin] = fabs( ZgRerrLow[ijet][ikin] - ZgRerrLowAv / ZgRAv );  // Frac. upper error on ZgR/<ZgR>
        // cout << "gPur err, nominal: " << gPurErr[ijet][ikin] << ", eff: " << gPurErrEff
        //      << ",    ZgR errUp nominal: " << ZgRerrUp[ijet][ikin] << ", eff: " << ZgRerrEffUp[ijet][ikin]
        //      << ",    ZgR errLow nominal: " << ZgRerrLow[ijet][ikin] << ", eff: " << ZgRerrEffLow[ijet][ikin] << endl;

	hzvvgJNobs->SetBinContent(bin, Ngobs[ijet][ikin]);
	Ngobs[ijet][ikin] > 0 ? hgJstat->SetBinContent(bin, 1+1/Sqrt(Ngobs[ijet][ikin])) : hgJstat->SetBinContent(bin, 1+0);
	hgJZgR->SetBinContent(bin, ZgR[ijet][ikin]);
	hZgDR->SetBinContent(bin, ZgDR[ijet][ikin]);
	hgJZgRerr->SetBinContent(bin, 1+ZgRerr[ijet][ikin]);  // uncorrelated, combined into ZgDRErrUp,Low
	hgJZgRerrUp->SetBinContent(bin, 1+ZgRerrEffUp[ijet][ikin]);  // correlated
	hgJZgRerrLow->SetBinContent(bin, 1-ZgRerrEffLow[ijet][ikin]);  // (but Eff error = 0 if error is bin-independent)
	hgJPur->SetBinContent(bin, gPur[ijet][ikin]);
	hzvvgJPurErr->SetBinContent(bin, 1+gPurErrEff[ijet][ikin]);  // correlated
	// Avoid degenerate factors by combining ZgR stat and ZgDR errors into ZgDRErrUp,Low
	// hzvvZgDRerrUp->SetBinContent(bin, 1+ZgDRerrUp[ijet][ikin]);
	// hzvvZgDRerrLow->SetBinContent(bin, 1-ZgDRerrLow[ijet][ikin]);
	hzvvZgDRerrUp->SetBinContent(bin, 1+Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(btagSFerr, 2) + Power(ZgDRerrUp[ijet][ikin], 2)));  // uncorrelated
	hzvvZgDRerrLow->SetBinContent(bin, 1-Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(btagSFerr, 2) + Power(ZgDRerrLow[ijet][ikin], 2)));
	hzvvScaleErr->SetBinContent(bin, 1+Sqrt(Power(DRscaleErr, 2) + Power(DY0bPurErr, 2) + Power(DYtrigEffErr, 2) + Power(LeptonSFerr, 2)));

  	hDYvalue->SetBinContent(bin, DYvalues[ijet][ib][ikinDY]);
  	hzvvDYstat->SetBinContent(bin, 1+DYstat[ijet][ib][ikinDY]);
  	hzvvDYMCstat->SetBinContent(bin, 1+DYMCstat[ijet][ib][ikinDY]);
  	hzvvDYsysNjUp->SetBinContent(bin, 1+DYsysNjUp[ijet][ib][ikinDY]);
  	hzvvDYsysNjLow->SetBinContent(bin, 1-DYsysNjLow[ijet][ib][ikinDY]);
  	hzvvDYsysKin->SetBinContent(bin, 1+DYsysKin[ijet][ib][ikinDY]);
  	hzvvDYsysPur->SetBinContent(bin, 1+DYsysPur[ijet][ib][ikinDY]);

  	// cout << bin
	//      << " " << hDYvalue->GetXaxis()->GetBinLabel(bin)
	//      << "  " << hDYvalue->GetBinContent(bin) << endl;

	Float_t transferFactor =  ZgDR[ijet][ikin] * ZgR[ijet][ikin] * gPur[ijet][ikin] * DYvalues[ijet][ib][ikinDY];
	// cout << Ngobs[ijet][ikin] << " " << ZgR[ijet][ikin] << " " << gPur[ijet][ikin] << " " << DYvalues[ijet][ib][ikin] << " " << ZgDR[ijet][ikin] << endl;
	Float_t thisNgobs = Ngobs[ijet][ikin] > 0 ? Ngobs[ijet][ikin] : 1.0;
	Float_t ZinvValue = thisNgobs * transferFactor;
	hzvvTF->SetBinContent(bin, transferFactor);
	Float_t wtStat = 1/thisNgobs;
	// statErr[bin-1] = ZinvValue*Sqrt(wtStat + Power(DYstat[ijet][ib][ikinDY], 2));
	statErr[bin-1] = ZinvValue*Sqrt(wtStat);
	sysUp[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				      + Power(ZgRerrEffUp[ijet][ikin], 2)
				      + Power(gPurErrEff[ijet][ikin], 2)
				      + Power(ZgDRerrUp[ijet][ikin], 2)
				      + Power(btagSFerr, 2)
				      + Power(DRscaleErr, 2)
				      + Power(DY0bPurErr, 2)
				      + Power(DYtrigEffErr, 2)
				      + Power(LeptonSFerr, 2)
				      + Power(DYstat[ijet][ib][ikinDY], 2)
				      + Power(DYMCstat[ijet][ib][ikinDY], 2)
				      + Power(DYsysNjUp[ijet][ib][ikinDY], 2)
				      + Power(DYsysKin[ijet][ib][ikinDY], 2)
				      + Power(DYsysPur[ijet][ib][ikinDY], 2));
	sysLow[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				       + Power(ZgRerrEffLow[ijet][ikin], 2)
				       + Power(gPurErrEff[ijet][ikin], 2)
				       + Power(ZgDRerrLow[ijet][ikin], 2)
				       + Power(btagSFerr, 2)
				       + Power(DRscaleErr, 2)
				       + Power(DY0bPurErr, 2)
				       + Power(DYtrigEffErr, 2)
				       + Power(LeptonSFerr, 2)
				       + Power(DYstat[ijet][ib][ikinDY], 2)
				       + Power(DYMCstat[ijet][ib][ikinDY], 2)
				       + Power(DYsysNjLow[ijet][ib][ikinDY], 2)
				       + Power(DYsysKin[ijet][ib][ikinDY], 2)
				       + Power(DYsysPur[ijet][ib][ikinDY], 2));
	cout << bin << "  " << ZinvValue << " +/- " << statErr[bin-1] << " + " << sysUp[bin-1]
	     << " - " << sysLow[bin-1] << "  TF*Ngobs = "
	     << hzvvTF->GetBinContent(bin) * hzvvgJNobs->GetBinContent(bin) << endl;
	if (Ngobs[ijet][ikin] > 0) {
	  ZinvBGpred->SetBinContent(bin, ZinvValue);
	  ZinvBGsysUp->SetBinContent(bin, sysUp[bin-1]);
	  ZinvBGsysLow->SetBinContent(bin, sysLow[bin-1]);
	} else {
	  ZinvBG0EVpred->SetBinContent(bin, ZinvValue);
	  ZinvBG0EVsysUp->SetBinContent(bin, sysUp[bin-1]);
	  ZinvBG0EVsysLow->SetBinContent(bin, sysLow[bin-1]);
	}

	if (ib == 0) {
          //  Write out LaTex for AN table
	  if (ikin == 0 || (ijet > 2 && ikin == 1)) tableFile << endl;
	  sprintf(buf, "%5.0f & %5.0f & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%6.1f\\pm%4.1f^{+%4.1f}_{-%4.1f}$ \\\\\n",
		  NgobsEB[ijet][ikin], NgobsEE[ijet][ikin],
		  ZgR[ijet][ikin], ZgR[ijet][ikin]*ZgRerr[ijet][ikin], ZgR[ijet][ikin]*ZgRerrUp[ijet][ikin],
		  ZgR[ijet][ikin]*ZgRerrLow[ijet][ikin], 
		  ZgDR[ijet][ikin], DRscaleErr, ZgDR[ijet][ikin]*ZgDRerrUp[ijet][ikin],
		  ZgDR[ijet][ikin]*ZgDRerrLow[ijet][ikin],
		  ZinvValue, statErr[bin-1], sysUp[bin-1], sysLow[bin-1]);
	  tableFile << buf;
	}
      }  // ikin
    }  // ib
  }  // ijet
  tableFile.close();

  TGraphAsymmErrors* ZinvBGsyst = new TGraphAsymmErrors(ZinvBGpred);
  ZinvBGsyst->SetLineColorAlpha(kRed, 0);
  ZinvBGsyst->SetLineWidth(-1);
  ZinvBGsyst->SetFillColor(kRed-10);
  ZinvBGsyst->SetFillStyle(1001);
  TGraphAsymmErrors* ZinvBG0EVsyst = new TGraphAsymmErrors(ZinvBG0EVpred);
  ZinvBG0EVsyst->SetLineColor(-1);
  // ZinvBG0EVsyst->SetLineStyle(1);
  // ZinvBG0EVsyst->SetLineWidth(2);
  // gStyle->SetEndErrorSize(3);
  for (Int_t bin0=0; bin0<MaxBins; ++bin0) {
    if (ZinvBGpred->GetBinContent(bin0+1) > 0) {
      ZinvBGpred->SetBinError(bin0+1, statErr[bin0]);
      ZinvBGsyst->SetPointError(bin0, 0.5, 0.5, sysLow[bin0], sysUp[bin0]);
    } else {
      ZinvBG0EVpred->SetBinError(bin0+1, statErr[bin0]);
      ZinvBG0EVsyst->SetPointError(bin0, 0.5, 0.5, sysLow[bin0], sysUp[bin0]);
    }
  }

  TH1F *hZinvMCbin = 0;
  if (doSample == Signal) {
    // Zinv MC 
    TH1F* hZinvMC = 0;
    TString MChistoname("plot_zinv_nj5_nb4_kin10_1");
    // TString MChistoname("plot_zinv_nj4_nb4_kin10_1");
    if (MCfile != 0) {
      if (!MCfile->Get(MChistoname)) {
	cout << "No MC histogram found in this file and directory" << endl;
      } else {
	hZinvMC = (TH1F*) MCfile->Get(MChistoname);
	hZinvMCbin = (TH1F*)hTemplate->Clone("hZinvMCbin");
	hZinvMCbin->GetYaxis()->SetTitle("MC Z#rightarrow#nu#bar{#nu}");
	for (Int_t bin=1; bin<=MaxBins; ++bin) {
	  hZinvMCbin->SetBinContent(bin, hZinvMC->GetBinContent(bin));
	  hZinvMCbin->SetBinError(bin, hZinvMC->GetBinError(bin));
	}
	hZinvMCbin->SetLineColor(kGreen+2);
	hZinvMCbin->Scale(MClumiScale);
      }
    }
  }

  TCanvas *Cvalue = new TCanvas("Cvalue", "DY ratio to zero b", 1000, 500);
  // Cvalue->SetTopMargin(.03);
  // Cvalue->SetBottomMargin(canvasBottomMargin);
  // Cvalue->SetRightMargin(.03);
  gStyle->SetOptStat("");

  ZinvBGpred->Draw("e0");
  ZinvBGsyst->Draw("2 same");
  if (hZinvMCbin != 0) hZinvMCbin->Draw("hist e0 same");
  ZinvBGpred->Draw("e0 same");
  Cvalue->SetLogy();
  Cvalue->Update();

  TLegend* Zinvlegend;
  if (doSample == Signal)
    Zinvlegend = new TLegend(.72, .78, .92, .9, "");
  else
    Zinvlegend = new TLegend(.72, .80, .92, .9, "");
  Zinvlegend->AddEntry(ZinvBGpred, "Yields, stat. errors", "pe");
  Zinvlegend->AddEntry(ZinvBGsyst, "Syst. errors", "f");
  if (hZinvMCbin != 0) 
    Zinvlegend->AddEntry(hZinvMCbin, "Z#rightarrow#nu#bar{#nu} simulation", "le");
  Zinvlegend->Draw();
  Cvalue->Update();
  Cvalue->SaveAs(output_plotFile);

  ZinvHistos->Write();


  TCanvas *Cerrors = new TCanvas("Cerrors", "Error components", 1000, 500);
  // Cerrors->SetTopMargin(.03);
  // Cerrors->SetBottomMargin(canvasBottomMargin);
  // Cerrors->SetRightMargin(.03);
  hgJstat->SetMinimum(-0.2);  hgJstat->SetMaximum(2.2);
  hgJstat->SetLineColor(1);  hgJstat->Draw();
  hgJZgRerr->SetLineColor(2);  hgJZgRerr->Draw("same");
  // hgJZgRerrUp->SetLineColor(36);  hgJZgRerrUp->Draw("same");
  // hgJZgRerrLow->SetLineColor(4);  hgJZgRerrLow->Draw("same");
  hzvvgJPurErr->SetLineColor(7);  hzvvgJPurErr->Draw("same");
  hzvvScaleErr->SetLineColor(6);  hzvvScaleErr->Draw("same");
  hzvvZgDRerrUp->SetLineColor(5);  hzvvZgDRerrUp->Draw("same");
  hzvvZgDRerrLow->SetLineColor(5);  hzvvZgDRerrLow->Draw("same");
  hzvvDYstat->SetLineColor(8);  hzvvDYstat->Draw("same");
  hzvvDYMCstat->SetLineColor(9);  hzvvDYMCstat->Draw("same");
  hzvvDYsysNjUp->SetLineColor(3);  hzvvDYsysNjUp->Draw("same");
  hzvvDYsysNjLow->SetLineColor(3);  hzvvDYsysNjLow->Draw("same");
  hzvvDYsysKin->SetLineColor(42);  hzvvDYsysKin->Draw("same");
  hzvvDYsysPur->SetLineColor(46);  hzvvDYsysPur->Draw("same");
  TLegend* ErrorsLegend;
  if (doSample == Signal)
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  else
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  ErrorsLegend->AddEntry(hgJstat, "gJstat");
  ErrorsLegend->AddEntry(hgJZgRerr, "RZgErr");
  // ErrorsLegend->AddEntry(hgJZgRerrUp, "RZgUp");
  // ErrorsLegend->AddEntry(hgJZgRerrLow, "RZgLow");
  ErrorsLegend->AddEntry(hzvvgJPurErr, "gJpur");
  ErrorsLegend->AddEntry(hzvvScaleErr, "Scale");
  ErrorsLegend->AddEntry(hzvvZgDRerrUp, "DR");
  // ErrorsLegend->AddEntry(hzvvZgDRerrLow, "DRlow");
  ErrorsLegend->AddEntry(hzvvDYstat, "DYstat");
  ErrorsLegend->AddEntry(hzvvDYMCstat, "DYMCstat");
  ErrorsLegend->AddEntry(hzvvDYsysNjUp, "DYNj");
  // ErrorsLegend->AddEntry(hzvvDYsysNjLow, "DYNjLow");
  ErrorsLegend->AddEntry(hzvvDYsysKin, "DYkin");
  ErrorsLegend->AddEntry(hzvvDYsysPur, "DYpur");
  ErrorsLegend->Draw();
  Cerrors->Update();
  if (doSample == Signal)  Cerrors->SaveAs("errorComponents.png");


}  // ---------------------------------------------------------------

Int_t getData_gJets(const char* fileName,
		    std::vector <std::vector<Float_t> >& Ngobs,
		    std::vector <std::vector<Float_t> >& NgobsEB,
		    std::vector <std::vector<Float_t> >& NgobsEE,
		    std::vector <std::vector<Float_t> >& ZgR,
		    std::vector <std::vector<Float_t> >& ZgRerr,
		    std::vector <std::vector<Float_t> >& ZgRerrUp,
		    std::vector <std::vector<Float_t> >& ZgRerrLow,
		    std::vector <std::vector<Float_t> >& gFrag,
		    std::vector <std::vector<Float_t> >& gFragErr,
		    std::vector <std::vector<Float_t> >& gPur,
		    std::vector <std::vector<Float_t> >& gPurErr,
		    std::vector <std::vector<Float_t> >& ZgDR,
		    std::vector <std::vector<Float_t> >& ZgDRerrUp,
		    std::vector <std::vector<Float_t> >& ZgDRerrLow) {
  //
  // Read data from the file, parse, and set data arrays.
  //
  /*
  Notes on correlations:
  gPurErr is correlated among the 40 bins for which Nb = 0.  1111 = 15
  ZgRerrUp,Low are correlated among the 40 bins.             1111 = 15
  ZgRerr is uncorrelated among the 40 bins.                  0100 = 4, combined with
  ZgDRerrUp, Low are uncorrelated among the 40 bins.         0100 = 4  this one
  As of 10 Jul 2016 we combine these last two in ZgDRerrUp, Low.
   */
  Int_t Nrow = Ngobs.size();
  Int_t Ncol = Ngobs[0].size();
  ifstream dataStream;
  cout << fileName << endl;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
  cout << "Open failed for file " << fileName << endl;
    return 1; // exit if file not found
  }
  char buf[512];  // line buffer
  do {
    dataStream.getline(buf, 512);  // Discard (pre-)header lines
    // cout << buf << endl;
    if (dataStream.eof()) {
      cout << "EOF found prematurely while searching for gJets header line" << endl;
      return 1;
    }
  }
  while (strncmp(buf, " binIndex", 9) != 0);
  // while (strncmp(buf, "no ", 3) != 0);
  // while (strncmp(buf, "no  Nobs", 8) != 0);
  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ikin=0; ikin<Ncol; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (Ncol > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      // cout << ijet << " " << ikin << endl;
      // Read next line of the stream and extract its data
      dataStream.getline(buf, 512);
      // cout << buf << endl;
      if (dataStream.eof()) {
	cout << "EOF found prematurely while searching for jJets data" << endl;
	return 1;
      }
      const char* token[20] = {}; // initialize to 0
      Int_t n = 0;
      // parse the line
      token[n] = strtok(buf, ":"); // first token
      if (token[n] == 0) {
	cout << n << endl;
	cout << "gj Empty line encountered" << endl;
	return 1;
      }
      n++; token[n] = strtok(0, "|");  // for file with MC columns
      n++; token[n] = strtok(0, "|");  // for file with MC EB, EE columns
      n++; token[n] = strtok(0, "|");  // for file with MC EB, EE columns
      n++; token[n] = strtok(0, "|");             // Nobs
      sscanf(token[n], "%f", &Ngobs[ijet][ikin]);
      n++; token[n] = strtok(0, "|");            // NgobsEB
      sscanf(token[n], "%f", &NgobsEB[ijet][ikin]);
      n++; token[n] = strtok(0, "|");
      n++; token[n] = strtok(0, "|");            // NgobsEE
      sscanf(token[n], "%f", &NgobsEE[ijet][ikin]);
      n++; token[n] = strtok(0, "|");
      n++; token[n] = strtok(0, "(");             // RZg
      sscanf(token[n], "%f", &ZgR[ijet][ikin]);
      n++; token[n] = strtok(0, ",");             // RZgErr  0100 = 8
      sscanf(token[n], "%f", &ZgRerr[ijet][ikin]);
      n++; token[n] = strtok(0, "-");             // RZgErrUp  1111 = 15
      sscanf(token[n], "%f", &ZgRerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // RZgErrLow  1111 = 15
      sscanf(token[n], "%f", &ZgRerrLow[ijet][ikin]);
      n++; token[n] = strtok(0, "|");
      n++; token[n] = strtok(0, "(");             // Fragmentation factor
      sscanf(token[n], "%f", &gFrag[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // FragmentationErr  1111 = 15?
      sscanf(token[n], "%f", &gFragErr[ijet][ikin]);
      n++; token[n] = strtok(0, "|");
      n++; token[n] = strtok(0, "(");             // Purity
      sscanf(token[n], "%f", &gPur[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // PurityErr  1111 = 15
      sscanf(token[n], "%f", &gPurErr[ijet][ikin]);
      n++; token[n] = strtok(0, " ");
      n++; token[n] = strtok(0, "(");             // ZgDR
      sscanf(token[n], "%f", &ZgDR[ijet][ikin]);
      n++; token[n] = strtok(0, "-");             // ZgDRErrUp  0000 = 0 (not used)
      sscanf(token[n], "%f", &ZgDRerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // ZgDRErrLow  0000 = 0 (not used)
      sscanf(token[n], "%f", &ZgDRerrLow[ijet][ikin]);
    }
  }
  return 0;
}  // ---------------------------------------------------------------


Int_t getData_DR(const char* fileName,
		 std::vector <std::vector<Float_t> >& ZgDR,
		 std::vector <std::vector<Float_t> >& ZgDRerrUp,
		 std::vector <std::vector<Float_t> >& ZgDRerrLow,
		 Float_t& DRscaleErr,
		 Float_t& DY0bPurErr,
		 Float_t& DYtrigEffErr,
		 Float_t& LeptonSFerr,
		 Float_t& btagSFerr) {
  //
  // Read data from the file, parse, and set data arrays.
  //
  /*
  Notes on correlations:
  gPurErr is correlated among the 40 bins for which Nb = 0.  1111 = 15
  ZgRerrUp,Low are correlated among the 40 bins.             1111 = 15
  ZgRerr is uncorrelated among the 40 bins.                  0100 = 4, combined with
  ZgDRerrUp, Low are uncorrelated among the 40 bins.         0100 = 4  this one
  As of 10 Jul 2016 we combine these last two in ZgDRerrUp, Low.
   */
  Int_t Nrow = ZgDR.size();
  Int_t Ncol = ZgDR[0].size();
  ifstream dataStream;
  cout << fileName << endl;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
  cout << "Open failed for file " << fileName << endl;
    return 1; // exit if file not found
  }
  char buf[512];  // line buffer
  do {
    dataStream.getline(buf, 512);  // Discard header line
    // cout << buf << endl;
    if (dataStream.eof()) {
      cout << "EOF found prematurely while searching for gJets header line" << endl;
      return 1;
    }
  }
  while (strncmp(buf, "DR", 2) != 0);
  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ikin=0; ikin<Ncol; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (Ncol > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      // cout << ijet << " " << ikin << endl;
      // Read next line of the stream and extract its data
       dataStream.getline(buf, 512);
      // cout << buf << endl;
      if (dataStream.eof()) {
	cout << "EOF found prematurely while searching for DR data" << endl;
	return 1;
      }
      const char* token[40] = {}; // initialize to 0
      Int_t n = 0;
      // parse the line
      token[n] = strtok(buf, "|");                // ZgDR
      if (token[n] == 0) {
	cout << n << endl;
	cout << "DR Empty line encountered" << endl;
	return 1;
      }
      sscanf(token[n], "%f", &ZgDR[ijet][ikin]);
      n++; token[n] = strtok(0, "|");  // Double ratio scale error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DRscaleErr);
      n++; token[n] = strtok(0, "|");             // ZgDRErrUp
      sscanf(token[n], "%f", &ZgDRerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // ZgDRErrLow
      sscanf(token[n], "%f", &ZgDRerrLow[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // 0b purity
      n++; token[n] = strtok(0, "|");             // 0b purity error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DY0bPurErr);
      n++; token[n] = strtok(0, "|");             // 0b trigger eff.
      n++; token[n] = strtok(0, "|");             // 0b trigger eff. error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DYtrigEffErr);
      n++; token[n] = strtok(0, "|");             // lepton SF
      n++; token[n] = strtok(0, "|");             // lepton SF error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &LeptonSFerr);
      n++; token[n] = strtok(0, "|");             // btag SF
      n++; token[n] = strtok(0, "|");             // btag SF error  (not used)
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &btagSFerr);
    }
  }
  return 0;
}  // ---------------------------------------------------------------

Int_t getData_DY(const char* fileName,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0stat,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0MCstat,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0sysUp,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0sysLow,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0sysKin,
		 std::vector <std::vector <std::vector<Float_t> > >& Rb0sysPur) {
  //
  // Read data from the file, parse, and set data arrays.
  //
/*
  Notes on correlations:
  All DY nuisances except DYsysKin are correlated across all kinematic bins.
  DYstat is correlated across Njets bins for Njets >= 7                1011 = 11, 2999
  Rb0MCstat is correlated only over kinematic bins                     0011 = 3
  DYsysNjUp, Low are uncorrelated in Njets, Nbjets                     0111 = 3
  DYsysPur is correlated across Njets bins for Njets >= 5,
    and across Nb bins for Nb >=2                                      1111 = 15, 2299        
  DYsysKin is uncorrelated (though all are zero for Nb = 0).           0000 = 0
*/
  Int_t Nrow = Rb0.size();
  Int_t Ncol = Rb0[0].size();
  Int_t Nplane = Rb0[0][0].size();
  ifstream dataStream;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
  cout << "Open failed for file " << fileName << endl;
    return 1; // exit if file not found
  }
  char buf[512];  // line buffer
  dataStream.getline(buf, 512);  // Discard header line

  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ib=0; ib<Ncol; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<Nplane; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (Nplane > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
        //
        // Read next line of the stream and extract its data
        //
	// cout << ijet << " " << ib << endl;
	dataStream.getline(buf, 512);
	if (dataStream.eof()) {
	  cout << "EOF found prematurely while searching for DY data" << endl;
	  return 1;
	}
	const char* token[20] = {}; // initialize to 0
	Int_t n = 0;
        // parse the line
	token[n] = strtok(buf, "|"); // first token
	if (token[n] == 0) {
	  cout << "Empty line encountered" << endl;
	  return 1;
	}
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0stat[ijet][ib][ikin]);  //   1011 = 11, 2999
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0MCstat[ijet][ib][ikin]);  // 0011 = 3
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysUp[ijet][ib][ikin]);  //  0111 = 3
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysLow[ijet][ib][ikin]);  // 0111 = 3
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysKin[ijet][ib][ikin]);  // 0000 = 0
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysPur[ijet][ib][ikin]);  // 1111 = 15, 2299
      }
    }
  }
  return 0;
}  // ---------------------------------------------------------------

void setCorrelationLabels(TH1* histo, const unsigned correlVar, const int minCorrelBin=9999) {
  TString hname = histo->GetName();
  cout << hname;
  hname.Remove(0,1);  // Strip off the "h" at start of histogram name
  for (Int_t bin=1; bin<=histo->GetNbinsX(); ++bin) {
    TString binLabel = histo->GetXaxis()->GetBinLabel(bin);
    // binLabel = "_NJets<n>_BTags<m>_MHT<p>_HT<q>"
    // Remove substring for variable with correlated uncertainties.
    // correlVar bit 0 HT; bit 1 MHT; bit 2 Nb; bit 3 Njets

    if (correlVar & 8) {
      Ssiz_t index = binLabel.Index("_NJets");
      if (minCorrelBin/1000 % 10 == 9) {
	binLabel.Remove(index, 7);  // Correlated across all NJets bins
      } else {
	TString theBin = binLabel(index+6,1);  // Correlated across some NJets bins
	char newBin[10];  sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), minCorrelBin/1000 % 10));
	binLabel.Replace(index+6, 1, newBin);
      }
    }

    if (correlVar & 4) {
      Ssiz_t index = binLabel.Index("_BTags");
      if (minCorrelBin/100 % 10 == 9) {
    	binLabel.Remove(index, 7);  // Correlated across all BTag bins
      } else {
	TString theBin = binLabel(index+6,1);  // Correlated across some BTag bins
	char newBin[10];  sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), minCorrelBin/100 % 10));
	binLabel.Replace(index+6, 1, newBin);
      }
    }

    if (correlVar & 2) {
      Ssiz_t index = binLabel.Index("_MHT");
      // cout << "minCorrelBin/10 % 10 = " << minCorrelBin/10 % 10 << endl;
      if (minCorrelBin/10 % 10 == 9) {
	// cout << "  removing _MHTx" << endl;
    	binLabel.Remove(index, 5);  // Correlated across all MHT bins
      } else {
	TString theBin = binLabel(index+4,1);  // Correlated across some MHT bins
	char newBin[10];  sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), minCorrelBin/10 % 10));
	binLabel.Replace(index+4, 1, newBin);
      }
    }

    if (correlVar & 1) {
      Ssiz_t index = binLabel.Index("_HT");
      if (minCorrelBin % 10 == 9) {
    	binLabel.Remove(index, 4);  // Correlated across all HT bins
      } else {
	TString theBin = binLabel(index+3,1);  // Correlated across some HT bins
	char newBin[10];  sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), minCorrelBin % 10));
	binLabel.Replace(index+3, 1, newBin);
      }
    }

    binLabel.Prepend(hname);
    if (bin == 135) cout << ",  binLabel(135) = " << binLabel;
    histo->GetXaxis()->SetBinLabel(bin, binLabel.Data());
  }
  cout << endl;
}  // ---------------------------------------------------------------
