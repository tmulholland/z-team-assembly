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
			const TString MCfileName = TString("../plots/histograms/ZinvMCttzMC160bin.root"),
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
		      std::vector <std::vector<Float_t> >& gPur,
		      std::vector <std::vector<Float_t> >& gPurErr,
		      std::vector <std::vector<Float_t> >& ZgRdataMC,
		      std::vector <std::vector<Float_t> >& ZgRdataMCerrUp,
		      std::vector <std::vector<Float_t> >& ZgRdataMCerrLow);
  Int_t getData_DR(const char* fileName,
		   std::vector <std::vector<Float_t> >& ZgRdataMC,
		   std::vector <std::vector<Float_t> >& ZgRdataMCerrUp,
		   std::vector <std::vector<Float_t> >& ZgRdataMCerrLow,
		   Float_t& DRscaleErr,
		   Float_t& DY0bPurErr,
		   Float_t& DYtrigEffErr,
		   Float_t& LeptonSFerr,
		   Float_t& btag_SFerr);
  Int_t getData_DR(const char* fileName,
		   std::vector <std::vector<Float_t> >& ZgRdataMC,
		   std::vector <std::vector<Float_t> >& ZgRdataMCerrUp,
		   std::vector <std::vector<Float_t> >& ZgRdataMCerrLow);
  Int_t getData_DY(const char* fileName,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0stat,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0MCstat,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysUp,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysLow,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysKin,
		   std::vector <std::vector <std::vector<Float_t> > >& Rb0sysPur);

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
    MaxNjets = 4;
    MaxNb = 4;
    MaxKin = 13;
    histoXlabelSize = 0.015;
    canvasBottomMargin = 0.15;
    if (doSample == LDP) {
      dataFile_gJets = gJetsFnRoot+TString("_ldp.dat");
      dataFile_DR = DRfnRoot+TString("_ldp.dat");
      dataFile_DY = DYfnRoot+TString( "_ldp.dat");
      output_rootFile = "../plots/histograms/ZinvHistos_ldp.root";
      output_plotFile = "../plots/png/ZinvBGpred_ldp.png";
    } else if (doSample == HDP) {
      dataFile_gJets = gJetsFnRoot+TString("_hdp.dat");
      dataFile_DR = DRfnRoot+TString("_hdp.dat");
      dataFile_DY = DYfnRoot+TString("_hdp.dat");
      output_rootFile = "../plots/histograms/ZinvHistos_hdp.root";
      output_plotFile = "../plots/png/ZinvBGpred_hdp.png";
    }
  } else {
    MaxNjets = 4;
    MaxNb = 4;
    MaxKin = 10;
    histoXlabelSize = 0.04;
    canvasBottomMargin = 0.35;
    // dataFile_gJets = "gJets_signal.dat";
    dataFile_gJets = gJetsFnRoot+TString("_signal.dat");
    dataFile_DR = DRfnRoot+TString("_signal.dat");
    dataFile_DY = DYfnRoot+TString("_signal.dat");
    output_rootFile = "../plots/histograms/ZinvHistos.root";
    output_plotFile = "../plots/png/ZinvBGpred.png";
  }
  MaxKinDY = expandDYkin ? MaxKin : 1;
  const Int_t MaxBins(MaxNjets*MaxNb*MaxKin);

  // N(b=0) from gamma+jets
  std::vector <std::vector<Float_t> > Ngobs(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > NgobsEB(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > NgobsEE(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgR(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerr(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrUp(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRerrLow(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gPur(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > gPurErr(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMC(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMCerrUp(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMCerrLow(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMC_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMCerrUp_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
  std::vector <std::vector<Float_t> > ZgRdataMCerrLow_from_gJets(MaxNjets, std::vector<Float_t>(MaxKin, 0));
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
  if (0 != getData_gJets(dataFile_gJets, Ngobs, NgobsEB, NgobsEE, ZgR, ZgRerr, ZgRerrUp, ZgRerrLow, gPur, gPurErr,
			 ZgRdataMC_from_gJets, ZgRdataMCerrUp_from_gJets, ZgRdataMCerrLow_from_gJets)) {
    cout << "Failed to get data." << endl;
    return;
  }
  if (0 != getData_DR(dataFile_DR, ZgRdataMC, ZgRdataMCerrUp, ZgRdataMCerrLow, DRscaleErr, DY0bPurErr, DYtrigEffErr, LeptonSFerr, btagSFerr)) {
    cout << "Failed to get data." << endl;
    return;
  }
  if (0 != getData_DY(dataFile_DY, DYvalues, DYstat, DYMCstat, DYsysNjUp, DYsysNjLow, DYsysKin, DYsysPur)) {
    cout << "Failed to get data." << endl;
    return;
  }

  Float_t gPurAv = 0, gPurErrAv = 0, NobsSum = 0, ZgRAv = 0, ZgRerrUpAv = 0, ZgRerrLowAv = 0;
  cout << "DRscaleErr = " << DRscaleErr << ", DY0bPurErr = " << DY0bPurErr << ", DYtrigEffErr = " << DYtrigEffErr << ", LeptonSFerr = " <<LeptonSFerr  << ", btagSFerr = " << btagSFerr << endl;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
      cout << MaxKin*ijet+ikin+1 << "  Ngobs = " << Ngobs[ijet][ikin] << "  ZgR = " << ZgR[ijet][ikin] 
           << "  ZgRerr = " << ZgRerr[ijet][ikin]
           << "  ZgRerrUp = " << ZgRerrUp[ijet][ikin]
           << "  ZgRerrLow = " << ZgRerrLow[ijet][ikin]
           << "  gPur = " << gPur[ijet][ikin]
           << "  gPurErr = " << gPurErr[ijet][ikin]
           << "  ZgRdataMC = " << ZgRdataMC[ijet][ikin]
           << "  ZgRdataMCerrUp = " << ZgRdataMCerrUp[ijet][ikin]
           << "  ZgRdataMCerrLow = " << ZgRdataMCerrLow[ijet][ikin] << endl;
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
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      for (Int_t ikin=0; ikin<MaxKinDY; ++ikin) {
	cout << 4*ijet+ib+1 << "  DYvalue = " << DYvalues[ijet][ib][ikin]
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

  TH1F*hTemplate = new TH1F("hTemplate", "template for Zinv histograms", MaxBins, 0.5, MaxBins+0.5);

  // Tune up the template histogram
  //hTemplate->Reset();
  hTemplate->GetXaxis()->SetLabelSize(histoXlabelSize);
  hTemplate->GetYaxis()->SetLabelSize(0.05);
  //hTemplate->GetYaxis()->SetTitleSize(0.05);
  hTemplate->SetOption("HIST");
  hTemplate->SetMarkerSize(0);
 // Set the x-axis labels
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
	Int_t bin = MaxKin*(MaxNb*ijet + ib) + ikin + 1;
	char label[25];
	Int_t iMHT = ikin/3;
	if(ikin==8) {
	  iMHT=3;
	}
	Int_t iHT = ikin;
	//if (doSample != Signal && iMHT == 3) iHT++;
	//if (doSample == Signal && ikin==5) iMHT = 2;
	//if (doSample == Signal) iHT = ikin;
	if(doSample == Signal)
	  sprintf(label, "%5s%1d%6s%1d%4s%1d%3s%1d", "NJets", ijet, "-BTags", ib, "-MHT", iMHT, "-HT", iHT);
	else {
	  if(ikin==0)
	    sprintf(label, "%5s%1d%6s%1d%11s", "NJets", ijet, "-BTags", ib, "-MHTC1-HTC1");
	  else if(ikin==1)
	    sprintf(label, "%5s%1d%6s%1d%11s", "NJets", ijet, "-BTags", ib, "-MHTC2-HTC2");
	  else if(ikin==2)
	    sprintf(label, "%5s%1d%6s%1d%11s", "NJets", ijet, "-BTags", ib, "-MHTC3-HTC3");
	  else {
	    iMHT = (ikin-3)/3;
	    if((ikin-3)==8) {
	      iMHT=3;
	    }
	    iHT = ikin-3;
	    sprintf(label, "%5s%1d%6s%1d%4s%1d%3s%1d", "NJets", ijet, "-BTags", ib, "-MHT", iMHT, "-HT", iHT);
	  }
	}
	// cout << label << endl;
	hTemplate->GetXaxis()->SetBinLabel(bin, label);
      }
    }
  }


  // Output file
  TFile *ZinvHistos = TFile::Open(output_rootFile, "RECREATE");

  gStyle->SetOptStat("n");

  // Book histograms for N(b=0) yield factors from gamma+jets
  //  NZvv_pred = Ngobs * ZgR * gPur * ZgRdataMC
  TH1F *hgJNobs = (TH1F*)hTemplate->Clone("hgJNobs");
  hgJNobs->GetYaxis()->SetTitle("gamma+jets yield");
  TH1F *hgJstat = (TH1F*)hTemplate->Clone("hgJstat");
  hgJstat->GetYaxis()->SetTitle("gamma+jets yield stat error");
  TH1F *hgJZgR = (TH1F*)hTemplate->Clone("hgJZgR");
  hgJZgR->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio");
  TH1F *hgJZgRerr = (TH1F*)hTemplate->Clone("hgJZgRerr");
  hgJZgRerr->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");
  TH1F *hgJZgRerrUp = (TH1F*)hTemplate->Clone("hgJZgRerrUp");
  hgJZgRerrUp->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");
  TH1F *hgJZgRerrLow = (TH1F*)hTemplate->Clone("hgJZgRerrLow");
  hgJZgRerrLow->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");
  TH1F *hgJPur = (TH1F*)hTemplate->Clone("hgJPur");
  hgJPur->GetYaxis()->SetTitle("gamma+jets purity");
  TH1F *hgJPurErr = (TH1F*)hTemplate->Clone("hgJPurErr");
  hgJPurErr->GetYaxis()->SetTitle("gamma+jets purity error");
  TH1F *hgJZgRdataMC = (TH1F*)hTemplate->Clone("hgJZgRdataMC");
  hgJZgRdataMC->GetYaxis()->SetTitle("gamma+jets double ratio");
  TH1F *hgJZgRdataMCerrUp = (TH1F*)hTemplate->Clone("hgJZgRdataMCerrUp");
  hgJZgRdataMCerrUp->GetYaxis()->SetTitle("gamma+jets double ratio error");
  TH1F *hgJZgRdataMCerrLow = (TH1F*)hTemplate->Clone("hgJZgRdataMCerrLow");
  hgJZgRdataMCerrLow->GetYaxis()->SetTitle("gamma+jets double ratio error");
  TH1F *hZinvScaleErr = (TH1F*)hTemplate->Clone("hZinvScaleErr");
  hZinvScaleErr->GetYaxis()->SetTitle("Zinv global scale error");

  // Book histograms for Nb/N0 ratio values and error components
  //  Nb/N0 pred = DYvalues
  TH1F *hDYvalue = (TH1F*)hTemplate->Clone("hDYvalue");
  hDYvalue->GetYaxis()->SetTitle("DY ratio to 0b value");
  TH1F *hDYstat = (TH1F*)hTemplate->Clone("hDYstat");
  hDYstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");
  TH1F *hDYMCstat = (TH1F*)hTemplate->Clone("hDYMCstat");
  hDYMCstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");
  TH1F *hDYsysNjUp = (TH1F*)hTemplate->Clone("hDYsysNjUp");
  hDYsysNjUp->GetYaxis()->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");
  TH1F *hDYsysNjLow = (TH1F*)hTemplate->Clone("hDYsysNjLow");
  hDYsysNjLow->GetYaxis()->SetTitle("DY ratio to 0b syst- error Nj extrapolation");
  TH1F *hDYsysKin = (TH1F*)hTemplate->Clone("hDYsysKin");
  hDYsysKin->GetYaxis()->SetTitle("DY ratio to 0b syst error kinematics dependence");
  TH1F *hDYsysPur = (TH1F*)hTemplate->Clone("hDYsysPur");
  hDYsysPur->GetYaxis()->SetTitle("DY ratio to 0b syst error purity");

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

  ofstream tableFile;  tableFile.open("../textFiles/table_for_AN.txt");
  char buf[256];

  // Fill histograms
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      for (Int_t ikin=0; ikin<MaxKin; ++ikin) {

	Int_t bin = MaxKin*(MaxNb*ijet + ib) + ikin + 1;
	Int_t ikinDY = MaxKinDY == MaxKin ? ikin : 0;

	if (ib == 0) {
	  hgJNobs->SetBinContent(bin, Ngobs[ijet][ikin]);
	  Ngobs[ijet][ikin] > 0 ? hgJstat->SetBinContent(bin, 1/Sqrt(Ngobs[ijet][ikin])) : hgJstat->SetBinContent(bin, 0);
	  hgJZgR->SetBinContent(bin, ZgR[ijet][ikin]);
	  hgJZgRdataMC->SetBinContent(bin, ZgRdataMC[ijet][ikin]);
	  // Avoid degenerate factors by combining ZgR stat and ZgRdataMC errors into ZgRdataMCErrUp,Low
	  // hgJZgRerr->SetBinContent(bin, ZgRerr[ijet][ikin]);
	  // hgJZgRdataMCerrUp->SetBinContent(bin, ZgRdataMCerrUp[ijet][ikin]);
	  // hgJZgRdataMCerrLow->SetBinContent(bin, ZgRdataMCerrLow[ijet][ikin]);
	  hgJZgRerrUp->SetBinContent(bin, ZgRerrUp[ijet][ikin]);  // correlated
	  hgJZgRerrLow->SetBinContent(bin, ZgRerrLow[ijet][ikin]);
	  hgJZgRerr->SetBinContent(bin, 0);  // uncorrelated, combined into ZgRdataMCErrUp,Low
	  hgJZgRdataMCerrUp->SetBinContent(bin, Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(btagSFerr, 2) + Power(ZgRdataMCerrUp[ijet][ikin], 2)));  // uncorrelated
	  hgJZgRdataMCerrLow->SetBinContent(bin, Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(btagSFerr, 2) + Power(ZgRdataMCerrLow[ijet][ikin], 2)));
	  // hgJZgRerrUp->SetBinContent(bin, Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(ZgRerrUp[ijet][ikin], 2) + Power(ZgRdataMCerrUp[ijet][ikin], 2)));
	  // hgJZgRerrLow->SetBinContent(bin, Sqrt(Power(ZgRerr[ijet][ikin], 2) + Power(ZgRerrLow[ijet][ikin], 2) + Power(ZgRdataMCerrLow[ijet][ikin], 2)));
	  hgJPur->SetBinContent(bin, gPur[ijet][ikin]);
	  hgJPurErr->SetBinContent(bin, gPurErr[ijet][ikin]);  // correlated
	  hZinvScaleErr->SetBinContent(bin, Sqrt(Power(DRscaleErr, 2) + Power(DY0bPurErr, 2) + Power(DYtrigEffErr, 2) + Power(LeptonSFerr, 2)));
	  //
	  // Here we take into account that photon purity and ZgR syst errors partly cancel in the prediction
	  //
	  Bool_t useEffectiveErr = true;
	  Float_t gPurErrEff = fabs( gPurErr[ijet][ikin] - gPurErrAv / gPurAv );  // Frac. error on gPur/<gPur>
	  Float_t ZgRerrUpEff = fabs( ZgRerrUp[ijet][ikin] - ZgRerrUpAv / ZgRAv );  // Frac. upper error on ZgR/<ZgR>
	  Float_t ZgRerrLowEff = fabs( ZgRerrLow[ijet][ikin] - ZgRerrLowAv / ZgRAv );  // Frac. upper error on ZgR/<ZgR>
	  // cout << "gPur err, nominal: " << gPurErr[ijet][ikin] << ", eff: " << gPurErrEff
	  //      << ",    ZgR errUp nominal: " << ZgRerrUp[ijet][ikin] << ", eff: " << ZgRerrUpEff
	  //      << ",    ZgR errLow nominal: " << ZgRerrLow[ijet][ikin] << ", eff: " << ZgRerrLowEff << endl;
	  if (useEffectiveErr) {
	    hgJZgRerrUp->SetBinContent(bin, ZgRerrUpEff);  // correlated
	    hgJZgRerrLow->SetBinContent(bin, ZgRerrLowEff);
	    hgJPurErr->SetBinContent(bin, gPurErrEff);  // correlated
	  }
	} else {
	  hgJNobs->SetBinContent(bin, ignoreBin);
	  hgJstat->SetBinContent(bin, ignoreBin);
	  hgJZgR->SetBinContent(bin, ignoreBin);
	  hgJZgRerr->SetBinContent(bin, ignoreBin);
	  hgJZgRdataMC->SetBinContent(bin, ignoreBin);
	  hgJZgRerrUp->SetBinContent(bin, ignoreBin);
	  hgJZgRerrLow->SetBinContent(bin, ignoreBin);
	  hgJZgRdataMCerrUp->SetBinContent(bin, ignoreBin);
	  hgJZgRdataMCerrLow->SetBinContent(bin, ignoreBin);
	  hgJPur->SetBinContent(bin, ignoreBin);
	  hgJPurErr->SetBinContent(bin, ignoreBin);
	  hZinvScaleErr->SetBinContent(bin, ignoreBin);
	}

  	hDYvalue->SetBinContent(bin, DYvalues[ijet][ib][ikinDY]);
  	hDYstat->SetBinContent(bin, DYstat[ijet][ib][ikinDY]);
  	hDYMCstat->SetBinContent(bin, DYMCstat[ijet][ib][ikinDY]);
  	hDYsysNjUp->SetBinContent(bin, DYsysNjUp[ijet][ib][ikinDY]);
  	hDYsysNjLow->SetBinContent(bin, DYsysNjLow[ijet][ib][ikinDY]);
  	hDYsysKin->SetBinContent(bin, DYsysKin[ijet][ib][ikinDY]);
  	hDYsysPur->SetBinContent(bin, DYsysPur[ijet][ib][ikinDY]);

  	// cout << bin
	//      << " " << hDYvalue->GetXaxis()->GetBinLabel(bin)
	//      << "  " << hDYvalue->GetBinContent(bin) << endl;

	Float_t thisNgobs = Ngobs[ijet][ikin] > 0 ? Ngobs[ijet][ikin] : 1.0;
	Float_t ZinvValue = thisNgobs * ZgR[ijet][ikin] * gPur[ijet][ikin] * DYvalues[ijet][ib][ikinDY];
	ZinvValue *= ZgRdataMC[ijet][ikin];  // Here apply the double ratio
	// cout << Ngobs[ijet][ikin] << " " << ZgR[ijet][ikin] << " " << gPur[ijet][ikin] << " " << DYvalues[ijet][ib][ikin] << " " << ZgRdataMC[ijet][ikin] << endl;
	Float_t wtStat = 1/thisNgobs;
	statErr[bin-1] = ZinvValue*Sqrt(wtStat + Power(DYstat[ijet][ib][ikinDY], 2));
	sysUp[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				      + Power(ZgRerrUp[ijet][ikin], 2)
				      + Power(gPurErr[ijet][ikin], 2)
				      + Power(ZgRdataMCerrUp[ijet][ikin], 2)
				      + Power(DYMCstat[ijet][ib][ikinDY], 2)
				      + Power(DYsysNjUp[ijet][ib][ikinDY], 2)
				      + Power(DYsysKin[ijet][ib][ikinDY], 2)
				      + Power(DYsysPur[ijet][ib][ikinDY], 2));
	sysLow[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				       + Power(ZgRerrLow[ijet][ikin], 2)
				       + Power(gPurErr[ijet][ikin], 2)
				       + Power(ZgRdataMCerrLow[ijet][ikin], 2)
				       + Power(DYMCstat[ijet][ib][ikinDY], 2)
				       + Power(DYsysNjLow[ijet][ib][ikinDY], 2)
				       + Power(DYsysKin[ijet][ib][ikinDY], 2)
				       + Power(DYsysPur[ijet][ib][ikinDY], 2));
	// cout << bin << "  " << ZinvValue << " +/- " << statErr[bin-1] << " + " << sysUp[bin-1]
	//      << " - " << sysLow[bin-1] << endl;
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
	  if (ikin == 0) tableFile << endl;
	  sprintf(buf, "%5.0f & %5.0f & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%6.1f\\pm%4.1f^{+%4.1f}_{-%4.1f}$ \\\\\n",
		  NgobsEB[ijet][ikin], NgobsEE[ijet][ikin],
		  ZgR[ijet][ikin], ZgR[ijet][ikin]*ZgRerr[ijet][ikin], ZgR[ijet][ikin]*ZgRerrUp[ijet][ikin],
		  ZgR[ijet][ikin]*ZgRerrLow[ijet][ikin], 
		  ZgRdataMC[ijet][ikin], DRscaleErr, ZgRdataMC[ijet][ikin]*ZgRdataMCerrUp[ijet][ikin],
		  ZgRdataMC[ijet][ikin]*ZgRdataMCerrLow[ijet][ikin],
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
    if (MCfile != 0) {
      if (!MCfile->Get("plot_zinv_nj4_nb4_kin10_1")) {
	cout << "No MC histogram found in this file and directory" << endl;
      } else {
	hZinvMC = (TH1F*) MCfile->Get("plot_zinv_nj4_nb4_kin10_1");
	hZinvMCbin = (TH1F*)hTemplate->Clone("hZinvMCbin");
	hZinvMCbin->GetYaxis()->SetTitle("MC Z#rightarrow#nu#bar{#nu}");
	for (Int_t bin=1; bin<=160; ++bin) {
	  hZinvMCbin->SetBinContent(bin, hZinvMC->GetBinContent(bin));
	  hZinvMCbin->SetBinError(bin, hZinvMC->GetBinError(bin));
	}
	hZinvMCbin->SetLineColor(kGreen+2);
	hZinvMCbin->Scale(MClumiScale);
      }
    }
  }

  TCanvas *Cvalue = new TCanvas("Cvalue", "DY ratio to zero b", 1000, 500);
  Cvalue->SetTopMargin(.03);
  Cvalue->SetBottomMargin(canvasBottomMargin);
  Cvalue->SetRightMargin(.03);
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
  Cerrors->SetTopMargin(.03);
  Cerrors->SetBottomMargin(canvasBottomMargin);
  Cerrors->SetRightMargin(.03);
  hgJstat->SetMinimum(0);  hgJstat->SetMaximum(1.5);
  hgJstat->SetLineColor(1);  hgJstat->Draw();
  hgJZgRerr->SetLineColor(2);  hgJZgRerr->Draw("same");
  // hgJZgRerrUp->SetLineColor(6);  hgJZgRerrUp->Draw("same");
  // hgJZgRerrLow->SetLineColor(4);  hgJZgRerrLow->Draw("same");
  hgJPurErr->SetLineColor(7);  hgJPurErr->Draw("same");
  hZinvScaleErr->SetLineColor(6);  hZinvScaleErr->Draw("same");
  hgJZgRdataMCerrUp->SetLineColor(3);  hgJZgRdataMCerrUp->Draw("same");
  hgJZgRdataMCerrLow->SetLineColor(5);  hgJZgRdataMCerrLow->Draw("same");
  hDYstat->SetLineColor(8);  hDYstat->Draw("same");
  hDYMCstat->SetLineColor(9);  hDYMCstat->Draw("same");
  hDYsysNjUp->SetLineColor(38);  hDYsysNjUp->Draw("same");
  hDYsysNjLow->SetLineColor(41);  hDYsysNjLow->Draw("same");
  hDYsysKin->SetLineColor(42);  hDYsysKin->Draw("same");
  hDYsysPur->SetLineColor(46);  hDYsysPur->Draw("same");
  TLegend* ErrorsLegend;
  if (doSample == Signal)
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  else
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  ErrorsLegend->AddEntry(hgJstat, "gJstat");
  ErrorsLegend->AddEntry(hgJZgRerr, "RZgErr");
  // ErrorsLegend->AddEntry(hgJZgRerrUp, "RZgU");
  // ErrorsLegend->AddEntry(hgJZgRerrLow, "RZgLow");
  ErrorsLegend->AddEntry(hgJPurErr, "gJpur");
  ErrorsLegend->AddEntry(hZinvScaleErr, "Scale");
  ErrorsLegend->AddEntry(hgJZgRdataMCerrUp, "DRup");
  ErrorsLegend->AddEntry(hgJZgRdataMCerrLow, "DRlow");
  ErrorsLegend->AddEntry(hDYstat, "DYstat");
  ErrorsLegend->AddEntry(hDYMCstat, "DYMCstat");
  ErrorsLegend->AddEntry(hDYsysNjUp, "DYNjUp");
  ErrorsLegend->AddEntry(hDYsysNjLow, "DYNjLow");
  ErrorsLegend->AddEntry(hDYsysKin, "DYkin");
  ErrorsLegend->AddEntry(hDYsysPur, "DYpur");
  ErrorsLegend->Draw();
  Cerrors->Update();
  if (doSample == Signal)  Cerrors->SaveAs("../plots/png/errorComponents.png");


}  // ---------------------------------------------------------------

Int_t getData_gJets(const char* fileName,
		    std::vector <std::vector<Float_t> >& Ngobs,
		    std::vector <std::vector<Float_t> >& NgobsEB,
		    std::vector <std::vector<Float_t> >& NgobsEE,
		    std::vector <std::vector<Float_t> >& ZgR,
		    std::vector <std::vector<Float_t> >& ZgRerr,
		    std::vector <std::vector<Float_t> >& ZgRerrUp,
		    std::vector <std::vector<Float_t> >& ZgRerrLow,
		    std::vector <std::vector<Float_t> >& gPur,
		    std::vector <std::vector<Float_t> >& gPurErr,
		    std::vector <std::vector<Float_t> >& ZgRdataMC,
		    std::vector <std::vector<Float_t> >& ZgRdataMCerrUp,
		    std::vector <std::vector<Float_t> >& ZgRdataMCerrLow) {
  //
  // Read data from the file, parse, and set data arrays.
  //
  /*
  Notes on correlations:
  gPurErr is correlated among the 40 bins for which Nb = 0.
  ZgRerrUp,Low are correlated among the 18 bins.
  ZgRerr is uncorrelated among the 40 bins.
  ZgRdataMCerrUp, Low are uncorrelated among the 40 bins.
  As of 10 Jul 2016 we combine these last two in ZgRdataMCerrUp, Low.
   */
  Float_t gPurErrIn;
  Float_t gFdirErr = 0.076 * 0;  //  This hard-wired number should be handled in a more robust way
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
  while (strncmp(buf, "no ", 3) != 0);
  // while (strncmp(buf, "no  Nobs", 8) != 0);
  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ikin=0; ikin<Ncol; ++ikin) {
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
      n++; token[n] = strtok(0, ",");             // RZgErr
      sscanf(token[n], "%f", &ZgRerr[ijet][ikin]);
      n++; token[n] = strtok(0, "-");             // RZgErrUp
      sscanf(token[n], "%f", &ZgRerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // RZgErrLow
      sscanf(token[n], "%f", &ZgRerrLow[ijet][ikin]);
      n++; token[n] = strtok(0, ")");
      n++; token[n] = strtok(0, "|");

      n++; token[n] = strtok(0, "(");             // Purity
      sscanf(token[n], "%f", &gPur[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // PurityErr
      // sscanf(token[n], "%f", &gPurErr[ijet][ikin]);
      sscanf(token[n], "%f", &gPurErrIn);
      gPurErr[ijet][ikin] = Sqrt(Power(gPurErrIn, 2) - Power(gFdirErr, 2));
      // printf("%s%6.3f%6.3f\n", "Purity error, modified = ", gPurErrIn, gPurErr[ijet][ikin]);
      n++; token[n] = strtok(0, " ");
      n++; token[n] = strtok(0, "(");             // ZgRdataMC
      sscanf(token[n], "%f", &ZgRdataMC[ijet][ikin]);
      n++; token[n] = strtok(0, "-");             // ZgRdataMCErrUp
      sscanf(token[n], "%f", &ZgRdataMCerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, ")");             // ZgRdataMCErrLow
      sscanf(token[n], "%f", &ZgRdataMCerrLow[ijet][ikin]);
    }
  }
  return 0;
}  // ---------------------------------------------------------------


Int_t getData_DR(const char* fileName,
		 std::vector <std::vector<Float_t> >& ZgRdataMC,
		 std::vector <std::vector<Float_t> >& ZgRdataMCerrUp,
		 std::vector <std::vector<Float_t> >& ZgRdataMCerrLow,
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
  gPurErr is correlated among the 40 bins for which Nb = 0.
  ZgRerrUp,Low are correlated among the 18 bins.
  ZgRerr is uncorrelated among the 40 bins.
  ZgRdataMCerrUp, Low are uncorrelated among the 40 bins.
  As of 10 Jul 2016 we combine these last two in ZgRdataMCerrUp, Low.
   */
  Int_t Nrow = ZgRdataMC.size();
  Int_t Ncol = ZgRdataMC[0].size();
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
      token[n] = strtok(buf, "|");  // ZgRdataMC
      if (token[n] == 0) {
	cout << n << endl;
	cout << "DR Empty line encountered" << endl;
	return 1;
      }
      sscanf(token[n], "%f", &ZgRdataMC[ijet][ikin]);
      n++; token[n] = strtok(0, "|");  // Double ratio scale error
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DRscaleErr);
      n++; token[n] = strtok(0, "|");             // ZgRdataMCErrUp
      sscanf(token[n], "%f", &ZgRdataMCerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // ZgRdataMCErrLow
      sscanf(token[n], "%f", &ZgRdataMCerrLow[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // 0b purity
      n++; token[n] = strtok(0, "|");             // 0b purity error
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DY0bPurErr);
      n++; token[n] = strtok(0, "|");             // 0b trigger eff.
      n++; token[n] = strtok(0, "|");             // 0b trigger eff. error
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DYtrigEffErr);
      n++; token[n] = strtok(0, "|");             // lepton SF
      n++; token[n] = strtok(0, "|");             // lepton SF error
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &LeptonSFerr);
      n++; token[n] = strtok(0, "|");             // btag SF
      n++; token[n] = strtok(0, "|");             // btag SF error
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
  DYstat is 3 independent errors, each common to bins with 
    Nb = 1, 2, >=3, respectively.
  DYsysPur is 2 independent errors, one common to bins with Nb = 1,
    the other common to bins with Nb > 1.
  DYsysNjUp, Low are correlated over the 6 bins Njets = 7-8, >=9; Nb = 1, 2, >=3;
    for the lowest Njets bin they vanish.
  The other errors are uncorrelated (though all are zero for Nb = 0).
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
      for (Int_t ikin=0; ikin<Nplane; ++ikin) {
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
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0stat[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0MCstat[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysUp[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysLow[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysKin[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysPur[ijet][ib][ikin]);
      }
    }
  }
  return 0;
}  // ---------------------------------------------------------------

