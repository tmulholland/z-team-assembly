{
  /*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  */

  bool doCCzvv = false;
  bool doCCttzvv = false;
  bool do1Dzvv = true;
  bool do1Dttzvv = true;

  RA2bZinvAnalysis* analyzer = new RA2bZinvAnalysis();

  if (doCCzvv || doCCttzvv) {
    // Output file
    TFile *CChistos = TFile::Open("hCCzinv.root", "RECREATE");

    if (doCCzvv) {
      TH1F* hCCzvv = analyzer->makeCChist("zinv");
      hCCzvv->SetName("hCCzvv");  hCCzvv->Draw();
    }
    if (doCCttzvv) {
      TH1F* hCCttzvv = analyzer->makeCChist("ttzvv");
      hCCttzvv->SetName("hCCttzvv");  hCCttzvv->Draw();
    }
    if (doCCzvv && doCCttzvv) {
      TH1F* hCCzinv = (TH1F*) hCCzvv->Clone();
      hCCzinv->Add(hCCttzvv);
      hCCzinv->SetName("hCCzinv");  hCCzinv->Draw();
    }
    CChistos->Write();
  }

  if (do1Dzvv || do1Dttzvv) {
    TFile *histos1D = TFile::Open("histsZjets.root", "RECREATE");
    if (do1Dzvv) {
      std::vector<TH1F*> h_zinv = analyzer->makeHistograms("zinv");
      for (auto& theHist : h_zinv) theHist->Draw();
    }
    if (do1Dttzvv) {
      std::vector<TH1F*> h_ttzvv = analyzer->makeHistograms("ttzvv");
      for (auto& theHist : h_ttzvv) theHist->Draw();
    }
    histos1D->Write();
  }

  // TH1F* hTrueNint = analyzer->makeCChist("zinv");
  // TH1F* hTrueNint_ttzvv = analyzer->makeCChist("ttzvv");
  // hTrueNint->Add(hTrueNint_ttzvv);
  // hTrueNint->SaveAs("hTrueNint.root");
  // TH1F* hmadHT = analyzer->makeCChist("zinv");
  // hmadHT->SaveAs("hmadHT.root");
  // TH2F* hKinBinvsMHT = analyzer->make2Dhist("zinv");
  //  hKinBinvsMHT->SaveAs(" hKinBinvsMHT.root");

  gApplication->Terminate(0);

}
