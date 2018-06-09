{
  /*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  */

  RA2bZinvAnalysis* analyzer = new RA2bZinvAnalysis();

  // Output file
  TFile *ZinvHistos = TFile::Open("hCCzinv.root", "RECREATE");

  TH1F* hCCzvv = analyzer->makeCChist("zinv");
  hCCzvv->SetName("hCCzvv");  hCCzvv->Draw();
  TH1F* hCCttzvv = analyzer->makeCChist("ttzvv");
  hCCttzvv->SetName("hCCttzvv");  hCCttzvv->Draw();
  TH1F* hCCzinv = (TH1F*) hCCzvv->Clone();
  hCCzinv->Add(hCCttzvv);
  hCCzinv->SetName("hCCzinv");  hCCzinv->Draw();

  ZinvHistos->Write();

  gApplication->Terminate(0);

}
