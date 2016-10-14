{
  gROOT->Reset();
  gROOT->ProcessLine(".L RA2bin_inputs_Zinv.C+");
  // gROOT->ProcessLine(".L RA2bin_inputs_Zinv_reclus.C+");

  RA2bin_inputs_Zinv(LDP, "../datFiles/gJets_pu", "../datFiles/DR_pu", "../datFiles/DY_pu", "../plots/histograms/ZinvMCttzMC160bin.root", 1.0, true);
  RA2bin_inputs_Zinv(HDP, "../datFiles/gJets_pu", "../datFiles/DR_pu", "../datFiles/DY_pu", "../datFiles/ZinvMCttzMC160bin.root", 1.0, true);
  RA2bin_inputs_Zinv(Signal, "../datFiles/gJets_pu", "../datFiles/DR_pu", "../datFiles/DY_pu", "../plots/histograms/ZinvMCttzMC160bin.root", 1.0, true);

  }
