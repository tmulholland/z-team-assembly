{
  gROOT->Reset();
  gROOT->ProcessLine(".L RA2bin_inputs_Zinv.C+");

  RA2bin_inputs_Zinv(LDP, "../datFiles/gJets", "../datFiles/DR", "../datFiles/DY", "../plots/histograms/ZinvMCttzMC160bin.root", 1.0, true);
  RA2bin_inputs_Zinv(HDP, "../datFiles/gJets", "../datFiles/DR", "../datFiles/DY", "../plots/histograms/ZinvMCttzMC160bin.root", 1.0, true);
  RA2bin_inputs_Zinv(Signal, "../datFiles/gJets", "../datFiles/DR", "../datFiles/DY", "../plots/histograms/ZinvMCttzMC174bin.root", 1.0, true);

  
}
