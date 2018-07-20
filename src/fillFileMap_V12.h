
void
RA2bZinvAnalysis::fillFileMap() {
  TString dir(treeLoc_);

  // Z->vv MC
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-100to200.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-200to400.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-400to600.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-600to800.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-800to1200.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-1200to2500.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-2500toInf.root"));

  // Z->ttzvv MC
  fileMap_["ttzvv"].push_back(dir+TString("/tree_signal/tree_TTZToLLNuNu.root"));

  // DY->mumu MC
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-600to800.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-800to1200.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-1200to2500.root"));
  fileMap_["dymm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-2500toInf.root"));

  //DY->ee MC
  fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root"));
  fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root"));
  fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root"));
  /* fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-600to800.root")); */
  fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-800to1200.root"));
  /* fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-1200to2500.root")); */
  fileMap_["dyee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-2500toInf.root"));

  // DY->mumu data
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016B.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016C.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016D.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016E.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016F.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016G.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016H2.root"));
  fileMap_["zmm"].push_back(dir+TString("/tree_DYm_CleanVars/tree_SingleMuon_re2016H3.root"));

  // DY->ee data
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016B.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016C.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016D.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016E.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016F.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016G.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016H2.root"));
  fileMap_["zee"].push_back(dir+TString("/tree_DYe_CleanVars/tree_SingleElectron_re2016H3.root"));

}
