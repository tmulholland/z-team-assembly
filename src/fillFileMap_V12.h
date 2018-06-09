
void
RA2bZinvAnalysis::fillFileMap() {
  TString dir(treeLoc_);

  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-100to200.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-200to400.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-400to600.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-600to800.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-800to1200.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-1200to2500.root"));
  fileMap_["zinv"].push_back(dir+TString("/tree_signal/tree_ZJetsToNuNu_HT-2500toInf.root"));

  fileMap_["ttzvv"].push_back(dir+TString("/tree_signal/tree_TTZToLLNuNu.root"));

}
