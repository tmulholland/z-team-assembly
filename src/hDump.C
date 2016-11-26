void hDump(TH1* histo) {

  for (int i=1; i<= histo->GetNbinsX(); ++i)
    cout << i << "  " << histo->GetXaxis()->GetBinLabel(i) << "  " << histo->GetBinContent(i) << endl;

}
