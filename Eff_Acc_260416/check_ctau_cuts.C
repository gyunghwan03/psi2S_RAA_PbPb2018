void check_ctau_cuts() {
  TFile *f = TFile::Open("roots_ctau3D_v2/ctau3D_cut_ptBin_PRMC_cent0-180.root", "READ");
  if(f && find /data -name>IsZombie()) {
    f->ls();
    TH1D *h_for = (TH1D*)f->Get("h_ctau_ptBin_for");
    TH1D *h_mid = (TH1D*)f->Get("h_ctau_ptBin_mid");
    if(h_for) {
      std::cout << "\n=== Forward pt bins ctau3D cuts ===" << std::endl;
      for(int i=1; i<=h_for->GetNbinsX(); i++) {
        std::cout << "Bin " << i << " (pt " << h_for->GetBinLowEdge(i) << "-" << h_for->GetBinLowEdge(i+1) << "): " << h_for->GetBinContent(i) << std::endl;
      }
    }
    if(h_mid) {
      std::cout << "\n=== Mid pt bins ctau3D cuts ===" << std::endl;
      for(int i=1; i<=h_mid->GetNbinsX(); i++) {
        std::cout << "Bin " << i << " (pt " << h_mid->GetBinLowEdge(i) << "-" << h_mid->GetBinLowEdge(i+1) << "): " << h_mid->GetBinContent(i) << std::endl;
      }
    }
    f->Close();
  }
}
