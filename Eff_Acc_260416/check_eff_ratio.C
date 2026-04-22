void check_eff_ratio() {
  TFile *f_new = TFile::Open("./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260114_test.root");
  TFile *f_nocut = TFile::Open("./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_260108.root");

  TH1D *h_new_for = (TH1D*)f_new->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");
  TH1D *h_nocut_for = (TH1D*)f_nocut->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");

  cout << "=== Forward Rapidity (1.6-2.4) Efficiency ===" << endl;
  cout << "pT Bin\t\tNo Cut\t\tWith Cut\tRatio" << endl;
  for(int i=1; i<=h_new_for->GetNbinsX(); i++) {
    double eff_nocut = h_nocut_for->GetBinContent(i);
    double eff_new = h_new_for->GetBinContent(i);
    double ratio = (eff_nocut > 0) ? eff_new/eff_nocut : 0;
    cout << h_new_for->GetXaxis()->GetBinLowEdge(i) << "-" << h_new_for->GetXaxis()->GetBinUpEdge(i) 
         << "\t" << Form("%.4f", eff_nocut) << "\t" << Form("%.4f", eff_new) 
         << "\t" << Form("%.4f (%.1f%%)", ratio, ratio*100) << endl;
  }

  cout << "\n=== Mid Rapidity (0-1.6) Efficiency ===" << endl;
  TH1D *h_new_mid = (TH1D*)f_new->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6  ");
  TH1D *h_nocut_mid = (TH1D*)f_nocut->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6  ");

  cout << "pT Bin\t\tNo Cut\t\tWith Cut\tRatio" << endl;
  for(int i=1; i<=h_new_mid->GetNbinsX(); i++) {
    double eff_nocut = h_nocut_mid->GetBinContent(i);
    double eff_new = h_new_mid->GetBinContent(i);
    double ratio = (eff_nocut > 0) ? eff_new/eff_nocut : 0;
    cout << h_new_mid->GetXaxis()->GetBinLowEdge(i) << "-" << h_new_mid->GetXaxis()->GetBinUpEdge(i) 
         << "\t" << Form("%.4f", eff_nocut) << "\t" << Form("%.4f", eff_new) 
         << "\t" << Form("%.4f (%.1f%%)", ratio, ratio*100) << endl;
  }
}
