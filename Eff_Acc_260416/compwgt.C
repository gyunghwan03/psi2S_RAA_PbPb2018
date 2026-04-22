void compwgt(){
  gROOT->Macro("~/rootlogon.C");

  TFile *inn1 = new TFile("ratioDataMC_AA_Jpsi_DATA_y0.0-2.4_210915.root","READ");
  TFile *inn2 = new TFile("ratioDataMC_AA_Jpsi_DATA_y1.6-2.4_210915.root","READ");
  TFile *ino1 = new TFile("ratioDataMC_AA_PromptJpsi_DATA_All_y.root","READ");
  TFile *ino2 = new TFile("ratioDataMC_AA_PromptJpsi_DATA_Forward_y.root","READ");

  TF1 *fn1 = (TF1*)inn1->Get("dataMC_Ratio1");
  TF1 *fn2 = (TF1*)inn2->Get("dataMC_Ratio1");
  TF1 *fo1 = (TF1*)ino1->Get("dataMC_Ratio1");
  TF1 *fo2 = (TF1*)ino2->Get("dataMC_Ratio1");

  fn1->SetLineColor(kRed);
  fn2->SetLineColor(kBlue);
  fo1->SetLineColor(kRed);
  fo2->SetLineColor(kBlue);


  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->Divide(2,1);
  c1->cd(1);
  fo1->Draw();
  fn1->Draw("same");


  c1->cd(2);
  fo2->Draw();
  fn2->Draw("same");
}
