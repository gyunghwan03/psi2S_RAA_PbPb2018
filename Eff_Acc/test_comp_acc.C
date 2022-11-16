void test_comp_acc(){
  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptFit(0);
  TFile *f1 = new TFile("acceptance_Prompt_psi2s_GenOnly_wgt1_20211218_SysUp0_20211218.root","READ");
  TFile *f2 = new TFile("acceptance_Prompt_psi2s_GenOnly_wgt1_20211218_SysUp1_20211218.root","READ");
  TFile *f3 = new TFile("acceptance_Prompt_psi2s_GenOnly_wgt1_20211218_SysUp2_20211218.root","READ");

  TH1F *h1 = (TH1F*)f1->Get("hAccPt_2021_ally");
  TH1F *h2 = (TH1F*)f2->Get("hAccPt_2021_ally");
  TH1F *h3 = (TH1F*)f3->Get("hAccPt_2021_ally");
  TH1F *h4 = (TH1F*)f1->Get("hAccPt_2021_Fory");
  TH1F *h5 = (TH1F*)f2->Get("hAccPt_2021_Fory");
  TH1F *h6 = (TH1F*)f3->Get("hAccPt_2021_Fory");

  h1->SetMarkerColor(kBlue);
  h2->SetMarkerColor(kRed);
  h3->SetMarkerColor(kViolet);
  h2->SetMarkerStyle(24);
  h3->SetMarkerStyle(25);


  h4->SetMarkerColor(kBlue);
  h5->SetMarkerColor(kRed);
  h6->SetMarkerColor(kViolet);
  h5->SetMarkerStyle(24);
  h6->SetMarkerStyle(25);

  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->Divide(2,1);
  c1->cd(1);

  TLegend *l = new TLegend(0.60,0.50,0.89,0.62);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetMargin(0.2);
  l->SetTextSize(0.029);
  l->SetTextFont(42);
  l->SetHeader("|y| < 2.4");
  l->AddEntry(h1,"Nominal","PL");
  l->AddEntry(h2,"Up","PL");
  l->AddEntry(h3,"Down","PL");

  TLegend *l1 = new TLegend(0.60,0.50,0.89,0.62);
  l1->SetFillColor(0);
  l1->SetBorderSize(0);
  l1->SetMargin(0.2);
  l1->SetTextSize(0.029);
  l1->SetTextFont(42);
  l1->SetHeader("1.6 < |y| < 2.4");
  l1->AddEntry(h1,"Nominal","PL");
  l1->AddEntry(h2,"Up","PL");
  l1->AddEntry(h3,"Down","PL");


  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle("Acceptance");
  h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h1->GetYaxis()->SetRangeUser(0,1.0);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");

  l->Draw("same");

  c1->cd(2);
  h4->GetYaxis()->CenterTitle();
  h4->GetXaxis()->CenterTitle();
  h4->GetYaxis()->SetTitle("Acceptance");
  h4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h4->GetYaxis()->SetRangeUser(0,1.0);
  h4->Draw();
  h5->Draw("same");
  h6->Draw("same");

  l1->Draw("same");


  c1->SaveAs("plot_acceptance_comp.png");

}
