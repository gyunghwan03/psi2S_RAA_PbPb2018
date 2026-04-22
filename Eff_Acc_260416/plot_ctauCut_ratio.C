#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

void plot_ctauCut_ratio()
{
  const TString fNoCutPath = "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_260123.root";
  const TString fCutPath   = "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260123.root";

  TFile *fNoCut = TFile::Open(fNoCutPath, "READ");
  TFile *fCut   = TFile::Open(fCutPath, "READ");
  if (!fNoCut || fNoCut->IsZombie()) {
    std::cerr << "Cannot open file: " << fNoCutPath << std::endl;
    return;
  }
  if (!fCut || fCut->IsZombie()) {
    std::cerr << "Cannot open file: " << fCutPath << std::endl;
    return;
  }

  const char* hnames[4] = {"hpt_reco_1", "hpt_reco_2", "hcent_reco_1", "hcent_reco_2"};

  for (int i = 0; i < 4; ++i) {
    TH1 *hNoCut = dynamic_cast<TH1*>(fNoCut->Get(hnames[i]));
    TH1 *hCut   = dynamic_cast<TH1*>(fCut->Get(hnames[i]));

    if (!hNoCut || !hCut) {
      std::cerr << "Missing histogram: " << hnames[i] << std::endl;
      continue;
    }

    TH1 *hRatio = dynamic_cast<TH1*>(hCut->Clone(Form("%s_ratio", hnames[i])));
    hRatio->SetDirectory(0);
    hRatio->Divide(hNoCut);
    hRatio->SetTitle(Form("%s (ctauCut / noCut)", hnames[i]));
    hRatio->GetYaxis()->SetTitle("ctauCut / noCut");
    hRatio->GetYaxis()->SetRangeUser(0.0, 1.2);
    hRatio->GetXaxis()->SetTitle(hCut->GetXaxis()->GetTitle());

    TCanvas *c = new TCanvas(Form("c_%s", hnames[i]), hnames[i], 700, 600);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.0);
    hRatio->SetLineWidth(2);
    hRatio->Draw("E1");

    c->SaveAs(Form("./figs_ratio_%s.pdf", hnames[i]));
  }

  fNoCut->Close();
  fCut->Close();
}
