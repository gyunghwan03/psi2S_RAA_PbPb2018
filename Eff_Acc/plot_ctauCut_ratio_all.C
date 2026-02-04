#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>

struct SampleCfg {
  TString tag;
  TString fNoCut;
  TString fCut;
};

void DrawRatioSet(const SampleCfg &cfg, const char* hname)
{
  TFile *fNoCut = TFile::Open(cfg.fNoCut, "READ");
  TFile *fCut   = TFile::Open(cfg.fCut, "READ");
  if (!fNoCut || fNoCut->IsZombie()) {
    std::cerr << "Cannot open file: " << cfg.fNoCut << std::endl;
    return;
  }
  if (!fCut || fCut->IsZombie()) {
    std::cerr << "Cannot open file: " << cfg.fCut << std::endl;
    fNoCut->Close();
    return;
  }

  TH1 *hNoCut = dynamic_cast<TH1*>(fNoCut->Get(hname));
  TH1 *hCut   = dynamic_cast<TH1*>(fCut->Get(hname));
  if (!hNoCut || !hCut) {
    std::cerr << "Missing histogram: " << hname << " in " << cfg.tag << std::endl;
    fNoCut->Close();
    fCut->Close();
    return;
  }

  TH1 *hRatio = dynamic_cast<TH1*>(hCut->Clone(Form("%s_%s_ratio", cfg.tag.Data(), hname)));
  hRatio->SetDirectory(0);
  hRatio->Divide(hNoCut);
  hRatio->SetTitle(Form("%s %s (ctauCut / noCut)", cfg.tag.Data(), hname));
  hRatio->GetYaxis()->SetTitle("ctauCut / noCut");
  hRatio->GetXaxis()->SetTitle(hCut->GetXaxis()->GetTitle());
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.0);
  hRatio->SetLineWidth(2);

  TCanvas *c = new TCanvas(Form("c_%s_%s", cfg.tag.Data(), hname), hname, 700, 600);
  hRatio->Draw("E1");
  c->SaveAs(Form("./figs_ctauCut_ratio/%s_%s_ratio.pdf", cfg.tag.Data(), hname));

  fNoCut->Close();
  fCut->Close();
}

void plot_ctauCut_ratio_all()
{
  gSystem->mkdir("./figs_ctauCut_ratio", kTRUE);

  SampleCfg samples[] = {
    {"pbpb_Jpsi",  "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_260123.root",
                    "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260203.root"},
    {"pp_Jpsi",    "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_251104.root",
                    "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260126.root"},
    {"pbpb_psi2S", "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_ctauCut_251103.root",
                    "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_ctauCut_260126.root"},
    {"pp_psi2S",   "./roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root",
                    "./roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_ctauCut_260126.root"}
  };

  const char* hnames[4] = {"hpt_reco_1", "hpt_reco_2", "hcent_reco_1", "hcent_reco_2"};

  for (const auto &cfg : samples) {
    for (int i = 0; i < 4; ++i) {
      DrawRatioSet(cfg, hnames[i]);
    }
  }
}
