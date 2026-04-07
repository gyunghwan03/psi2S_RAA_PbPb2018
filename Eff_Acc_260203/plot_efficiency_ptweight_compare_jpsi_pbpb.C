#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../commonUtility.h"

struct PlotItem {
  TString outTag;
  TString histPtW;
  TString histNoPtW;
  TString xTitle;
  TString rapidityText;
  double xMin;
  double xMax;
};

namespace {

const TString kOutDir = "plot_outputs_eff_ptw_compare_jpsi";

void EnsureOutDir() {
  if (gSystem) gSystem->mkdir(kOutDir, true);
}

TString OutPath(const TString& name, const TString& ext) {
  return Form("%s/%s.%s", kOutDir.Data(), name.Data(), ext.Data());
}

TH1* LoadHistWithFallback(TFile* f, const TString& baseName, const TString& cloneName) {
  if (!f) return nullptr;

  std::vector<TString> candidates = {
      baseName,
      baseName + " ",
      baseName + "  "
  };

  for (const auto& key : candidates) {
    TH1* h = dynamic_cast<TH1*>(f->Get(key));
    if (!h) continue;
    TH1* out = dynamic_cast<TH1*>(h->Clone(cloneName));
    if (!out) return nullptr;
    out->SetDirectory(nullptr);
    return out;
  }

  std::cerr << "[ERROR] histogram not found: " << baseName << std::endl;
  return nullptr;
}

TH1* RescaleXAxis(const TH1* h, double scale, const TString& nameSuffix) {
  if (!h) return nullptr;

  const int nBins = h->GetNbinsX();
  const TArrayD* xbins = h->GetXaxis()->GetXbins();
  TH1D* hScaled = nullptr;

  if (xbins && xbins->GetSize() == nBins + 1) {
    std::vector<double> newEdges(nBins + 1, 0.0);
    for (int i = 0; i <= nBins; ++i) {
      newEdges[i] = xbins->GetAt(i) * scale;
    }
    hScaled = new TH1D(Form("%s_%s", h->GetName(), nameSuffix.Data()), "", nBins, newEdges.data());
  } else {
    const double xMin = h->GetXaxis()->GetXmin() * scale;
    const double xMax = h->GetXaxis()->GetXmax() * scale;
    hScaled = new TH1D(Form("%s_%s", h->GetName(), nameSuffix.Data()), "", nBins, xMin, xMax);
  }

  if (!hScaled) return nullptr;
  hScaled->SetDirectory(nullptr);
  hScaled->Sumw2();
  for (int i = 0; i <= nBins + 1; ++i) {
    hScaled->SetBinContent(i, h->GetBinContent(i));
    hScaled->SetBinError(i, h->GetBinError(i));
  }
  return hScaled;
}

void DrawSingleCompare(TFile* fPtW,
                       TFile* fNoPtW,
                       const PlotItem& cfg,
                       const TString& sampleLabel,
                       bool useCentAxis) {
  TH1* hPtW = LoadHistWithFallback(fPtW, cfg.histPtW, Form("h_ptw_%s", cfg.outTag.Data()));
  TH1* hNoPtW = LoadHistWithFallback(fNoPtW, cfg.histNoPtW, Form("h_noptw_%s", cfg.outTag.Data()));
  if (!hPtW || !hNoPtW) {
    delete hPtW;
    delete hNoPtW;
    return;
  }

  if (useCentAxis) {
    TH1* hPtWScaled = RescaleXAxis(hPtW, 0.5, "centScaled_ptw");
    TH1* hNoPtWScaled = RescaleXAxis(hNoPtW, 0.5, "centScaled_noptw");
    if (hPtWScaled && hNoPtWScaled) {
      delete hPtW;
      delete hNoPtW;
      hPtW = hPtWScaled;
      hNoPtW = hNoPtWScaled;
    } else {
      delete hPtWScaled;
      delete hNoPtWScaled;
    }
  }

  TCanvas* c = new TCanvas(Form("c_%s", cfg.outTag.Data()), "", 800, 700);
  TPad* padUp = new TPad(Form("padUp_%s", cfg.outTag.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padDn = new TPad(Form("padDn_%s", cfg.outTag.Data()), "", 0.0, 0.0, 1.0, 0.28);
  padUp->SetBottomMargin(0.03);
  padUp->SetTicks(1, 1);
  padDn->SetTopMargin(0.02);
  padDn->SetBottomMargin(0.33);
  padDn->SetTicks(1, 1);
  padUp->Draw();
  padDn->Draw();

  padUp->cd();
  hPtW->SetLineColor(kBlue + 1);
  hPtW->SetMarkerColor(kBlue + 1);
  hPtW->SetMarkerStyle(20);
  hPtW->SetLineWidth(2);
  hNoPtW->SetLineColor(kRed + 1);
  hNoPtW->SetMarkerColor(kRed + 1);
  hNoPtW->SetMarkerStyle(21);
  hNoPtW->SetLineWidth(2);

  hPtW->SetTitle("");
  hPtW->GetYaxis()->SetTitle("Efficiency");
  hPtW->GetYaxis()->SetRangeUser(0.0, 1.2);
  hPtW->GetXaxis()->SetTitle(cfg.xTitle);
  hPtW->GetXaxis()->SetTitleSize(0.0);
  hPtW->GetXaxis()->SetLabelSize(0.0);

  if (cfg.xMin < cfg.xMax) {
    hPtW->GetXaxis()->SetRangeUser(cfg.xMin, cfg.xMax);
    hNoPtW->GetXaxis()->SetRangeUser(cfg.xMin, cfg.xMax);
  }

  hPtW->Draw("PE");
  hNoPtW->Draw("PE SAME");

  TLegend* leg = new TLegend(0.44, 0.68, 0.84, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->AddEntry(hPtW, "With p_{T} weight", "lep");
  leg->AddEntry(hNoPtW, "No p_{T} weight", "lep");
  leg->Draw();

  drawText(sampleLabel.Data(), 0.20, 0.84, kBlack, 18);
  drawText(cfg.rapidityText.Data(), 0.20, 0.77, kBlack, 18);
  if (!useCentAxis) drawText("Cent. 0-90%", 0.20, 0.70, kBlack, 18);

  padDn->cd();
  TH1* hRatio = dynamic_cast<TH1*>(hNoPtW->Clone(Form("h_ratio_%s", cfg.outTag.Data())));
  hRatio->SetDirectory(nullptr);
  hRatio->Divide(hPtW);
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitle(cfg.xTitle);
  hRatio->GetYaxis()->SetTitle("NoPtW/PtW");
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
  hRatio->GetXaxis()->SetLabelSize(0.11);
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetTitleSize(0.10);
  hRatio->GetYaxis()->SetTitleOffset(0.42);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(21);
  hRatio->Draw("PE");

  TLine* l1 = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  l1->SetLineStyle(2);
  l1->SetLineColor(kGray + 2);
  l1->Draw("same");

  c->SaveAs(OutPath(cfg.outTag, "png"));
  c->SaveAs(OutPath(cfg.outTag, "pdf"));

  delete hRatio;
  delete l1;
  delete leg;
  delete hPtW;
  delete hNoPtW;
  delete c;
}

}  // namespace

void plot_efficiency_ptweight_compare_jpsi_pbpb(int state = 1) {
  gStyle->SetOptStat(0);
  EnsureOutDir();

  const bool isPrompt = (state == 1);
  const TString sample = isPrompt ? "Prompt J/#psi PbPb" : "Nonprompt J/#psi PbPb";

  const TString filePtW = isPrompt
                              ? "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310_2exp.root"
                              : "./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310_2exp.root";

  const TString fileNoPtW = isPrompt
                                ? "./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW1_tnp1_ctauCut_260310_noPtW.root"
                                : "./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW1_tnp1_ctauCut_260310_noPtW.root";

  TFile* fPtW = TFile::Open(filePtW, "READ");
  TFile* fNoPtW = TFile::Open(fileNoPtW, "READ");
  if (!fPtW || fPtW->IsZombie() || !fNoPtW || fNoPtW->IsZombie()) {
    std::cerr << "[ERROR] failed to open input files" << std::endl;
    std::cerr << "  PtW   : " << filePtW << std::endl;
    std::cerr << "  noPtW : " << fileNoPtW << std::endl;
    if (fPtW) {
      fPtW->Close();
      delete fPtW;
    }
    if (fNoPtW) {
      fNoPtW->Close();
      delete fNoPtW;
    }
    return;
  }

  std::vector<PlotItem> jobs = {
      {
          Form("ptw_compare_%s_pt_fwd", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
          "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
          "p_{T} (GeV/c)",
          "1.6 < |y| < 2.4",
          0.0,
          40.0
      },
      {
          Form("ptw_compare_%s_pt_mid", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
          "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
          "p_{T} (GeV/c)",
          "|y| < 1.6",
          0.0,
          40.0
      },
      {
          Form("ptw_compare_%s_cent_fwd", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
          "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
          "Centrality (%)",
          "1.6 < |y| < 2.4",
          0.0,
          90.0
      },
      {
          Form("ptw_compare_%s_cent_mid", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
          "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
          "Centrality (%)",
          "|y| < 1.6",
          0.0,
          90.0
      }
  };

  for (const auto& job : jobs) {
    const bool isCent = job.xTitle.Contains("Centrality");
    DrawSingleCompare(fPtW, fNoPtW, job, sample, isCent);
  }

  fPtW->Close();
  fNoPtW->Close();
  delete fPtW;
  delete fNoPtW;

  std::cout << "Saved plots under: " << kOutDir << std::endl;
}

void run_ptweight_compare_jpsi_pbpb_all() {
  plot_efficiency_ptweight_compare_jpsi_pbpb(1);
  plot_efficiency_ptweight_compare_jpsi_pbpb(2);
}

void plot_efficiency_ptweight_compare_jpsi_pbpb(const char* mode) {
  TString opt = mode ? mode : "";
  opt.ToLower();
  opt.ReplaceAll(" ", "");

  if (opt.IsNull() || opt == "all" || opt.Contains("run_ptweight_compare_jpsi_pbpb_all")) {
    run_ptweight_compare_jpsi_pbpb_all();
    return;
  }

  if (opt == "prompt" || opt == "1") {
    plot_efficiency_ptweight_compare_jpsi_pbpb(1);
    return;
  }

  if (opt == "nonprompt" || opt == "nprompt" || opt == "2") {
    plot_efficiency_ptweight_compare_jpsi_pbpb(2);
    return;
  }

  std::cerr << "[ERROR] Unknown mode: " << opt << std::endl;
  std::cerr << "Use one of: all, prompt, nonprompt, 1, 2" << std::endl;
}
