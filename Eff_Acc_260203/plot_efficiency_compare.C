#include <iostream>
#include <vector>
#include <cmath>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_v2mass.C"


struct EffInput {
  TString filePath;
  TString histPath;
  TString label;
  Color_t color;
  Style_t marker;
};

namespace {

const TString kPlotOutputDir = "plot_outputs_eff_compare_jpsi";

void EnsureOutputDir() {
  if (gSystem) gSystem->mkdir(kPlotOutputDir, true);
}

TString BuildOutputPath(const TString& baseName, const TString& ext) {
  return Form("%s/%s.%s", kPlotOutputDir.Data(), baseName.Data(), ext.Data());
}

TH1* LoadHistogram(const EffInput& in, TFile*& fileHandle) {
  fileHandle = TFile::Open(in.filePath, "READ");
  if (!fileHandle || fileHandle->IsZombie()) {
    std::cerr << "[ERROR] Cannot open file: " << in.filePath << std::endl;
    return nullptr;
  }

  TH1* h = dynamic_cast<TH1*>(fileHandle->Get(in.histPath));
  if (!h) {
    std::cerr << "[ERROR] Cannot find histogram: " << in.histPath
              << " in file " << in.filePath << std::endl;
    return nullptr;
  }

  TH1* hClone = dynamic_cast<TH1*>(h->Clone(Form("h_%s", in.label.Data())));
  if (!hClone) {
    std::cerr << "[ERROR] Failed to clone histogram: " << in.histPath << std::endl;
    return nullptr;
  }
  hClone->SetDirectory(nullptr);
  return hClone;
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

void StyleHistogram(TH1* h, const EffInput& in, const TString& xTitle, const TString& yTitle) {
  h->SetLineColor(in.color);
  h->SetMarkerColor(in.color);
  h->SetMarkerStyle(in.marker);
  h->SetMarkerSize(1.1);
  h->SetLineWidth(2);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
}

void GetAutoRatioYRange(const TH1* h,
                        double& yMin,
                        double& yMax,
                        double fallbackMin = 0.8,
                        double fallbackMax = 1.2) {
  if (!h) {
    yMin = fallbackMin;
    yMax = fallbackMax;
    return;
  }

  bool hasValid = false;
  double minVal = 1e9;
  double maxVal = -1e9;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double v = h->GetBinContent(i);
    const double e = h->GetBinError(i);
    if (!std::isfinite(v) || !std::isfinite(e)) continue;
    if (v <= 0.0 && e <= 0.0) continue;

    const double lo = v - e;
    const double hi = v + e;
    if (!std::isfinite(lo) || !std::isfinite(hi)) continue;
    if (hi <= 0.0) continue;

    hasValid = true;
    if (lo < minVal) minVal = lo;
    if (hi > maxVal) maxVal = hi;
  }

  if (!hasValid) {
    yMin = fallbackMin;
    yMax = fallbackMax;
    return;
  }

  if (minVal > 1.0) minVal = 1.0;
  if (maxVal < 1.0) maxVal = 1.0;

  double span = maxVal - minVal;
  if (span < 0.08) span = 0.08;
  const double pad = 0.25 * span;

  double yMin = minVal - pad;
  double yMax = maxVal + pad;
  if (yMin < 0.0) yMin = 0.0;
  if (yMax <= yMin) {
    yMin = fallbackMin;
    yMax = fallbackMax;
  }
}

void SetAutoRatioYAxis(TH1* h, double fallbackMin = 0.8, double fallbackMax = 1.2) {
  if (!h) return;
  double yMin = fallbackMin;
  double yMax = fallbackMax;
  GetAutoRatioYRange(h, yMin, yMax, fallbackMin, fallbackMax);
  h->GetYaxis()->SetRangeUser(yMin, yMax);
}

}  // namespace

void draw_efficiency_compare(const std::vector<EffInput>& inputs,
                             const TString& outName = "plot_efficiency_compare",
                             const TString& xTitle = "p_{T} (GeV/c)",
                             const TString& yTitle = "Efficiency",
                             double xMin = -1,
                             double xMax = -1,
                             double yMin = 0.0,
                             double yMax = 1.05,
                             const TString& drawOpt = "PE",
                             const TString& rapidityText = "") {
  if (inputs.size() < 2) {
    std::cerr << "[ERROR] Need at least 2 histograms to compare." << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  EnsureOutputDir();

  std::vector<TFile*> files;
  std::vector<TH1*> hists;
  files.reserve(inputs.size());
  hists.reserve(inputs.size());

  for (const auto& in : inputs) {
    TFile* f = nullptr;
    TH1* h = LoadHistogram(in, f);
    files.push_back(f);
    hists.push_back(h);
  }

  for (TH1* h : hists) {
    if (!h) {
      std::cerr << "[ERROR] Aborting draw due to invalid histogram input." << std::endl;
      for (TFile* f : files) {
        if (f) {
          f->Close();
          delete f;
        }
      }
      for (TH1* hh : hists) delete hh;
      return;
    }
  }

  if (xTitle.Contains("Centrality")) {
    for (size_t i = 0; i < hists.size(); ++i) {
      TH1* hScaled = RescaleXAxis(hists[i], 0.5, Form("centScaled_%zu", i));
      if (!hScaled) continue;
      delete hists[i];
      hists[i] = hScaled;
    }
  }

  TCanvas* c = new TCanvas(Form("c_%s", outName.Data()), "", 800, 700);
  TPad* padUp = new TPad(Form("padUp_%s", outName.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padDn = new TPad(Form("padDn_%s", outName.Data()), "", 0.0, 0.0, 1.0, 0.28);
  padUp->SetBottomMargin(0.03);
  padUp->SetTicks(1, 1);
  padDn->SetTopMargin(0.02);
  padDn->SetBottomMargin(0.33);
  padDn->SetTicks(1, 1);
  padUp->Draw();
  padDn->Draw();

  padUp->cd();

  for (size_t i = 0; i < hists.size(); ++i) {
    StyleHistogram(hists[i], inputs[i], xTitle, yTitle);
    if (xMin < xMax) hists[i]->GetXaxis()->SetRangeUser(xMin, xMax);
    hists[i]->GetYaxis()->SetRangeUser(yMin, yMax);
    hists[i]->GetXaxis()->SetLabelSize(0.0);
    hists[i]->GetXaxis()->SetTitleSize(0.0);

    if (i == 0) {
      hists[i]->Draw(drawOpt);
    } else {
      hists[i]->Draw(drawOpt + " SAME");
    }
  }

  TLegend* leg = new TLegend(0.44, 0.68, 0.84, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->SetMargin(0.22);
  for (size_t i = 0; i < hists.size(); ++i) {
    leg->AddEntry(hists[i], inputs[i].label, "lep");
  }
  leg->Draw();

  if (!rapidityText.IsNull()) {
    const float posX = 0.20;
    const float posY = 0.84;
    const int textColor = kBlack;
    const int textSize = 18;
    drawText(rapidityText.Data(), posX, posY, textColor, textSize);

    if (xTitle.Contains("p_{T}")) {
      drawText("Cent. 0-90%", posX, posY - 0.07, textColor, textSize);
    } else if (xTitle.Contains("Centrality")) {
      TString ptText = "";
      if (rapidityText.Contains("1.6 < |y| < 2.4")) {
        ptText = "3.5 < p_{T} < 40 GeV/c";
      } else if (rapidityText.Contains("|y| < 1.6")) {
        ptText = "6.5 < p_{T} < 40 GeV/c";
      }
      if (!ptText.IsNull()) {
        drawText(ptText.Data(), posX, posY - 0.07, textColor, textSize);
      }
    }
  }

  if (hists.size() >= 2) {
    padDn->cd();
    TH1* hRatio = dynamic_cast<TH1*>(hists[1]->Clone(Form("hRatio_%s", outName.Data())));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hists[0]);
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(xTitle);
    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetNdivisions(505);
    SetAutoRatioYAxis(hRatio);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.42);
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(hists[1]->GetMarkerStyle());
    hRatio->SetMarkerSize(1.0);
    hRatio->Draw("PE");

    TLine* l1 = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    l1->SetLineStyle(2);
    l1->SetLineColor(kGray + 2);
    l1->Draw("same");
  }

  c->SaveAs(BuildOutputPath(outName, "png"));
  c->SaveAs(BuildOutputPath(outName, "pdf"));

  for (TH1* h : hists) delete h;
  for (TFile* f : files) {
    if (f) {
      f->Close();
      delete f;
    }
  }
}

// Example 1) PbPb J/#psi, forward rapidity, compare before/after ctau cut
void run_compare_jpsi_pbpb_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
      "Prompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
      "Prompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_prompt_fwd",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 2) PbPb J/#psi, mid rapidity, compare before/after ctau cut
void run_compare_jpsi_pbpb_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
      "Prompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
      "Prompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_prompt_mid",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

// Example 3) PbPb J/#psi, forward rapidity, efficiency vs centrality
void run_compare_jpsi_pbpb_cent_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
      "Prompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
      "Prompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_prompt_cent_fwd",
                          "Centrality (%)",
                          "Efficiency",
                          0.0,
                          90.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 4) PbPb J/#psi, mid rapidity, efficiency vs centrality
void run_compare_jpsi_pbpb_cent_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
      "Prompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
      "Prompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_prompt_cent_mid",
                          "Centrality (%)",
                          "Efficiency",
                          0.0,
                          90.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

// Example 4-1) PbPb nonprompt J/#psi, forward rapidity, efficiency vs pT
void run_compare_jpsi_pbpb_nonprompt_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
        "Nonprompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4",
        "Nonprompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_nonprompt_fwd",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 4-2) PbPb nonprompt J/#psi, mid rapidity, efficiency vs pT
void run_compare_jpsi_pbpb_nonprompt_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
        "Nonprompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6",
        "Nonprompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_nonprompt_mid",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

// Example 4-3) PbPb nonprompt J/#psi, forward rapidity, efficiency vs centrality
void run_compare_jpsi_pbpb_nonprompt_cent_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
        "Nonprompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4",
        "Nonprompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_nonprompt_cent_fwd",
                          "Centrality (%)",
                          "Efficiency",
                          0.0,
                          90.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 4-4) PbPb nonprompt J/#psi, mid rapidity, efficiency vs centrality
void run_compare_jpsi_pbpb_nonprompt_cent_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
        "Nonprompt J/#psi PbPb (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6",
        "Nonprompt J/#psi PbPb (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pbpb_nonprompt_cent_mid",
                          "Centrality (%)",
                          "Efficiency",
                          0.0,
                          90.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

// Example 5) pp J/#psi, forward rapidity, efficiency vs pT
void run_compare_jpsi_pp_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
       "Prompt J/#psi pp (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
       "Prompt J/#psi pp (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pp_prompt_fwd",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 6) pp J/#psi, mid rapidity, efficiency vs pT
void run_compare_jpsi_pp_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
       "Prompt J/#psi pp (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
       "Prompt J/#psi pp (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pp_prompt_mid",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

// Example 6-1) pp nonprompt J/#psi, forward rapidity, efficiency vs pT
void run_compare_jpsi_pp_nonprompt_fwd() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
        "Nonprompt J/#psi pp (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
        "Nonprompt J/#psi pp (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pp_nonprompt_fwd",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "1.6 < |y| < 2.4");
}

// Example 6-2) pp nonprompt J/#psi, mid rapidity, efficiency vs pT
void run_compare_jpsi_pp_nonprompt_mid() {
  std::vector<EffInput> inputs = {
      {"./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
        "Nonprompt J/#psi pp (260219)", kBlue + 1, 20},
      {"./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
       "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
        "Nonprompt J/#psi pp (260310)", kRed + 1, 21},
  };

  draw_efficiency_compare(inputs,
                          "plot_eff_compare_jpsi_pp_nonprompt_mid",
                          "p_{T} (GeV/c)",
                          "Efficiency",
                          0.0,
                          40.0,
                          0.0,
                          1.2,
                          "PE",
                          "|y| < 1.6");
}

void draw_pp_with_single_integrated_rightpad(const TString& ptHistName,
                                             const TString& intHistName,
                                             const TString& rightPadLabel,
                                             double xMin,
                                             const TString& outName,
                                             const TString& fileA,
                                             const TString& fileB,
                                             const TString& leftLabelA,
                                             const TString& leftLabelB,
                                             Color_t c1,
                                             Color_t c2,
                                             Style_t m1,
                                             Style_t m2) {

  TFile* fA = TFile::Open(fileA, "READ");
  TFile* fB = TFile::Open(fileB, "READ");
  if (!fA || fA->IsZombie() || !fB || fB->IsZombie()) {
    std::cerr << "[ERROR] Cannot open pp input files for integrated-right-pad plot." << std::endl;
    if (fA) {
      fA->Close();
      delete fA;
    }
    if (fB) {
      fB->Close();
      delete fB;
    }
    return;
  }

  TH1* hPtA_src = dynamic_cast<TH1*>(fA->Get(ptHistName));
  TH1* hPtB_src = dynamic_cast<TH1*>(fB->Get(ptHistName));
  TH1* hIntA_src = dynamic_cast<TH1*>(fA->Get(intHistName));
  TH1* hIntB_src = dynamic_cast<TH1*>(fB->Get(intHistName));

  if (!hPtA_src || !hPtB_src || !hIntA_src || !hIntB_src) {
    std::cerr << "[ERROR] Missing one or more pp histograms needed for integrated-right-pad plot." << std::endl;
    fA->Close();
    fB->Close();
    delete fA;
    delete fB;
    return;
  }

  TH1* hPtA = dynamic_cast<TH1*>(hPtA_src->Clone(Form("hPtA_%s", outName.Data())));
  TH1* hPtB = dynamic_cast<TH1*>(hPtB_src->Clone(Form("hPtB_%s", outName.Data())));
  hPtA->SetDirectory(nullptr);
  hPtB->SetDirectory(nullptr);

  TH1D* hIntA = new TH1D(Form("hIntA_%s", outName.Data()), "", 1, 0.5, 1.5);
  TH1D* hIntB = new TH1D(Form("hIntB_%s", outName.Data()), "", 1, 0.5, 1.5);
  hIntA->SetDirectory(nullptr);
  hIntB->SetDirectory(nullptr);
  hIntA->SetBinContent(1, hIntA_src->GetBinContent(1));
  hIntA->SetBinError(1, hIntA_src->GetBinError(1));
  hIntB->SetBinContent(1, hIntB_src->GetBinContent(1));
  hIntB->SetBinError(1, hIntB_src->GetBinError(1));

  gStyle->SetOptStat(0);
  EnsureOutputDir();
  const double leftFrac = 800.0 / 940.0;  // Keep left panel same pixel width as standard 800px plots.
  TCanvas* c = new TCanvas(Form("c_%s", outName.Data()), "", 940, 700);
  TPad* padL = new TPad(Form("padL_%s", outName.Data()), "", 0.00, 0.00, leftFrac, 1.00);
  TPad* padR = new TPad(Form("padR_%s", outName.Data()), "", leftFrac, 0.00, 1.00, 1.00);
  padL->SetTicks(1, 1);
  padR->SetTicks(1, 1);
  padL->SetRightMargin(0.015);
  padR->SetLeftMargin(0.0);
  padL->Draw();
  padR->Draw();

  padL->cd();
  TPad* padLUp = new TPad(Form("padLUp_%s", outName.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padLDn = new TPad(Form("padLDn_%s", outName.Data()), "", 0.0, 0.0, 1.0, 0.28);
  padLUp->SetBottomMargin(0.03);
  padLUp->SetRightMargin(0.015);
  padLUp->SetTicks(1, 1);
  padLDn->SetTopMargin(0.02);
  padLDn->SetBottomMargin(0.33);
  padLDn->SetRightMargin(0.015);
  padLDn->SetTicks(1, 1);
  padLUp->Draw();
  padLDn->Draw();

  padLUp->cd();
  hPtA->SetTitle("");
  hPtA->SetLineColor(c1);
  hPtA->SetMarkerColor(c1);
  hPtA->SetMarkerStyle(m1);
  hPtA->SetLineWidth(2);
  hPtB->SetLineColor(c2);
  hPtB->SetMarkerColor(c2);
  hPtB->SetMarkerStyle(m2);
  hPtB->SetLineWidth(2);
  hPtA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtA->GetYaxis()->SetTitle("Efficiency");
  hPtA->GetXaxis()->SetRangeUser(xMin, 40.0);
  hPtA->GetYaxis()->SetRangeUser(0.0, 1.2);
  hPtA->GetXaxis()->SetLabelSize(0.0);
  hPtA->GetXaxis()->SetTitleSize(0.0);
  hPtA->Draw("PE");
  hPtB->Draw("PE SAME");

  TLegend* legL = new TLegend(0.46, 0.72, 0.86, 0.88);
  legL->SetBorderSize(0);
  legL->SetFillStyle(0);
  legL->SetTextSize(0.030);
  legL->SetMargin(0.22);
  legL->AddEntry(hPtA, leftLabelA, "lep");
  legL->AddEntry(hPtB, leftLabelB, "lep");
  legL->Draw();
  drawText("Cent. 0-90%", 0.20, 0.84, kBlack, 18);

  padLDn->cd();
  TH1* hRatioL = dynamic_cast<TH1*>(hPtB->Clone(Form("hRatioL_%s", outName.Data())));
  hRatioL->SetDirectory(nullptr);
  hRatioL->Divide(hPtA);
  hRatioL->SetTitle("");
  hRatioL->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatioL->GetYaxis()->SetTitle("Ratio");
  hRatioL->GetYaxis()->CenterTitle();
  hRatioL->GetYaxis()->SetNdivisions(505);
  hRatioL->GetXaxis()->SetRangeUser(xMin, 40.0);
  hRatioL->GetXaxis()->SetLabelSize(0.11);
  hRatioL->GetXaxis()->SetTitleSize(0.12);
  hRatioL->GetXaxis()->SetTitleOffset(1.0);
  hRatioL->GetYaxis()->SetLabelSize(0.10);
  hRatioL->GetYaxis()->SetTitleSize(0.10);
  hRatioL->GetYaxis()->SetTitleOffset(0.42);
  hRatioL->SetLineColor(kBlack);
  hRatioL->SetMarkerColor(kBlack);
  hRatioL->SetMarkerStyle(m2);
  hRatioL->SetMarkerSize(1.0);
  hRatioL->Draw("PE");

  TLine* l1L = new TLine(xMin, 1.0, 40.0, 1.0);
  l1L->SetLineStyle(2);
  l1L->SetLineColor(kGray + 2);
  l1L->Draw("same");

  padR->cd();
  TPad* padRUp = new TPad(Form("padRUp_%s", outName.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padRDn = new TPad(Form("padRDn_%s", outName.Data()), "", 0.0, 0.0, 1.0, 0.28);
  padRUp->SetBottomMargin(0.03);
  padRUp->SetLeftMargin(0.015);
  padRUp->SetTicks(1, 1);
  padRDn->SetTopMargin(0.02);
  padRDn->SetBottomMargin(0.33);
  padRDn->SetLeftMargin(0.015);
  padRDn->SetTicks(1, 1);
  padRUp->Draw();
  padRDn->Draw();

  padRUp->cd();
  TH1D* hFrame = new TH1D(Form("hFrame_%s", outName.Data()), "", 1, 0.5, 1.5);
  hFrame->SetDirectory(nullptr);
  hFrame->SetTitle("");
  hFrame->GetYaxis()->SetTitle("");
  hFrame->GetXaxis()->SetBinLabel(1, rightPadLabel);
  hFrame->GetYaxis()->SetRangeUser(0.0, 1.2);
  hFrame->GetXaxis()->SetLabelSize(0.0);
  hFrame->GetXaxis()->SetTitleSize(0.0);
  hFrame->GetYaxis()->SetLabelSize(0.0);
  hFrame->GetYaxis()->SetTitleSize(0.0);
  hFrame->GetYaxis()->SetTickLength(0.080);
  hFrame->Draw("AXIS");

  hIntA->SetLineColor(c1);
  hIntA->SetMarkerColor(c1);
  hIntA->SetMarkerStyle(m1);
  hIntA->SetLineWidth(2);
  hIntB->SetLineColor(c2);
  hIntB->SetMarkerColor(c2);
  hIntB->SetMarkerStyle(m2);
  hIntB->SetLineWidth(2);
  hIntA->Draw("PE SAME");
  hIntB->Draw("PE SAME");

  padRDn->cd();
  TH1D* hRatioR = new TH1D(Form("hRatioR_%s", outName.Data()), "", 1, 0.5, 1.5);
  hRatioR->SetDirectory(nullptr);
  hRatioR->SetBinContent(1, (hIntA->GetBinContent(1) != 0.0) ? hIntB->GetBinContent(1) / hIntA->GetBinContent(1) : 0.0);
  hRatioR->SetBinError(1, 0.0);
  hRatioR->SetTitle("");
  hRatioR->GetXaxis()->SetBinLabel(1, rightPadLabel);
  hRatioR->GetYaxis()->SetTitle("");
  hRatioR->GetYaxis()->CenterTitle();
  hRatioR->GetYaxis()->SetNdivisions(505);
  double yMinL = 0.8, yMaxL = 1.2;
  double yMinR = 0.8, yMaxR = 1.2;
  GetAutoRatioYRange(hRatioL, yMinL, yMaxL);
  GetAutoRatioYRange(hRatioR, yMinR, yMaxR);
  const double yMinCommon = (yMinL < yMinR) ? yMinL : yMinR;
  const double yMaxCommon = (yMaxL > yMaxR) ? yMaxL : yMaxR;
  hRatioL->GetYaxis()->SetRangeUser(yMinCommon, yMaxCommon);
  hRatioR->GetYaxis()->SetRangeUser(yMinCommon, yMaxCommon);

  padLDn->cd();
  hRatioL->Draw("PE");
  l1L->Draw("same");

  padRDn->cd();
  hRatioR->GetXaxis()->SetLabelSize(0.20);
  hRatioR->GetXaxis()->SetTitleSize(0.16);
  hRatioR->GetYaxis()->SetLabelSize(0.0);
  hRatioR->GetYaxis()->SetTitleSize(0.0);
  hRatioR->GetXaxis()->SetTickLength(0.16);
  hRatioR->GetYaxis()->SetTickLength(0.080);
  hRatioR->SetLineColor(kBlack);
  hRatioR->SetLineWidth(2);
  hRatioR->SetMarkerColor(kBlack);
  hRatioR->SetMarkerStyle(m2);
  hRatioR->SetMarkerSize(1.0);
  hRatioR->Draw("HIST");
  hRatioR->Draw("PE SAME");

  TLine* l1R = new TLine(0.5, 1.0, 1.5, 1.0);
  l1R->SetLineStyle(2);
  l1R->SetLineColor(kGray + 2);
  l1R->Draw("same");

  c->SaveAs(BuildOutputPath(outName, "png"));
  c->SaveAs(BuildOutputPath(outName, "pdf"));

  delete hFrame;
  delete hRatioR;
  delete hIntA;
  delete hIntB;
  delete hPtA;
  delete hPtB;
  fA->Close();
  fB->Close();
  delete fA;
  delete fB;
}

// Example 7) pp J/#psi forward pT + right pad with integrated forward bin only
void run_compare_jpsi_pp_fwd_with_integrated_rightpad() {
  draw_pp_with_single_integrated_rightpad("mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
                                          "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
                                          "1.6 < |y| < 2.4",
                                          0.0,
                                          "plot_eff_compare_jpsi_pp_prompt_fwd_with_integrated_rightpad",
                                          "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
                                          "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
                                          "Prompt J/#psi pp (260219)",
                                          "Prompt J/#psi pp (260310)",
                                          kBlue + 1,
                                          kRed + 1,
                                          20,
                                          21);
}

// Example 8) pp J/#psi mid pT + right pad with integrated mid bin only
void run_compare_jpsi_pp_mid_with_integrated_rightpad() {
  draw_pp_with_single_integrated_rightpad("mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
                                          "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
                                          "|y| < 1.6",
                                          0.0,
                                          "plot_eff_compare_jpsi_pp_prompt_mid_with_integrated_rightpad",
                                          "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
                                          "./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
                                          "Prompt J/#psi pp (260219)",
                                          "Prompt J/#psi pp (260310)",
                                          kBlue + 1,
                                          kRed + 1,
                                          20,
                                          21);
}

// Example 9) pp nonprompt forward pT + right pad with integrated forward bin only
void run_compare_jpsi_pp_nonprompt_fwd_with_integrated_rightpad() {
  draw_pp_with_single_integrated_rightpad("mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
                                          "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
                                          "1.6 < |y| < 2.4",
                                          0.0,
                                          "plot_eff_compare_jpsi_pp_nonprompt_fwd_with_integrated_rightpad",
                                          "./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
                                          "./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
                                          "Nonprompt J/#psi pp (260219)",
                                          "Nonprompt J/#psi pp (260310)",
                                          kBlue + 1,
                                          kRed + 1,
                                          20,
                                          21);
}

// Example 10) pp nonprompt mid pT + right pad with integrated mid bin only
void run_compare_jpsi_pp_nonprompt_mid_with_integrated_rightpad() {
  draw_pp_with_single_integrated_rightpad("mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
                                          "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
                                          "|y| < 1.6",
                                          0.0,
                                          "plot_eff_compare_jpsi_pp_nonprompt_mid_with_integrated_rightpad",
                                          "./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260219.root",
                                          "./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260310.root",
                                          "Nonprompt J/#psi pp (260219)",
                                          "Nonprompt J/#psi pp (260310)",
                                          kBlue + 1,
                                          kRed + 1,
                                          20,
                                          21);
}

// Quick entry point
void run_compare_all_default() {
  run_compare_jpsi_pbpb_fwd();
  run_compare_jpsi_pbpb_mid();
  run_compare_jpsi_pbpb_cent_fwd();
  run_compare_jpsi_pbpb_cent_mid();
  run_compare_jpsi_pbpb_nonprompt_fwd();
  run_compare_jpsi_pbpb_nonprompt_mid();
  run_compare_jpsi_pbpb_nonprompt_cent_fwd();
  run_compare_jpsi_pbpb_nonprompt_cent_mid();
  run_compare_jpsi_pp_fwd();
  run_compare_jpsi_pp_mid();
  run_compare_jpsi_pp_nonprompt_fwd();
  run_compare_jpsi_pp_nonprompt_mid();
  run_compare_jpsi_pp_fwd_with_integrated_rightpad();
  run_compare_jpsi_pp_mid_with_integrated_rightpad();
  run_compare_jpsi_pp_nonprompt_fwd_with_integrated_rightpad();
  run_compare_jpsi_pp_nonprompt_mid_with_integrated_rightpad();
}

// Entry point for `rq plot_efficiency_compare.C`
void plot_efficiency_compare() {
  run_compare_all_default();
}
