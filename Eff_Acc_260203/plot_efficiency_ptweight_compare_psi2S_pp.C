#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
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

const TString kOutDir = "plot_outputs_eff_ptw_compare_psi2S_pp";

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

void DrawSingleCompare(TFile* fPtW,
                       TFile* fNoPtW,
                       const PlotItem& cfg,
                       const TString& sampleLabel) {
  TH1* hPtW = LoadHistWithFallback(fPtW, cfg.histPtW, Form("h_ptw_%s", cfg.outTag.Data()));
  TH1* hNoPtW = LoadHistWithFallback(fNoPtW, cfg.histNoPtW, Form("h_noptw_%s", cfg.outTag.Data()));
  if (!hPtW || !hNoPtW) {
    delete hPtW;
    delete hNoPtW;
    return;
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

void DrawPtAndIntegratedPanel(TFile* fPtW,
                              TFile* fNoPtW,
                              const PlotItem& ptCfg,
                              const PlotItem& intCfg,
                              const TString& sampleLabel) {
  TH1* hPtW = LoadHistWithFallback(fPtW, ptCfg.histPtW, Form("h_ptw_%s", ptCfg.outTag.Data()));
  TH1* hNoPtW = LoadHistWithFallback(fNoPtW, ptCfg.histNoPtW, Form("h_noptw_%s", ptCfg.outTag.Data()));
  TH1* hIntPtW = LoadHistWithFallback(fPtW, intCfg.histPtW, Form("h_ptw_%s", intCfg.outTag.Data()));
  TH1* hIntNoPtW = LoadHistWithFallback(fNoPtW, intCfg.histNoPtW, Form("h_noptw_%s", intCfg.outTag.Data()));

  if (!hPtW || !hNoPtW || !hIntPtW || !hIntNoPtW) {
    delete hPtW;
    delete hNoPtW;
    delete hIntPtW;
    delete hIntNoPtW;
    return;
  }

  const TString rightPadLabel = intCfg.rapidityText.Contains(",")
                                    ? intCfg.rapidityText(0, intCfg.rapidityText.First(','))
                                    : intCfg.rapidityText;

  gStyle->SetOptStat(0);
  const double leftFrac = 800.0 / 940.0;
  TCanvas* c = new TCanvas(Form("c_pair_%s", ptCfg.outTag.Data()), "", 940, 700);
  TPad* padL = new TPad(Form("padL_%s", ptCfg.outTag.Data()), "", 0.00, 0.00, leftFrac, 1.00);
  TPad* padR = new TPad(Form("padR_%s", ptCfg.outTag.Data()), "", leftFrac, 0.00, 1.00, 1.00);
  padL->SetTicks(1, 1);
  padR->SetTicks(1, 1);
  padL->SetRightMargin(0.015);
  padR->SetLeftMargin(0.0);
  padL->Draw();
  padR->Draw();

  padL->cd();
  TPad* padLUp = new TPad(Form("padLUp_%s", ptCfg.outTag.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padLDn = new TPad(Form("padLDn_%s", ptCfg.outTag.Data()), "", 0.0, 0.0, 1.0, 0.28);
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
  hPtW->SetTitle("");
  hPtW->SetLineColor(kBlue + 1);
  hPtW->SetMarkerColor(kBlue + 1);
  hPtW->SetMarkerStyle(20);
  hPtW->SetLineWidth(2);
  hNoPtW->SetLineColor(kRed + 1);
  hNoPtW->SetMarkerColor(kRed + 1);
  hNoPtW->SetMarkerStyle(21);
  hNoPtW->SetLineWidth(2);
  hPtW->GetXaxis()->SetTitle(ptCfg.xTitle);
  hPtW->GetYaxis()->SetTitle("Efficiency");
  hPtW->GetXaxis()->SetRangeUser(ptCfg.xMin, ptCfg.xMax);
  hPtW->GetYaxis()->SetRangeUser(0.0, 1.2);
  hPtW->GetXaxis()->SetLabelSize(0.0);
  hPtW->GetXaxis()->SetTitleSize(0.0);
  hPtW->Draw("PE");
  hNoPtW->Draw("PE SAME");

  TLegend* legL = new TLegend(0.46, 0.72, 0.86, 0.88);
  legL->SetBorderSize(0);
  legL->SetFillStyle(0);
  legL->SetTextSize(0.030);
  legL->SetMargin(0.22);
  legL->AddEntry(hPtW, "With p_{T} weight", "lep");
  legL->AddEntry(hNoPtW, "No p_{T} weight", "lep");
  legL->Draw();
  drawText(sampleLabel.Data(), 0.20, 0.84, kBlack, 18);
  drawText(ptCfg.rapidityText.Data(), 0.20, 0.77, kBlack, 18);

  padLDn->cd();
  TH1* hRatioL = dynamic_cast<TH1*>(hNoPtW->Clone(Form("hRatioL_%s", ptCfg.outTag.Data())));
  hRatioL->SetDirectory(nullptr);
  hRatioL->Divide(hPtW);
  hRatioL->SetTitle("");
  hRatioL->GetXaxis()->SetTitle(ptCfg.xTitle);
  hRatioL->GetYaxis()->SetTitle("Ratio");
  hRatioL->GetYaxis()->CenterTitle();
  hRatioL->GetYaxis()->SetNdivisions(505);
  hRatioL->GetXaxis()->SetRangeUser(ptCfg.xMin, ptCfg.xMax);
  hRatioL->GetXaxis()->SetLabelSize(0.11);
  hRatioL->GetXaxis()->SetTitleSize(0.12);
  hRatioL->GetXaxis()->SetTitleOffset(1.0);
  hRatioL->GetYaxis()->SetLabelSize(0.10);
  hRatioL->GetYaxis()->SetTitleSize(0.10);
  hRatioL->GetYaxis()->SetTitleOffset(0.42);
  hRatioL->SetLineColor(kBlack);
  hRatioL->SetMarkerColor(kBlack);
  hRatioL->SetMarkerStyle(21);
  hRatioL->SetMarkerSize(1.0);

  padR->cd();
  TPad* padRUp = new TPad(Form("padRUp_%s", ptCfg.outTag.Data()), "", 0.0, 0.28, 1.0, 1.0);
  TPad* padRDn = new TPad(Form("padRDn_%s", ptCfg.outTag.Data()), "", 0.0, 0.0, 1.0, 0.28);
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
  TH1D* hFrame = new TH1D(Form("hFrame_%s", ptCfg.outTag.Data()), "", 1, 0.5, 1.5);
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

  TH1D* hIntA = new TH1D(Form("hIntA_%s", ptCfg.outTag.Data()), "", 1, 0.5, 1.5);
  TH1D* hIntB = new TH1D(Form("hIntB_%s", ptCfg.outTag.Data()), "", 1, 0.5, 1.5);
  hIntA->SetDirectory(nullptr);
  hIntB->SetDirectory(nullptr);
  hIntA->SetBinContent(1, hIntPtW->GetBinContent(1));
  hIntA->SetBinError(1, hIntPtW->GetBinError(1));
  hIntB->SetBinContent(1, hIntNoPtW->GetBinContent(1));
  hIntB->SetBinError(1, hIntNoPtW->GetBinError(1));
  hIntA->SetLineColor(kBlue + 1);
  hIntA->SetMarkerColor(kBlue + 1);
  hIntA->SetMarkerStyle(20);
  hIntA->SetLineWidth(2);
  hIntB->SetLineColor(kRed + 1);
  hIntB->SetMarkerColor(kRed + 1);
  hIntB->SetMarkerStyle(21);
  hIntB->SetLineWidth(2);
  hIntA->Draw("PE SAME");
  hIntB->Draw("PE SAME");

  padRDn->cd();
  TH1D* hRatioR = new TH1D(Form("hRatioR_%s", ptCfg.outTag.Data()), "", 1, 0.5, 1.5);
  hRatioR->SetDirectory(nullptr);
  const double ratioVal = (hIntA->GetBinContent(1) != 0.0) ? hIntB->GetBinContent(1) / hIntA->GetBinContent(1) : 0.0;
  hRatioR->SetBinContent(1, ratioVal);
  hRatioR->SetBinError(1, 0.0);
  hRatioR->SetTitle("");
  hRatioR->GetXaxis()->SetBinLabel(1, rightPadLabel);
  hRatioR->GetYaxis()->SetTitle("");
  hRatioR->GetYaxis()->CenterTitle();
  hRatioR->GetYaxis()->SetNdivisions(505);

  double yMinL = 1e9, yMaxL = -1e9;
  for (int i = 1; i <= hRatioL->GetNbinsX(); ++i) {
    const double v = hRatioL->GetBinContent(i);
    if (v <= 0.0) continue;
    if (v < yMinL) yMinL = v;
    if (v > yMaxL) yMaxL = v;
  }
  if (yMaxL < yMinL) {
    yMinL = 0.8;
    yMaxL = 1.2;
  }
  const double yMinR = (ratioVal > 0.0) ? ratioVal : 1.0;
  const double yMaxR = (ratioVal > 0.0) ? ratioVal : 1.0;
  double yMinCommon = (yMinL < yMinR) ? yMinL : yMinR;
  double yMaxCommon = (yMaxL > yMaxR) ? yMaxL : yMaxR;
  const double pad = 0.08;
  yMinCommon -= pad;
  yMaxCommon += pad;
  if (yMinCommon < 0.0) yMinCommon = 0.0;
  if (yMaxCommon <= yMinCommon) {
    yMinCommon = 0.8;
    yMaxCommon = 1.2;
  }
  hRatioL->GetYaxis()->SetRangeUser(yMinCommon, yMaxCommon);
  hRatioR->GetYaxis()->SetRangeUser(yMinCommon, yMaxCommon);

  padLDn->cd();
  hRatioL->Draw("PE");
  TLine* l1L = new TLine(ptCfg.xMin, 1.0, ptCfg.xMax, 1.0);
  l1L->SetLineStyle(2);
  l1L->SetLineColor(kGray + 2);
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
  hRatioR->SetMarkerStyle(21);
  hRatioR->SetMarkerSize(1.0);
  hRatioR->Draw("HIST");
  hRatioR->Draw("PE SAME");
  TLine* l1R = new TLine(0.5, 1.0, 1.5, 1.0);
  l1R->SetLineStyle(2);
  l1R->SetLineColor(kGray + 2);
  l1R->Draw("same");

  c->SaveAs(OutPath(ptCfg.outTag + "_withIntegrated", "png"));
  c->SaveAs(OutPath(ptCfg.outTag + "_withIntegrated", "pdf"));

  delete legL;
  delete hFrame;
  delete hRatioL;
  delete hRatioR;
  delete hIntA;
  delete hIntB;
  delete l1L;
  delete l1R;
  delete hPtW;
  delete hNoPtW;
  delete hIntPtW;
  delete hIntNoPtW;
  delete padLUp;
  delete padLDn;
  delete padRUp;
  delete padRDn;
  delete padL;
  delete padR;
  delete c;
}

}  // namespace

void plot_efficiency_ptweight_compare_psi2S_pp(int state = 1) {
  gStyle->SetOptStat(0);
  EnsureOutDir();

  const bool isPrompt = (state == 1);
  const TString sample = isPrompt ? "Prompt #psi(2S) pp" : "Nonprompt #psi(2S) pp";

  const TString filePtW = isPrompt
                              ? "./roots/mc_eff_vs_pt_rap_prompt_pp_psi2S_PtW1_tnp1_ctauCut_260323_2exp.root"
                              : "./roots/mc_eff_vs_pt_rap_nprompt_pp_psi2S_PtW1_tnp1_ctauCut_260323_2exp.root";

  const TString fileNoPtW = isPrompt
                                ? "./roots/mc_eff_vs_pt_rap_prompt_pp_psi2S_PtW1_tnp1_ctauCut_260323_noPtW.root"
                                : "./roots/mc_eff_vs_pt_rap_nprompt_pp_psi2S_PtW1_tnp1_ctauCut_260323_noPtW.root";

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
        Form("ptw_compare_%s_ptInt_fwd", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
        "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
        "Integrated p_{T} bin",
        "1.6 < |y| < 2.4, 3.5 < p_{T} < 40",
        0.0,
        2.0
      },
      {
        Form("ptw_compare_%s_ptInt_mid", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
        "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
        "Integrated p_{T} bin",
        "|y| < 1.6, 6.5 < p_{T} < 40",
        0.0,
        2.0
      },
      {
          Form("ptw_compare_%s_pt_mid", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
          "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
          "p_{T} (GeV/c)",
          "|y| < 1.6",
          0.0,
          40.0
      },
      {
          Form("ptw_compare_%s_pt_fwd", isPrompt ? "prompt" : "nprompt"),
          "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
          "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
          "p_{T} (GeV/c)",
          "1.6 < |y| < 2.4",
          0.0,
          40.0
      }
  };

  for (const auto& job : jobs) {
    DrawSingleCompare(fPtW, fNoPtW, job, sample);
  }

      const PlotItem ptMid = {
        Form("ptw_compare_%s_pt_mid", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
        "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6",
        "p_{T} (GeV/c)",
        "|y| < 1.6",
        0.0,
        40.0
      };
      const PlotItem intMid = {
        Form("ptw_compare_%s_ptInt_mid", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
        "mc_eff_Integrated_TnP1_PtW1_absy0_1p6",
        "Integrated p_{T} bin",
        "|y| < 1.6, 6.5 < p_{T} < 40",
        0.0,
        2.0
      };
      const PlotItem ptFwd = {
        Form("ptw_compare_%s_pt_fwd", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
        "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4",
        "p_{T} (GeV/c)",
        "1.6 < |y| < 2.4",
        0.0,
        40.0
      };
      const PlotItem intFwd = {
        Form("ptw_compare_%s_ptInt_fwd", isPrompt ? "prompt" : "nprompt"),
        "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
        "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4",
        "Integrated p_{T} bin",
        "1.6 < |y| < 2.4, 3.5 < p_{T} < 40",
        0.0,
        2.0
      };

      DrawPtAndIntegratedPanel(fPtW, fNoPtW, ptMid, intMid, sample);
      DrawPtAndIntegratedPanel(fPtW, fNoPtW, ptFwd, intFwd, sample);

  fPtW->Close();
  fNoPtW->Close();
  delete fPtW;
  delete fNoPtW;

  std::cout << "Saved plots under: " << kOutDir << std::endl;
}

void run_ptweight_compare_psi2S_pp_all() {
  plot_efficiency_ptweight_compare_psi2S_pp(1);
  plot_efficiency_ptweight_compare_psi2S_pp(2);
}

void plot_efficiency_ptweight_compare_psi2S_pp(const char* mode) {
  TString opt = mode ? mode : "";
  opt.ToLower();
  opt.ReplaceAll(" ", "");

  if (opt.IsNull() || opt == "all" || opt.Contains("run_ptweight_compare_psi2s_pp_all")) {
    run_ptweight_compare_psi2S_pp_all();
    return;
  }

  if (opt == "prompt" || opt == "1") {
    plot_efficiency_ptweight_compare_psi2S_pp(1);
    return;
  }

  if (opt == "nonprompt" || opt == "nprompt" || opt == "2") {
    plot_efficiency_ptweight_compare_psi2S_pp(2);
    return;
  }

  std::cerr << "[ERROR] Unknown mode: " << opt << std::endl;
  std::cerr << "Use one of: all, prompt, nonprompt, 1, 2" << std::endl;
}
