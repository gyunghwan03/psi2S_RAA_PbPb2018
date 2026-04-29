#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
const TString kRefDir = "/data/users/pjgwak/work/daily_code_tracker/2026/04/11_run2_eff_study/acc_eff_skim/skim_roots";
const char *kRefLegendLabel = "Gwak";
const char *kNewLegendLabel = "Bak";

TString BaseDir()
{
  return gSystem->DirName(__FILE__);
}

TString OutDir(const char *subdir)
{
  return Form("%s/plot_outputs_acc_eff_compare_20260411/%s", BaseDir().Data(), subdir);
}

TString RefAccFile(bool isPbPb, bool isPrompt)
{
  return Form("%s/acc_%s_ppInput_isMC1_%s_ncollW0_genW1_ptW1.root",
              kRefDir.Data(), isPbPb ? "PbPb2018" : "pp2018", isPrompt ? "PR" : "NP");
}

TString NewAccFile(bool isPbPb, bool isPrompt)
{
  return Form("%s/roots/acceptance_%s_GenOnly_wgt1_%s_SysUp0_260429_1d.root",
              BaseDir().Data(), isPrompt ? "PromptJpsi" : "BtoJpsi", isPbPb ? "PbPb" : "pp");
}

TString RefEffFile(bool isPbPb, bool isPrompt)
{
  if (isPbPb)
    return Form("%s/eff_PbPb2018_isMC1_%s_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
                kRefDir.Data(), isPrompt ? "PR" : "NP");
  return Form("%s/eff_pp5p02TeV_isMC1_%s_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
              kRefDir.Data(), isPrompt ? "PR" : "NP");
}

TString NewEffFile(bool isPbPb, bool isPrompt)
{
  if (isPbPb)
    return Form("%s/roots/mc_eff_vs_pt_cent_0_to_180_rap_%s_pbpb_JPsi_PtWnomi_tnp1_260429_1d.root",
                BaseDir().Data(), isPrompt ? "prompt" : "nprompt");
  return Form("%s/roots/mc_eff_vs_pt_rap_%s_pp_Jpsi_PtWnomi_tnp1_260429_1d.root",
              BaseDir().Data(), isPrompt ? "prompt" : "nprompt");
}

TString CollTag(bool isPbPb)
{
  return isPbPb ? "pbpb" : "pp";
}

TString StateTag(bool isPrompt)
{
  return isPrompt ? "prompt" : "nonprompt";
}

TString StateLabel(bool isPrompt)
{
  return isPrompt ? "Prompt J/#psi" : "Nonprompt J/#psi";
}

TH1 *LoadHist(TFile *file, const char *key, const char *name, bool required = true)
{
  if (!file)
    return nullptr;
  TH1 *hist = dynamic_cast<TH1 *>(file->Get(key));
  if (!hist)
  {
    if (required)
      std::cout << "[ERROR] missing histogram: " << key << "\n";
    else
      std::cout << "[WARN] missing histogram: " << key << "\n";
    return nullptr;
  }
  TH1 *clone = static_cast<TH1 *>(hist->Clone(name));
  clone->SetDirectory(nullptr);
  return clone;
}

bool SameBinning(const TH1 *a, const TH1 *b)
{
  if (!a || !b || a->GetNbinsX() != b->GetNbinsX())
    return false;
  for (int i = 1; i <= a->GetNbinsX() + 1; ++i)
  {
    if (std::fabs(a->GetXaxis()->GetBinLowEdge(i) - b->GetXaxis()->GetBinLowEdge(i)) > 1e-6)
      return false;
  }
  return true;
}

void Style(TH1 *hist, int color, int marker)
{
  if (!hist)
    return;
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
}

TLegend *BuildLegend(TH1 *ref, TH1 *cur)
{
  TLegend *leg = new TLegend(0.62, 0.67, 0.91, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
  leg->SetEntrySeparation(0.25);
  leg->SetMargin(0.28);
  leg->AddEntry(ref, kRefLegendLabel, "lep");
  leg->AddEntry(cur, kNewLegendLabel, "lep");
  return leg;
}

TH1D *BuildIntegratedRatioHist(const TH1 *num, const TH1 *den, const char *name)
{
  if (!num || !den)
    return nullptr;

  double sumNum = 0.0;
  double err2Num = 0.0;
  double sumDen = 0.0;
  double err2Den = 0.0;
  for (int ib = 1; ib <= num->GetNbinsX(); ++ib)
  {
    sumNum += num->GetBinContent(ib);
    err2Num += std::pow(num->GetBinError(ib), 2);
  }
  for (int ib = 1; ib <= den->GetNbinsX(); ++ib)
  {
    sumDen += den->GetBinContent(ib);
    err2Den += std::pow(den->GetBinError(ib), 2);
  }

  TH1D *out = new TH1D(name, "", 1, 0.5, 1.5);
  out->SetDirectory(nullptr);
  if (sumDen <= 0.0)
    return out;

  const double value = sumNum / sumDen;
  double error = 0.0;
  if (sumNum > 0.0)
  {
    const double relNum = std::sqrt(err2Num) / sumNum;
    const double relDen = std::sqrt(err2Den) / sumDen;
    error = value * std::sqrt(relNum * relNum + relDen * relDen);
  }
  out->SetBinContent(1, value);
  out->SetBinError(1, error);
  return out;
}

TH1D *CloneAsIntegrated(const TH1 *hist, const char *name)
{
  if (!hist)
    return nullptr;
  TH1D *out = new TH1D(name, "", 1, 0.5, 1.5);
  out->SetDirectory(nullptr);
  out->SetBinContent(1, hist->GetBinContent(1));
  out->SetBinError(1, hist->GetBinError(1));
  return out;
}

void DrawRatioPanel(TH1 *ref, TH1 *cur, const char *xTitle)
{
  if (!ref || !cur)
    return;

  if (SameBinning(ref, cur))
  {
    TH1 *ratio = static_cast<TH1 *>(cur->Clone(Form("%s_ratio", cur->GetName())));
    ratio->SetDirectory(nullptr);
    ratio->Divide(ref);
    ratio->SetTitle("");
    ratio->GetXaxis()->SetTitle(xTitle);
    ratio->GetYaxis()->SetTitle("Bak/Gwak");
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetLabelSize(0.09);
    ratio->GetYaxis()->SetTitleOffset(0.55);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetRangeUser(0.70, 1.30);
    Style(ratio, kBlack, 20);
    ratio->Draw("PE");
    TLine *one = new TLine(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2);
    one->Draw("SAME");
    return;
  }

  TH1 *frame = static_cast<TH1 *>(ref->Clone(Form("%s_ratio_frame", ref->GetName())));
  frame->SetDirectory(nullptr);
  frame->Reset("ICES");
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetYaxis()->SetTitle("Bak/Gwak");
  frame->GetXaxis()->SetTitleSize(0.12);
  frame->GetXaxis()->SetLabelSize(0.10);
  frame->GetYaxis()->SetTitleSize(0.10);
  frame->GetYaxis()->SetLabelSize(0.09);
  frame->GetYaxis()->SetTitleOffset(0.55);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetRangeUser(0.70, 1.30);
  frame->Draw();
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.08);
  text.DrawLatex(0.20, 0.42, "bin mismatch");
}

void DrawBasicCompare(const char *canvasName,
                      const char *label,
                      const char *yTitle,
                      const char *xTitle,
                      TH1 *ref,
                      TH1 *cur,
                      const TString &outDir)
{
  if (!ref || !cur)
    return;

  TCanvas *canvas = new TCanvas(canvasName, canvasName, 820, 820);
  TPad *top = new TPad(Form("%s_top", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *bot = new TPad(Form("%s_bot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  top->SetBottomMargin(0.03);
  top->SetLeftMargin(0.13);
  top->SetRightMargin(0.04);
  top->SetTicks(1, 1);
  bot->SetTopMargin(0.03);
  bot->SetBottomMargin(0.33);
  bot->SetLeftMargin(0.13);
  bot->SetRightMargin(0.04);
  bot->SetTicks(1, 1);
  top->Draw();
  bot->Draw();

  top->cd();
  Style(ref, kBlue + 1, 20);
  Style(cur, kRed + 1, 21);
  ref->SetTitle("");
  ref->GetXaxis()->SetLabelSize(0);
  ref->GetYaxis()->SetTitle(yTitle);
  ref->GetYaxis()->SetTitleSize(0.055);
  ref->GetYaxis()->SetLabelSize(0.045);
  ref->GetYaxis()->SetTitleOffset(1.1);
  ref->GetYaxis()->SetRangeUser(0.0, 1.2);
  ref->Draw("PE");
  cur->Draw("PE SAME");

  TLegend *leg = BuildLegend(ref, cur);
  leg->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.043);
  latex.DrawLatex(0.15, 0.82, label);

  bot->cd();
  DrawRatioPanel(ref, cur, xTitle);

  gSystem->mkdir(outDir, true);
  canvas->SaveAs(Form("%s/%s.pdf", outDir.Data(), canvasName));
}

void DrawCompareWithRightIntegrated(const char *canvasName,
                                    const char *label,
                                    const char *yTitle,
                                    const char *xTitle,
                                    TH1 *ref,
                                    TH1 *cur,
                                    TH1 *refIntSource,
                                    TH1 *curIntSource,
                                    const TString &outDir)
{
  if (!refIntSource || !curIntSource)
  {
    DrawBasicCompare(canvasName, label, yTitle, xTitle, ref, cur, outDir);
    return;
  }

  TH1D *refInt = CloneAsIntegrated(refIntSource, Form("%s_refInt", canvasName));
  TH1D *curInt = CloneAsIntegrated(curIntSource, Form("%s_curInt", canvasName));
  if (!refInt || !curInt)
  {
    delete refInt;
    delete curInt;
    DrawBasicCompare(canvasName, label, yTitle, xTitle, ref, cur, outDir);
    return;
  }

  const double leftFrac = 0.84;
  TCanvas *canvas = new TCanvas(canvasName, canvasName, 980, 820);
  TPad *left = new TPad(Form("%s_left", canvasName), "", 0.0, 0.0, leftFrac, 1.0);
  TPad *right = new TPad(Form("%s_right", canvasName), "", leftFrac, 0.0, 1.0, 1.0);
  left->SetRightMargin(0.01);
  right->SetLeftMargin(0.0);
  left->Draw();
  right->Draw();

  left->cd();
  TPad *ltop = new TPad(Form("%s_ltop", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *lbot = new TPad(Form("%s_lbot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  ltop->SetBottomMargin(0.03);
  ltop->SetLeftMargin(0.13);
  ltop->SetRightMargin(0.015);
  ltop->SetTicks(1, 1);
  lbot->SetTopMargin(0.03);
  lbot->SetBottomMargin(0.33);
  lbot->SetLeftMargin(0.13);
  lbot->SetRightMargin(0.015);
  lbot->SetTicks(1, 1);
  ltop->Draw();
  lbot->Draw();

  right->cd();
  TPad *rtop = new TPad(Form("%s_rtop", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *rbot = new TPad(Form("%s_rbot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  rtop->SetBottomMargin(0.03);
  rtop->SetLeftMargin(0.02);
  rtop->SetRightMargin(0.08);
  rtop->SetTicks(1, 1);
  rbot->SetTopMargin(0.03);
  rbot->SetBottomMargin(0.33);
  rbot->SetLeftMargin(0.02);
  rbot->SetRightMargin(0.08);
  rbot->SetTicks(1, 1);
  rtop->Draw();
  rbot->Draw();

  ltop->cd();
  Style(ref, kBlue + 1, 20);
  Style(cur, kRed + 1, 21);
  ref->SetTitle("");
  ref->GetXaxis()->SetLabelSize(0);
  ref->GetYaxis()->SetTitle(yTitle);
  ref->GetYaxis()->SetTitleSize(0.055);
  ref->GetYaxis()->SetLabelSize(0.045);
  ref->GetYaxis()->SetTitleOffset(1.1);
  ref->GetYaxis()->SetRangeUser(0.0, 1.2);
  ref->Draw("PE");
  cur->Draw("PE SAME");
  TLegend *leg = BuildLegend(ref, cur);
  leg->Draw();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.043);
  latex.DrawLatex(0.15, 0.82, label);

  rtop->cd();
  Style(refInt, kBlue + 1, 20);
  Style(curInt, kRed + 1, 21);
  TH1D *intFrame = new TH1D(Form("%s_intFrame", canvasName), "", 1, 0.5, 1.5);
  intFrame->SetDirectory(nullptr);
  intFrame->GetXaxis()->SetBinLabel(1, "Int.");
  intFrame->GetYaxis()->SetRangeUser(0.0, 1.2);
  intFrame->GetXaxis()->SetLabelSize(0.0);
  intFrame->GetYaxis()->SetLabelSize(0.0);
  intFrame->Draw("AXIS");
  refInt->Draw("PE SAME");
  curInt->Draw("PE SAME");

  lbot->cd();
  DrawRatioPanel(ref, cur, xTitle);

  rbot->cd();
  TH1D *ratio = new TH1D(Form("%s_intRatio", canvasName), "", 1, 0.5, 1.5);
  ratio->SetDirectory(nullptr);
  ratio->GetXaxis()->SetBinLabel(1, "Int.");
  ratio->GetYaxis()->SetRangeUser(0.70, 1.30);
  ratio->GetXaxis()->SetLabelSize(0.20);
  ratio->GetYaxis()->SetLabelSize(0.0);
  Style(ratio, kBlack, 20);
  if (refInt->GetBinContent(1) != 0.0)
    ratio->SetBinContent(1, curInt->GetBinContent(1) / refInt->GetBinContent(1));
  ratio->Draw("PE");
  TLine *one = new TLine(0.5, 1.0, 1.5, 1.0);
  one->SetLineStyle(2);
  one->Draw("SAME");

  gSystem->mkdir(outDir, true);
  canvas->SaveAs(Form("%s/%s.pdf", outDir.Data(), canvasName));

  delete refInt;
  delete curInt;
}

void CompareOneAcc(bool isPbPb, bool isPrompt)
{
  TFile *fRef = TFile::Open(RefAccFile(isPbPb, isPrompt), "READ");
  TFile *fNew = TFile::Open(NewAccFile(isPbPb, isPrompt), "READ");
  if (!fRef || fRef->IsZombie() || !fNew || fNew->IsZombie())
  {
    std::cout << "[ERROR] failed to open acceptance files: " << RefAccFile(isPbPb, isPrompt)
              << " / " << NewAccFile(isPbPb, isPrompt) << "\n";
    return;
  }

  const TString outDir = OutDir("acceptance");
  const TString tag = Form("acc_%s_%s", CollTag(isPbPb).Data(), StateTag(isPrompt).Data());
  const TString labelBase = Form("%s %s", StateLabel(isPrompt).Data(), isPbPb ? "PbPb" : "pp");
  //const TString labelBase = Form("%s %s acceptance", StateLabel(isPrompt).Data(), isPbPb ? "PbPb" : "pp");

  TH1 *refMid = LoadHist(fRef, "hist_acc_mid", "refAccMid");
  TH1 *refFwd = LoadHist(fRef, "hist_acc_fwd", "refAccFwd");
  TH1 *newMid = LoadHist(fNew, "hAccPt_2021_midy", "newAccMid");
  TH1 *newFwd = LoadHist(fNew, "hAccPt_2021_Fory", "newAccFwd");

  TH1D *refMidInt = BuildIntegratedRatioHist(LoadHist(fRef, "hist_acc_num_mid", "refAccNumMid", false),
                                             LoadHist(fRef, "hist_acc_den_mid", "refAccDenMid", false),
                                             "refAccMidInt");
  TH1D *refFwdInt = BuildIntegratedRatioHist(LoadHist(fRef, "hist_acc_num_fwd", "refAccNumFwd", false),
                                             LoadHist(fRef, "hist_acc_den_fwd", "refAccDenFwd", false),
                                             "refAccFwdInt");
  TH1 *newMidInt = LoadHist(fNew, "hAccPt_2021_midy_Int", "newAccMidInt", false);
  TH1 *newFwdInt = LoadHist(fNew, "hAccPt_2021_Fory_Int", "newAccFwdInt", false);

  if (isPbPb)
  {
    DrawBasicCompare(Form("%s_mid", tag.Data()), Form("%s: |y| < 1.6", labelBase.Data()),
                     "Acceptance", "p_{T} (GeV/c)", refMid, newMid, outDir);
    DrawBasicCompare(Form("%s_fwd", tag.Data()), Form("%s: 1.6 < |y| < 2.4", labelBase.Data()),
                     "Acceptance", "p_{T} (GeV/c)", refFwd, newFwd, outDir);
  }
  else
  {
    DrawCompareWithRightIntegrated(Form("%s_mid", tag.Data()), Form("%s: |y| < 1.6", labelBase.Data()),
                                   "Acceptance", "p_{T} (GeV/c)", refMid, newMid, refMidInt, newMidInt, outDir);
    DrawCompareWithRightIntegrated(Form("%s_fwd", tag.Data()), Form("%s: 1.6 < |y| < 2.4", labelBase.Data()),
                                   "Acceptance", "p_{T} (GeV/c)", refFwd, newFwd, refFwdInt, newFwdInt, outDir);
  }

  fRef->Close();
  fNew->Close();
}

void CompareOneEff(bool isPbPb, bool isPrompt)
{
  TFile *fRef = TFile::Open(RefEffFile(isPbPb, isPrompt), "READ");
  TFile *fNew = TFile::Open(NewEffFile(isPbPb, isPrompt), "READ");
  if (!fRef || fRef->IsZombie() || !fNew || fNew->IsZombie())
  {
    std::cout << "[ERROR] failed to open efficiency files: " << RefEffFile(isPbPb, isPrompt)
              << " / " << NewEffFile(isPbPb, isPrompt) << "\n";
    return;
  }

  const TString outDir = OutDir("efficiency");
  const TString tag = Form("eff_%s_%s", CollTag(isPbPb).Data(), StateTag(isPrompt).Data());
  const TString labelBase = Form("%s %s", StateLabel(isPrompt).Data(), isPbPb ? "PbPb" : "pp");
  //const TString labelBase = Form("%s %s efficiency", StateLabel(isPrompt).Data(), isPbPb ? "PbPb" : "pp");

  TH1 *refMid = LoadHist(fRef, "hist_eff_mid", "refEffMid");
  TH1 *refFwd = LoadHist(fRef, "hist_eff_fwd", "refEffFwd");
  TH1 *newMid = isPbPb ? LoadHist(fNew, "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6", "newEffMid")
                       : LoadHist(fNew, "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6", "newEffMid");
  TH1 *newFwd = isPbPb ? LoadHist(fNew, "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4", "newEffFwd")
                       : LoadHist(fNew, "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4", "newEffFwd");

  if (isPbPb)
  {
    DrawBasicCompare(Form("%s_mid", tag.Data()), Form("%s: |y| < 1.6", labelBase.Data()),
                     "Efficiency", "p_{T} (GeV/c)", refMid, newMid, outDir);
    DrawBasicCompare(Form("%s_fwd", tag.Data()), Form("%s: 1.6 < |y| < 2.4", labelBase.Data()),
                     "Efficiency", "p_{T} (GeV/c)", refFwd, newFwd, outDir);

    TH1 *refCentMid = LoadHist(fRef, "hist_eff_cent_mid", "refEffCentMid", false);
    TH1 *refCentFwd = LoadHist(fRef, "hist_eff_cent_fwd", "refEffCentFwd", false);
    TH1 *newCentMid = LoadHist(fNew, "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6", "newEffCentMid", false);
    TH1 *newCentFwd = LoadHist(fNew, "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4", "newEffCentFwd", false);
    DrawBasicCompare(Form("%s_cent_mid", tag.Data()), Form("%s: |y| < 1.6", labelBase.Data()),
                     "Efficiency", "Centrality", refCentMid, newCentMid, outDir);
    DrawBasicCompare(Form("%s_cent_fwd", tag.Data()), Form("%s: 1.6 < |y| < 2.4", labelBase.Data()),
                     "Efficiency", "Centrality", refCentFwd, newCentFwd, outDir);
  }
  else
  {
    TH1D *refMidInt = BuildIntegratedRatioHist(LoadHist(fRef, "hist_eff_num_mid", "refEffNumMid", false),
                                               LoadHist(fRef, "hist_eff_den_mid", "refEffDenMid", false),
                                               "refEffMidInt");
    TH1D *refFwdInt = BuildIntegratedRatioHist(LoadHist(fRef, "hist_eff_num_fwd", "refEffNumFwd", false),
                                               LoadHist(fRef, "hist_eff_den_fwd", "refEffDenFwd", false),
                                               "refEffFwdInt");
    TH1D *newMidInt = BuildIntegratedRatioHist(LoadHist(fNew, "hist_eff_num_mid", "newEffNumMid", false),
                                               LoadHist(fNew, "hist_eff_den_mid", "newEffDenMid", false),
                                               "newEffMidInt");
    TH1D *newFwdInt = BuildIntegratedRatioHist(LoadHist(fNew, "hist_eff_num_fwd", "newEffNumFwd", false),
                                               LoadHist(fNew, "hist_eff_den_fwd", "newEffDenFwd", false),
                                               "newEffFwdInt");
    DrawCompareWithRightIntegrated(Form("%s_mid", tag.Data()), Form("%s: |y| < 1.6", labelBase.Data()),
                                   "Efficiency", "p_{T} (GeV/c)", refMid, newMid, refMidInt, newMidInt, outDir);
    DrawCompareWithRightIntegrated(Form("%s_fwd", tag.Data()), Form("%s: 1.6 < |y| < 2.4", labelBase.Data()),
                                   "Efficiency", "p_{T} (GeV/c)", refFwd, newFwd, refFwdInt, newFwdInt, outDir);
  }

  fRef->Close();
  fNew->Close();
}
} // namespace

void plot_compare_acc_eff_20260411()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  CompareOneAcc(false, true);
  CompareOneAcc(false, false);
  CompareOneAcc(true, true);
  CompareOneAcc(true, false);

  CompareOneEff(false, true);
  CompareOneEff(false, false);
  CompareOneEff(true, true);
  CompareOneEff(true, false);
}
