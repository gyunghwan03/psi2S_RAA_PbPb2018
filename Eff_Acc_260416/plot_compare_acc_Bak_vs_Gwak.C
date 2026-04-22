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
TString ResolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}

TString NormalizePath(const TString &path)
{
  if (path.IsNull())
    return "";
  if (gSystem->IsAbsoluteFileName(path))
    return path;
  if (!gSystem->AccessPathName(path))
    return path;
  return Form("%s/%s", ResolveBaseDir().Data(), path.Data());
}

TString DefaultStudyRemark(bool isPbPb, bool isPrompt)
{
  if (isPbPb && isPrompt)
    return "20260420";
  if (!isPbPb && isPrompt)
    return "20260420";
  if (isPbPb && !isPrompt)
    return "20260420";
  if (!isPbPb && !isPrompt)
    return "20260420";
  return "";
}

TString BuildBakPath(bool isPbPb,
                     bool isPrompt,
                     const TString &studyRemark,
                     int studyWeightOption,
                     int studyPtWeightSystematic)
{
  if (studyRemark.IsNull())
    return "";

  const TString flavor = isPrompt ? "PromptJpsi" : "BtoJpsi";
  const TString collLabel = isPbPb ? "PbPb" : "pp";
  return Form("%s/roots/acceptance_%s_GenOnly_wgt%d_%s_SysUp%d_%s.root",
              ResolveBaseDir().Data(),
              flavor.Data(),
              studyWeightOption,
              collLabel.Data(),
              studyPtWeightSystematic,
              studyRemark.Data());
}

TString BuildGwakPath(bool isPbPb,
                      bool isPrompt,
                      int acceptanceGenWeight,
                      int acceptancePtWeight)
{
  const TString collLabel = isPbPb ? "PbPb2018" : "pp2018";
  const TString stateLabel = isPrompt ? "PR" : "NP";
  return Form("%s/skim_roots/acc_%s_ppInput_isMC1_%s_ncollW0_genW%d_ptW%d.root",
              ResolveBaseDir().Data(),
              collLabel.Data(),
              stateLabel.Data(),
              acceptanceGenWeight,
              acceptancePtWeight);
}

TString BuildOutDir(const char *outDirPath)
{
  if (outDirPath && outDirPath[0])
    return NormalizePath(outDirPath);
  return Form("%s/plot_outputs_acc_compare_BakGwak", ResolveBaseDir().Data());
}

TString CollisionTag(bool isPbPb)
{
  return isPbPb ? "pbpb" : "pp";
}

TString StateTag(bool isPrompt)
{
  return isPrompt ? "prompt" : "nonprompt";
}

TString BuildCanvasName(const char *regionTag,
                        bool isPbPb,
                        bool isPrompt,
                        const TString &studyRemark)
{
  TString name = Form("compare_acc_pt_%s_Bak_vs_Gwak_%s_%s",
                      regionTag,
                      CollisionTag(isPbPb).Data(),
                      StateTag(isPrompt).Data());
  if (!studyRemark.IsNull())
    name += Form("_%s", studyRemark.Data());
  return name;
}

TString BuildLabelText(bool isPbPb,
                       bool isPrompt,
                       const char *regionText)
{
  const TString stateText = isPrompt ? "Prompt J/#psi" : "Nonprompt J/#psi";
  const TString collText = isPbPb ? "PbPb" : "pp";
  return Form("%s %s: %s", stateText.Data(), collText.Data(), regionText);
}

TH1 *LoadHistClone(TFile *fIn, const char *key, const char *newName, bool required = true)
{
  if (!fIn)
    return nullptr;

  TH1 *h = dynamic_cast<TH1 *>(fIn->Get(key));
  if (!h)
  {
    if (required)
      std::cout << "[ERROR] missing histogram: " << key << "\n";
    else
      std::cout << "[WARN] missing histogram (skip related plot): " << key << "\n";
    return nullptr;
  }

  TH1 *hClone = static_cast<TH1 *>(h->Clone(newName));
  hClone->SetDirectory(nullptr);
  return hClone;
}

bool HasSameBinning(const TH1 *a, const TH1 *b)
{
  if (!a || !b)
    return false;
  if (a->GetNbinsX() != b->GetNbinsX())
    return false;

  for (int i = 1; i <= a->GetNbinsX() + 1; ++i)
  {
    const double edgeA = a->GetXaxis()->GetBinLowEdge(i);
    const double edgeB = b->GetXaxis()->GetBinLowEdge(i);
    if (std::fabs(edgeA - edgeB) > 1e-6)
      return false;
  }
  return true;
}

TH1 *DropLeadingBins(const TH1 *hIn, int nDrop, const char *newName)
{
  if (!hIn)
    return nullptr;
  if (nDrop <= 0)
  {
    TH1 *hClone = static_cast<TH1 *>(hIn->Clone(newName));
    hClone->SetDirectory(nullptr);
    return hClone;
  }
  if (hIn->GetNbinsX() <= nDrop)
    return nullptr;

  const int nOutBins = hIn->GetNbinsX() - nDrop;
  std::vector<double> edges(nOutBins + 1, 0.0);
  for (int i = 0; i <= nOutBins; ++i)
    edges[i] = hIn->GetXaxis()->GetBinLowEdge(i + nDrop + 1);

  TH1D *hOut = new TH1D(newName, hIn->GetTitle(), nOutBins, edges.data());
  hOut->SetDirectory(nullptr);
  hOut->Sumw2();
  for (int i = 1; i <= nOutBins; ++i)
  {
    hOut->SetBinContent(i, hIn->GetBinContent(i + nDrop));
    hOut->SetBinError(i, hIn->GetBinError(i + nDrop));
  }
  return hOut;
}

bool MatchesAfterDroppingLeadingBins(const TH1 *hBak, const TH1 *hGwak, int nDrop)
{
  if (!hBak || !hGwak)
    return false;
  if (nDrop < 0)
    return false;
  if (hBak->GetNbinsX() - nDrop != hGwak->GetNbinsX())
    return false;

  for (int i = 1; i <= hGwak->GetNbinsX() + 1; ++i)
  {
    const double bakEdge = hBak->GetXaxis()->GetBinLowEdge(i + nDrop);
    const double gwakEdge = hGwak->GetXaxis()->GetBinLowEdge(i);
    if (std::fabs(bakEdge - gwakEdge) > 1e-6)
      return false;
  }
  return true;
}

TH1 *AlignBakBinningToGwak(const TH1 *hBak, const TH1 *hGwak, const char *newName)
{
  if (!hBak || !hGwak)
    return nullptr;

  if (HasSameBinning(hBak, hGwak))
    return nullptr;

  if (MatchesAfterDroppingLeadingBins(hBak, hGwak, 1))
  {
    std::cout << "[INFO] auto-align binning by dropping first Bak bin: "
              << hBak->GetName() << " -> " << newName << "\n";
    return DropLeadingBins(hBak, 1, newName);
  }

  return nullptr;
}

void StyleHist(TH1 *h, int color, int markerStyle)
{
  if (!h)
    return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerSize(1.0);
  h->SetLineWidth(2);
}

TH1D *BuildIntegratedViewHist(const TH1 *hSrc, const char *name)
{
  if (!hSrc)
    return nullptr;

  TH1D *hOut = new TH1D(name, "", 1, 0.5, 1.5);
  hOut->SetDirectory(nullptr);
  hOut->SetBinContent(1, hSrc->GetBinContent(1));
  hOut->SetBinError(1, hSrc->GetBinError(1));
  return hOut;
}

void DrawComparisonCanvas(const char *canvasName,
                          const char *labelText,
                          const char *xAxisTitle,
                          TH1 *hBak,
                          TH1 *hGwak,
                          const TString &outDir)
{
  if (!hBak || !hGwak)
    return;

  const bool canDrawRatio = HasSameBinning(hBak, hGwak);

  TCanvas *c = new TCanvas(canvasName, canvasName, 800, 820);
  TPad *padTop = new TPad(Form("%s_top", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *padBot = new TPad(Form("%s_bot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  padTop->SetBottomMargin(0.03);
  padTop->SetLeftMargin(0.13);
  padTop->SetRightMargin(0.04);
  padTop->SetTicks(1, 1);
  padBot->SetTopMargin(0.03);
  padBot->SetBottomMargin(0.33);
  padBot->SetLeftMargin(0.13);
  padBot->SetRightMargin(0.04);
  padBot->SetTicks(1, 1);

  padTop->Draw();
  padBot->Draw();

  padTop->cd();
  StyleHist(hBak, kBlue + 1, 20);
  StyleHist(hGwak, kRed + 1, 21);

  hBak->SetTitle("");
  hBak->GetXaxis()->SetLabelSize(0);
  hBak->GetYaxis()->SetTitle("Acceptance");
  hBak->GetYaxis()->SetTitleSize(0.055);
  hBak->GetYaxis()->SetLabelSize(0.045);
  hBak->GetYaxis()->SetTitleOffset(1.1);
  hBak->GetYaxis()->SetRangeUser(0.0, 1.2);

  hBak->Draw("PE");
  hGwak->Draw("PE SAME");

  TLegend *leg = new TLegend(0.62, 0.72, 0.90, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBak, "Bak (JpsiaccStudy)", "lep");
  leg->AddEntry(hGwak, "Gwak (acceptance_1d)", "lep");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.15, 0.80, labelText);

  padBot->cd();
  if (canDrawRatio)
  {
    TH1 *hRatio = static_cast<TH1 *>(hGwak->Clone(Form("%s_ratio", canvasName)));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hBak);
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(xAxisTitle);
    hRatio->GetYaxis()->SetTitle("Gwak/Bak");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleOffset(0.55);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetRangeUser(0.70, 1.30);
    StyleHist(hRatio, kBlack, 20);
    hRatio->Draw("PE");

    TLine *lineOne = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                               hRatio->GetXaxis()->GetXmax(), 1.0);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(2);
    lineOne->Draw("SAME");
  }
  else
  {
    TH1 *frame = static_cast<TH1 *>(hBak->Clone(Form("%s_ratioFrame", canvasName)));
    frame->Reset("ICES");
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(xAxisTitle);
    frame->GetYaxis()->SetTitle("Gwak/Bak");
    frame->GetXaxis()->SetTitleSize(0.12);
    frame->GetXaxis()->SetLabelSize(0.10);
    frame->GetYaxis()->SetTitleSize(0.10);
    frame->GetYaxis()->SetLabelSize(0.09);
    frame->GetYaxis()->SetTitleOffset(0.55);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetRangeUser(0.70, 1.30);
    frame->Draw();

    TLatex latRatio;
    latRatio.SetNDC();
    latRatio.SetTextSize(0.08);
    latRatio.DrawLatex(0.20, 0.42, "Ratio skipped (bin mismatch)");
  }

  gSystem->mkdir(outDir, true);
  c->SaveAs(Form("%s/%s.pdf", outDir.Data(), canvasName));
  c->SaveAs(Form("%s/%s.png", outDir.Data(), canvasName));
}

void DrawComparisonCanvasWithIntegratedRightPad(const char *canvasName,
                                                const char *labelText,
                                                const char *xAxisTitle,
                                                TH1 *hBakPt,
                                                TH1 *hGwakPt,
                                                TH1 *hBakIntSrc,
                                                TH1 *hGwakIntSrc,
                                                const TString &outDir)
{
  if (!hBakPt || !hGwakPt)
    return;

  if (!hBakIntSrc || !hGwakIntSrc)
  {
    DrawComparisonCanvas(canvasName, labelText, xAxisTitle, hBakPt, hGwakPt, outDir);
    return;
  }

  TH1D *hBakInt = BuildIntegratedViewHist(hBakIntSrc, Form("%s_hBakInt", canvasName));
  TH1D *hGwakInt = BuildIntegratedViewHist(hGwakIntSrc, Form("%s_hGwakInt", canvasName));
  if (!hBakInt || !hGwakInt)
  {
    delete hBakInt;
    delete hGwakInt;
    DrawComparisonCanvas(canvasName, labelText, xAxisTitle, hBakPt, hGwakPt, outDir);
    return;
  }

  const bool canDrawRatioLeft = HasSameBinning(hBakPt, hGwakPt);
  const bool canDrawRatioRight = (hBakInt->GetBinContent(1) != 0.0);
  const double leftFrac = 800.0 / 940.0;

  TCanvas *c = new TCanvas(canvasName, canvasName, 940, 820);
  TPad *padL = new TPad(Form("%s_padL", canvasName), "", 0.00, 0.00, leftFrac, 1.00);
  TPad *padR = new TPad(Form("%s_padR", canvasName), "", leftFrac, 0.00, 1.00, 1.00);
  padL->SetTicks(1, 1);
  padR->SetTicks(1, 1);
  padL->SetRightMargin(0.015);
  padR->SetLeftMargin(0.0);
  padL->Draw();
  padR->Draw();

  padL->cd();
  TPad *padLTop = new TPad(Form("%s_padLTop", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *padLBot = new TPad(Form("%s_padLBot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  padLTop->SetBottomMargin(0.03);
  padLTop->SetLeftMargin(0.13);
  padLTop->SetRightMargin(0.015);
  padLTop->SetTicks(1, 1);
  padLBot->SetTopMargin(0.03);
  padLBot->SetBottomMargin(0.33);
  padLBot->SetLeftMargin(0.13);
  padLBot->SetRightMargin(0.015);
  padLBot->SetTicks(1, 1);
  padLTop->Draw();
  padLBot->Draw();

  padR->cd();
  TPad *padRTop = new TPad(Form("%s_padRTop", canvasName), "", 0.0, 0.30, 1.0, 1.0);
  TPad *padRBot = new TPad(Form("%s_padRBot", canvasName), "", 0.0, 0.0, 1.0, 0.30);
  padRTop->SetBottomMargin(0.03);
  padRTop->SetLeftMargin(0.015);
  padRTop->SetRightMargin(0.06);
  padRTop->SetTicks(1, 1);
  padRBot->SetTopMargin(0.03);
  padRBot->SetBottomMargin(0.33);
  padRBot->SetLeftMargin(0.015);
  padRBot->SetRightMargin(0.06);
  padRBot->SetTicks(1, 1);
  padRTop->Draw();
  padRBot->Draw();

  StyleHist(hBakPt, kBlue + 1, 20);
  StyleHist(hGwakPt, kRed + 1, 21);
  StyleHist(hBakInt, kBlue + 1, 20);
  StyleHist(hGwakInt, kRed + 1, 21);

  padLTop->cd();
  hBakPt->SetTitle("");
  hBakPt->GetXaxis()->SetLabelSize(0);
  hBakPt->GetYaxis()->SetTitle("Acceptance");
  hBakPt->GetYaxis()->SetTitleSize(0.055);
  hBakPt->GetYaxis()->SetLabelSize(0.045);
  hBakPt->GetYaxis()->SetTitleOffset(1.1);
  hBakPt->GetYaxis()->SetRangeUser(0.0, 1.2);
  hBakPt->Draw("PE");
  hGwakPt->Draw("PE SAME");

  TLegend *leg = new TLegend(0.62, 0.72, 0.90, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBakPt, "Bak (JpsiaccStudy)", "lep");
  leg->AddEntry(hGwakPt, "Gwak (acceptance_1d)", "lep");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.15, 0.80, labelText);

  padRTop->cd();
  TH1D *hIntFrame = new TH1D(Form("%s_hIntFrame", canvasName), "", 1, 0.5, 1.5);
  hIntFrame->SetDirectory(nullptr);
  hIntFrame->GetXaxis()->SetBinLabel(1, "Int.");
  hIntFrame->GetYaxis()->SetRangeUser(0.0, 1.2);
  hIntFrame->GetXaxis()->SetLabelSize(0.0);
  hIntFrame->GetXaxis()->SetTitleSize(0.0);
  hIntFrame->GetYaxis()->SetLabelSize(0.0);
  hIntFrame->GetYaxis()->SetTitleSize(0.0);
  hIntFrame->GetYaxis()->SetTickLength(0.08);
  hIntFrame->Draw("AXIS");
  hBakInt->Draw("PE SAME");
  hGwakInt->Draw("PE SAME");

  padLBot->cd();
  if (canDrawRatioLeft)
  {
    TH1 *hRatioL = static_cast<TH1 *>(hGwakPt->Clone(Form("%s_hRatioL", canvasName)));
    hRatioL->SetDirectory(nullptr);
    hRatioL->Divide(hBakPt);
    hRatioL->SetTitle("");
    hRatioL->GetXaxis()->SetTitle(xAxisTitle);
    hRatioL->GetYaxis()->SetTitle("Gwak/Bak");
    hRatioL->GetXaxis()->SetTitleSize(0.12);
    hRatioL->GetXaxis()->SetLabelSize(0.10);
    hRatioL->GetYaxis()->SetTitleSize(0.10);
    hRatioL->GetYaxis()->SetLabelSize(0.09);
    hRatioL->GetYaxis()->SetTitleOffset(0.55);
    hRatioL->GetYaxis()->SetNdivisions(505);
    hRatioL->GetYaxis()->SetRangeUser(0.70, 1.30);
    StyleHist(hRatioL, kBlack, 20);
    hRatioL->Draw("PE");
    TLine *lineOneL = new TLine(hRatioL->GetXaxis()->GetXmin(), 1.0,
                                hRatioL->GetXaxis()->GetXmax(), 1.0);
    lineOneL->SetLineStyle(2);
    lineOneL->SetLineWidth(2);
    lineOneL->Draw("SAME");
  }
  else
  {
    TH1 *frameL = static_cast<TH1 *>(hBakPt->Clone(Form("%s_hRatioLFrame", canvasName)));
    frameL->Reset("ICES");
    frameL->SetTitle("");
    frameL->GetXaxis()->SetTitle(xAxisTitle);
    frameL->GetYaxis()->SetTitle("Gwak/Bak");
    frameL->GetXaxis()->SetTitleSize(0.12);
    frameL->GetXaxis()->SetLabelSize(0.10);
    frameL->GetYaxis()->SetTitleSize(0.10);
    frameL->GetYaxis()->SetLabelSize(0.09);
    frameL->GetYaxis()->SetTitleOffset(0.55);
    frameL->GetYaxis()->SetNdivisions(505);
    frameL->GetYaxis()->SetRangeUser(0.70, 1.30);
    frameL->Draw();
    TLatex latRatio;
    latRatio.SetNDC();
    latRatio.SetTextSize(0.08);
    latRatio.DrawLatex(0.20, 0.42, "Ratio skipped (bin mismatch)");
  }

  padRBot->cd();
  TH1D *hRatioR = new TH1D(Form("%s_hRatioR", canvasName), "", 1, 0.5, 1.5);
  hRatioR->SetDirectory(nullptr);
  hRatioR->GetXaxis()->SetBinLabel(1, "Int.");
  hRatioR->GetYaxis()->SetRangeUser(0.70, 1.30);
  hRatioR->GetXaxis()->SetLabelSize(0.20);
  hRatioR->GetXaxis()->SetTitleSize(0.16);
  hRatioR->GetXaxis()->SetTickLength(0.16);
  hRatioR->GetYaxis()->SetLabelSize(0.0);
  hRatioR->GetYaxis()->SetTitleSize(0.0);
  hRatioR->GetYaxis()->SetTickLength(0.08);
  hRatioR->SetLineColor(kBlack);
  hRatioR->SetLineWidth(2);
  hRatioR->SetMarkerColor(kBlack);
  hRatioR->SetMarkerStyle(20);
  hRatioR->SetMarkerSize(1.0);
  if (canDrawRatioRight)
  {
    hRatioR->SetBinContent(1, hGwakInt->GetBinContent(1) / hBakInt->GetBinContent(1));
    hRatioR->SetBinError(1, 0.0);
  }
  hRatioR->Draw("HIST");
  if (canDrawRatioRight)
    hRatioR->Draw("PE SAME");
  TLine *lineOneR = new TLine(0.5, 1.0, 1.5, 1.0);
  lineOneR->SetLineStyle(2);
  lineOneR->SetLineWidth(2);
  lineOneR->Draw("SAME");
  if (!canDrawRatioRight)
  {
    TLatex latNoRatio;
    latNoRatio.SetNDC();
    latNoRatio.SetTextSize(0.17);
    latNoRatio.DrawLatex(0.18, 0.42, "n/a");
  }

  gSystem->mkdir(outDir, true);
  c->SaveAs(Form("%s/%s.pdf", outDir.Data(), canvasName));
  c->SaveAs(Form("%s/%s.png", outDir.Data(), canvasName));

  delete hRatioR;
  delete hIntFrame;
  delete hBakInt;
  delete hGwakInt;
}
} // namespace

bool plot_compare_acc_Bak_vs_Gwak(bool isPbPb = true,
                                  bool isPrompt = true,
                                  const char *studyRemark = "",
                                  int acceptanceGenWeight = 1,
                                  int acceptancePtWeight = 1,
                                  int studyWeightOption = 1,
                                  int studyPtWeightSystematic = 0,
                                  const char *bakFilePath = "",
                                  const char *gwakFilePath = "",
                                  const char *outDirPath = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TString remark = (studyRemark && studyRemark[0]) ? TString(studyRemark) : DefaultStudyRemark(isPbPb, isPrompt);
  const TString bakPath = (bakFilePath && bakFilePath[0])
                              ? NormalizePath(bakFilePath)
                              : BuildBakPath(isPbPb, isPrompt, remark, studyWeightOption, studyPtWeightSystematic);
  const TString gwakPath = (gwakFilePath && gwakFilePath[0])
                               ? NormalizePath(gwakFilePath)
                               : BuildGwakPath(isPbPb, isPrompt, acceptanceGenWeight, acceptancePtWeight);
  const TString outDir = BuildOutDir(outDirPath);

  if (bakPath.IsNull())
  {
    std::cout << "[ERROR] Bak file path is empty. Pass studyRemark or bakFilePath.\n";
    return false;
  }

  TFile *fBak = TFile::Open(bakPath, "READ");
  TFile *fGwak = TFile::Open(gwakPath, "READ");
  if (!fBak || fBak->IsZombie())
  {
    std::cout << "[ERROR] failed to open Bak file: " << bakPath << "\n";
    return false;
  }
  if (!fGwak || fGwak->IsZombie())
  {
    std::cout << "[ERROR] failed to open Gwak file: " << gwakPath << "\n";
    fBak->Close();
    delete fBak;
    return false;
  }

  TH1 *hBakPtAll = LoadHistClone(fBak, "hAccPt_2021_ally", "hBakPtAll", false);
  TH1 *hBakPtFwd = LoadHistClone(fBak, "hAccPt_2021_Fory", "hBakPtFwd");
  TH1 *hBakPtMid = LoadHistClone(fBak, "hAccPt_2021_midy", "hBakPtMid");
  TH1 *hBakPtAllInt = LoadHistClone(fBak, "hAccPt_2021_ally_Int", "hBakPtAllInt", false);
  TH1 *hBakPtFwdInt = LoadHistClone(fBak, "hAccPt_2021_Fory_Int", "hBakPtFwdInt", false);
  TH1 *hBakPtMidInt = LoadHistClone(fBak, "hAccPt_2021_midy_Int", "hBakPtMidInt", false);
  TH1 *hGwakPtAll = LoadHistClone(fGwak, "hist_acc_all", "hGwakPtAll", false);
  TH1 *hGwakPtFwd = LoadHistClone(fGwak, "hist_acc_fwd", "hGwakPtFwd");
  TH1 *hGwakPtMid = LoadHistClone(fGwak, "hist_acc_mid", "hGwakPtMid");
  TH1 *hGwakPtAllInt = LoadHistClone(fGwak, "hist_acc_all_Int", "hGwakPtAllInt", false);
  TH1 *hGwakPtFwdInt = LoadHistClone(fGwak, "hist_acc_fwd_Int", "hGwakPtFwdInt", false);
  TH1 *hGwakPtMidInt = LoadHistClone(fGwak, "hist_acc_mid_Int", "hGwakPtMidInt", false);

  bool ok = (hBakPtFwd && hBakPtMid && hGwakPtFwd && hGwakPtMid);
  if (!ok)
  {
    fBak->Close();
    fGwak->Close();
    delete fBak;
    delete fGwak;
    delete hBakPtAll;
    delete hBakPtFwd;
    delete hBakPtMid;
    delete hBakPtAllInt;
    delete hBakPtFwdInt;
    delete hBakPtMidInt;
    delete hGwakPtAll;
    delete hGwakPtFwd;
    delete hGwakPtMid;
    delete hGwakPtAllInt;
    delete hGwakPtFwdInt;
    delete hGwakPtMidInt;
    return false;
  }

  if (hBakPtAll && hGwakPtAll)
  {
    if (TH1 *hAligned = AlignBakBinningToGwak(hBakPtAll, hGwakPtAll, "hBakPtAll_aligned"))
    {
      delete hBakPtAll;
      hBakPtAll = hAligned;
    }
  }
  if (TH1 *hAligned = AlignBakBinningToGwak(hBakPtFwd, hGwakPtFwd, "hBakPtFwd_aligned"))
  {
    delete hBakPtFwd;
    hBakPtFwd = hAligned;
  }
  if (TH1 *hAligned = AlignBakBinningToGwak(hBakPtMid, hGwakPtMid, "hBakPtMid_aligned"))
  {
    delete hBakPtMid;
    hBakPtMid = hAligned;
  }

  const TString canvasNameAll = BuildCanvasName("all", isPbPb, isPrompt, remark);
  const TString canvasNameFwd = BuildCanvasName("fwd", isPbPb, isPrompt, remark);
  const TString canvasNameMid = BuildCanvasName("mid", isPbPb, isPrompt, remark);
  if (hBakPtAll && hGwakPtAll)
  {
    DrawComparisonCanvasWithIntegratedRightPad(canvasNameAll.Data(),
                                               BuildLabelText(isPbPb, isPrompt, "|y| < 2.4").Data(),
                                               "p_{T} (GeV/c)",
                                               hBakPtAll, hGwakPtAll,
                                               hBakPtAllInt, hGwakPtAllInt,
                                               outDir);
  }

  DrawComparisonCanvasWithIntegratedRightPad(canvasNameFwd.Data(),
                                             BuildLabelText(isPbPb, isPrompt, "1.6 < |y| < 2.4").Data(),
                                             "p_{T} (GeV/c)",
                                             hBakPtFwd, hGwakPtFwd,
                                             hBakPtFwdInt, hGwakPtFwdInt,
                                             outDir);

  DrawComparisonCanvasWithIntegratedRightPad(canvasNameMid.Data(),
                                             BuildLabelText(isPbPb, isPrompt, "|y| < 1.6").Data(),
                                             "p_{T} (GeV/c)",
                                             hBakPtMid, hGwakPtMid,
                                             hBakPtMidInt, hGwakPtMidInt,
                                             outDir);

  fBak->Close();
  fGwak->Close();
  delete fBak;
  delete fGwak;
  delete hBakPtAll;
  delete hBakPtFwd;
  delete hBakPtMid;
  delete hBakPtAllInt;
  delete hBakPtFwdInt;
  delete hBakPtMidInt;
  delete hGwakPtAll;
  delete hGwakPtFwd;
  delete hGwakPtMid;
  delete hGwakPtAllInt;
  delete hGwakPtFwdInt;
  delete hGwakPtMidInt;

  std::cout << "[INFO] Bak  -> JpsiaccStudy   : " << bakPath << "\n";
  std::cout << "[INFO] Gwak -> acceptance_1d : " << gwakPath << "\n";
  std::cout << "[INFO] saved comparison plots to: " << outDir << "\n";
  return true;
}

void plot_compare_acc_pbpb_prompt_Bak_vs_Gwak(const char *studyRemark = "20260420")
{
  plot_compare_acc_Bak_vs_Gwak(true, true, studyRemark);
}

void plot_compare_acc_pp_prompt_Bak_vs_Gwak(const char *studyRemark = "20260419")
{
  plot_compare_acc_Bak_vs_Gwak(false, true, studyRemark);
}

void plot_compare_acc_pbpb_nonprompt_Bak_vs_Gwak(const char *studyRemark,
                                                 int acceptanceGenWeight = 1,
                                                 int acceptancePtWeight = 1,
                                                 int studyWeightOption = 1,
                                                 int studyPtWeightSystematic = 0)
{
  plot_compare_acc_Bak_vs_Gwak(true,
                               false,
                               studyRemark,
                               acceptanceGenWeight,
                               acceptancePtWeight,
                               studyWeightOption,
                               studyPtWeightSystematic);
}

void plot_compare_acc_pp_nonprompt_Bak_vs_Gwak(const char *studyRemark,
                                               int acceptanceGenWeight = 1,
                                               int acceptancePtWeight = 1,
                                               int studyWeightOption = 1,
                                               int studyPtWeightSystematic = 0)
{
  plot_compare_acc_Bak_vs_Gwak(false,
                               false,
                               studyRemark,
                               acceptanceGenWeight,
                               acceptancePtWeight,
                               studyWeightOption,
                               studyPtWeightSystematic);
}
