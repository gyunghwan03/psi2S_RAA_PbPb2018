#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
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
#include "../CMS_lumi_v2mass.C"

namespace
{
TString ResolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}

TH1 *LoadHistClone(TFile *fIn, const char *key, const char *newName)
{
  if (!fIn)
    return nullptr;

  TH1 *h = dynamic_cast<TH1 *>(fIn->Get(key));
  if (!h)
  {
    std::cout << "[ERROR] missing histogram: " << key << "\n";
    return nullptr;
  }

  TH1 *hClone = static_cast<TH1 *>(h->Clone(newName));
  hClone->SetDirectory(nullptr);
  return hClone;
}

TH1 *LoadHistClone(TFile *fIn, const char *key, const char *newName, bool required)
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

bool MatchesAfterDroppingLeadingBins(const TH1 *hBak, const TH1 *hRef, int nDrop)
{
  if (!hBak || !hRef)
    return false;
  if (nDrop < 0)
    return false;
  if (hBak->GetNbinsX() - nDrop != hRef->GetNbinsX())
    return false;

  for (int i = 1; i <= hRef->GetNbinsX() + 1; ++i)
  {
    const double bakEdge = hBak->GetXaxis()->GetBinLowEdge(i + nDrop);
    const double refEdge = hRef->GetXaxis()->GetBinLowEdge(i);
    if (std::fabs(bakEdge - refEdge) > 1e-6)
      return false;
  }
  return true;
}

TH1 *AlignBakBinningToGwak(const TH1 *hBak, const TH1 *hGwak, const char *newName)
{
  if (!hBak || !hGwak)
    return nullptr;

  if (HasSameBinning(hBak, hGwak))
  {
    TH1 *hClone = static_cast<TH1 *>(hBak->Clone(newName));
    hClone->SetDirectory(nullptr);
    return hClone;
  }

  // Common mismatch here: Bak histogram keeps one low-pt leading bin
  // (e.g. 0-6.5 or 0-3.5) while Gwak starts from physics fiducial lower edge.
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

  if (hSrc->GetNbinsX() == 1)
  {
    hOut->SetBinContent(1, hSrc->GetBinContent(1));
    hOut->SetBinError(1, hSrc->GetBinError(1));
    return hOut;
  }

  double sum = 0.0;
  double err2 = 0.0;
  for (int ib = 1; ib <= hSrc->GetNbinsX(); ++ib)
  {
    sum += hSrc->GetBinContent(ib);
    err2 += std::pow(hSrc->GetBinError(ib), 2);
  }
  hOut->SetBinContent(1, sum);
  hOut->SetBinError(1, std::sqrt(err2));
  return hOut;
}

TH1D *BuildIntegratedEfficiencyHist(const TH1 *hNum, const TH1 *hDen, const char *name)
{
  if (!hNum || !hDen)
    return nullptr;

  TH1D *hOut = new TH1D(name, "", 1, 0.5, 1.5);
  hOut->SetDirectory(nullptr);

  double sumNum = 0.0;
  double err2Num = 0.0;
  for (int ib = 1; ib <= hNum->GetNbinsX(); ++ib)
  {
    sumNum += hNum->GetBinContent(ib);
    err2Num += std::pow(hNum->GetBinError(ib), 2);
  }

  double sumDen = 0.0;
  double err2Den = 0.0;
  for (int ib = 1; ib <= hDen->GetNbinsX(); ++ib)
  {
    sumDen += hDen->GetBinContent(ib);
    err2Den += std::pow(hDen->GetBinError(ib), 2);
  }

  if (sumDen <= 0.0)
  {
    hOut->SetBinContent(1, 0.0);
    hOut->SetBinError(1, 0.0);
    return hOut;
  }

  const double eff = sumNum / sumDen;
  double errEff = 0.0;
  if (sumNum > 0.0)
  {
    const double relNum = std::sqrt(err2Num) / sumNum;
    const double relDen = std::sqrt(err2Den) / sumDen;
    errEff = eff * std::sqrt(relNum * relNum + relDen * relDen);
  }
  hOut->SetBinContent(1, eff);
  hOut->SetBinError(1, errEff);
  return hOut;
}

void DrawComparisonCanvas(const char *canvasName,
                          const char *labelText,
                          const char *xAxisTitle,
                          TH1 *hBak,
                          TH1 *hGwak,
                          const TString &outDir,
                          int iPeriod)
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
  hBak->GetYaxis()->SetTitle("Efficiency");
  hBak->GetYaxis()->SetTitleSize(0.055);
  hBak->GetYaxis()->SetLabelSize(0.045);
  hBak->GetYaxis()->SetTitleOffset(1.1);

  const double maxY = std::max(hBak->GetMaximum(), hGwak->GetMaximum());
  //hBak->GetYaxis()->SetRangeUser(0.0, std::max(1.2 * maxY, 0.05));
  hBak->GetYaxis()->SetRangeUser(0.0, 1.2);

  hBak->Draw("PE");
  hGwak->Draw("PE SAME");

  TLegend *leg = new TLegend(0.21, 0.62, 0.35, 0.76);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBak, "Bak v7", "lep");
  leg->AddEntry(hGwak, "Gwak", "lep");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.15, 0.8, labelText);

  int iPos = 33;
  CMS_lumi_v2mass(padTop,iPeriod,iPos);

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
  //c->SaveAs(Form("%s/%s.png", outDir.Data(), canvasName));
}

void DrawComparisonCanvasWithIntegratedRightPad(const char *canvasName,
                                                const char *labelText,
                                                const char *xAxisTitle,
                                                TH1 *hBakPt,
                                                TH1 *hGwakPt,
                                                TH1 *hBakIntSrc,
                                                TH1 *hGwakIntSrc,
                                                const TString &outDir,
                                                int iPeriod)
{
  if (!hBakPt || !hGwakPt)
    return;

  if (!hBakIntSrc || !hGwakIntSrc)
  {
    DrawComparisonCanvas(canvasName, labelText, xAxisTitle, hBakPt, hGwakPt, outDir, iPeriod);
    return;
  }

  TH1D *hBakInt = BuildIntegratedViewHist(hBakIntSrc, Form("%s_hBakInt", canvasName));
  TH1D *hGwakInt = BuildIntegratedViewHist(hGwakIntSrc, Form("%s_hGwakInt", canvasName));
  if (!hBakInt || !hGwakInt)
  {
    delete hBakInt;
    delete hGwakInt;
    DrawComparisonCanvas(canvasName, labelText, xAxisTitle, hBakPt, hGwakPt, outDir, iPeriod);
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
  hBakPt->GetYaxis()->SetTitle("Efficiency");
  hBakPt->GetYaxis()->SetTitleSize(0.055);
  hBakPt->GetYaxis()->SetLabelSize(0.045);
  hBakPt->GetYaxis()->SetTitleOffset(1.1);
  hBakPt->GetYaxis()->SetRangeUser(0.0, 1.2);
  hBakPt->Draw("PE");
  hGwakPt->Draw("PE SAME");

  TLegend *leg = new TLegend(0.21, 0.62, 0.35, 0.76);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBakPt, "Bak v7", "lep");
  leg->AddEntry(hGwakPt, "Gwak", "lep");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.15, 0.8, labelText);

  int iPos = 33;
  CMS_lumi_v2mass(padLTop, iPeriod, iPos);

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
  //c->SaveAs(Form("%s/%s.png", outDir.Data(), canvasName));

  delete hRatioR;
  delete hIntFrame;
  delete hBakInt;
  delete hGwakInt;
}
} // namespace

void plot_compare_eff_Bak_vs_Gwak(
    const char *bakFilePath = "",
    const char *gwakFilePath = "",
    const char *outDirPath = "",
    const char *ppBakFilePath = "",
    const char *ppGwakFilePath = "",
    int state = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const TString baseDir = ResolveBaseDir();
  const bool isPrompt = (state == 1);
  const TString sampleTag = isPrompt ? "prompt" : "nprompt";
  const TString sampleLabel = isPrompt ? "Prompt J/#psi PbPb" : "Nonprompt J/#psi PbPb";
  const TString gwakTag = isPrompt ? "PR" : "NP";
  const TString bakPath = (bakFilePath && bakFilePath[0])
                              ? TString(bakFilePath)
                              : Form("%s/roots/mc_eff_vs_pt_cent_0_to_180_rap_%s_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260422_v7.root",
                                     baseDir.Data(), sampleTag.Data());
  const TString gwakPath = (gwakFilePath && gwakFilePath[0])
                               ? TString(gwakFilePath)
                               : Form("%s/skim_roots/eff_PbPb2018_isMC1_%s_ncollW1_genW1_ptW1_tnpW1_ctauW1_Dimuon_MiniAOD_2exp.root",
                                      baseDir.Data(), gwakTag.Data());
  const TString outDir = (outDirPath && outDirPath[0])
                             ? TString(outDirPath)
                             : Form("%s/plot_outputs_eff_compare_v7_vs_Gwak_%s", baseDir.Data(), gwakTag.Data());
  const TString ppBakPath = (ppBakFilePath && ppBakFilePath[0])
                                ? TString(ppBakFilePath)
                                 : Form("%s/roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_260420_2exp.root",
                                       baseDir.Data());
  const TString ppGwakPath = (ppGwakFilePath && ppGwakFilePath[0])
                                 ? TString(ppGwakFilePath)
                                 : Form("%s/skim_roots/eff_pp5p02TeV_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_ctauW1_Dimuon_MiniAOD.root",
                                        baseDir.Data());

  TFile *fBak = TFile::Open(bakPath, "READ");
  TFile *fGwak = TFile::Open(gwakPath, "READ");
  if (!fBak || fBak->IsZombie())
  {
    std::cout << "[ERROR] failed to open Bak file: " << bakPath << "\n";
    return;
  }
  if (!fGwak || fGwak->IsZombie())
  {
    std::cout << "[ERROR] failed to open Gwak file: " << gwakPath << "\n";
    fBak->Close();
    return;
  }

  TH1 *hBakPtFwd = LoadHistClone(fBak, "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4", "hBakPtFwd");
  TH1 *hBakPtMid = LoadHistClone(fBak, "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6", "hBakPtMid");
  TH1 *hBakCentFwd = LoadHistClone(fBak, "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4", "hBakCentFwd");
  TH1 *hBakCentMid = LoadHistClone(fBak, "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6", "hBakCentMid");

  TH1 *hGwakPtFwd = LoadHistClone(fGwak, "hist_eff_fwd", "hGwakPtFwd");
  TH1 *hGwakPtMid = LoadHistClone(fGwak, "hist_eff_mid", "hGwakPtMid");
  TH1 *hGwakCentFwd = LoadHistClone(fGwak, "hist_eff_cent_fwd", "hGwakCentFwd");
  TH1 *hGwakCentMid = LoadHistClone(fGwak, "hist_eff_cent_mid", "hGwakCentMid");

  // Auto-align Bak pT histograms when they contain an extra leading low-pT bin.
  if (TH1 *hAligned = AlignBakBinningToGwak(hBakPtFwd, hGwakPtFwd, "hBakPtFwd_aligned"))
    hBakPtFwd = hAligned;
  if (TH1 *hAligned = AlignBakBinningToGwak(hBakPtMid, hGwakPtMid, "hBakPtMid_aligned"))
    hBakPtMid = hAligned;

  DrawComparisonCanvas("compare_eff_pt_fwd_v7_vs_Gwak",
                       Form("%s: 1.6 < |y| < 2.4", sampleLabel.Data()),
                       "p_{T} (GeV/c)",
                       hBakPtFwd, hGwakPtFwd, outDir, 2);

  DrawComparisonCanvas("compare_eff_pt_mid_v7_vs_Gwak",
                       Form("%s: |y| < 1.6", sampleLabel.Data()),
                       "p_{T} (GeV/c)",
                       hBakPtMid, hGwakPtMid, outDir, 2);

  DrawComparisonCanvas("compare_eff_cent_fwd_v7_vs_Gwak",
                       Form("%s: 1.6 < |y| < 2.4, 3.5 < p_{T} < 40", sampleLabel.Data()),
                       "Centrality",
                       hBakCentFwd, hGwakCentFwd, outDir, 2);

  DrawComparisonCanvas("compare_eff_cent_mid_v7_vs_Gwak",
                       Form("%s: |y| < 1.6, 6.5 < p_{T} < 40", sampleLabel.Data()),
                       "Centrality",
                       hBakCentMid, hGwakCentMid, outDir, 2);

  if (isPrompt)
  {
    TFile *fBakPp = TFile::Open(ppBakPath, "READ");
    TFile *fGwakPp = TFile::Open(ppGwakPath, "READ");
    if (!fBakPp || fBakPp->IsZombie())
    {
      std::cout << "[WARN] failed to open pp Bak file. Skip pp comparison: " << ppBakPath << "\n";
    }
    else if (!fGwakPp || fGwakPp->IsZombie())
    {
      std::cout << "[WARN] failed to open pp Gwak file. Skip pp comparison: " << ppGwakPath << "\n";
      fBakPp->Close();
    }
    else
    {
      TH1 *hBakPpPtFwd = LoadHistClone(fBakPp, "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4", "hBakPpPtFwd");
      TH1 *hBakPpPtMid = LoadHistClone(fBakPp, "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6", "hBakPpPtMid");
      TH1 *hBakPpIntFwd = LoadHistClone(fBakPp, "mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4", "hBakPpIntFwd", false);
      TH1 *hBakPpIntMid = LoadHistClone(fBakPp, "mc_eff_Integrated_TnP1_PtW1_absy0_1p6", "hBakPpIntMid", false);

      TH1 *hGwakPpPtFwd = LoadHistClone(fGwakPp, "hist_eff_fwd", "hGwakPpPtFwd");
      TH1 *hGwakPpPtMid = LoadHistClone(fGwakPp, "hist_eff_mid", "hGwakPpPtMid");
      TH1 *hGwakPpIntFwd = LoadHistClone(fGwakPp, "hist_eff_fwd_Int", "hGwakPpIntFwd", false);
      TH1 *hGwakPpIntMid = LoadHistClone(fGwakPp, "hist_eff_mid_Int", "hGwakPpIntMid", false);
      TH1 *hGwakPpNumFwd = nullptr;
      TH1 *hGwakPpDenFwd = nullptr;
      TH1 *hGwakPpNumMid = nullptr;
      TH1 *hGwakPpDenMid = nullptr;
      TH1D *hGwakPpIntFwdBuilt = nullptr;
      TH1D *hGwakPpIntMidBuilt = nullptr;

      if (!hGwakPpIntFwd || !hGwakPpIntMid)
      {
        hGwakPpNumFwd = LoadHistClone(fGwakPp, "hist_eff_num_fwd", "hGwakPpNumFwd", false);
        hGwakPpDenFwd = LoadHistClone(fGwakPp, "hist_eff_den_fwd", "hGwakPpDenFwd", false);
        hGwakPpNumMid = LoadHistClone(fGwakPp, "hist_eff_num_mid", "hGwakPpNumMid", false);
        hGwakPpDenMid = LoadHistClone(fGwakPp, "hist_eff_den_mid", "hGwakPpDenMid", false);

        if (!hGwakPpIntFwd)
        {
          hGwakPpIntFwdBuilt = BuildIntegratedEfficiencyHist(hGwakPpNumFwd, hGwakPpDenFwd, "hGwakPpIntFwdBuilt");
          hGwakPpIntFwd = hGwakPpIntFwdBuilt;
        }
        if (!hGwakPpIntMid)
        {
          hGwakPpIntMidBuilt = BuildIntegratedEfficiencyHist(hGwakPpNumMid, hGwakPpDenMid, "hGwakPpIntMidBuilt");
          hGwakPpIntMid = hGwakPpIntMidBuilt;
        }
      }

      if (TH1 *hAligned = AlignBakBinningToGwak(hBakPpPtFwd, hGwakPpPtFwd, "hBakPpPtFwd_aligned"))
        hBakPpPtFwd = hAligned;
      if (TH1 *hAligned = AlignBakBinningToGwak(hBakPpPtMid, hGwakPpPtMid, "hBakPpPtMid_aligned"))
        hBakPpPtMid = hAligned;

      DrawComparisonCanvasWithIntegratedRightPad("compare_eff_pp_pt_fwd_v7_vs_Gwak",
                                                 "Prompt J/#psi pp: 1.6 < |y| < 2.4",
                                                 "p_{T} (GeV/c)",
                                                 hBakPpPtFwd, hGwakPpPtFwd,
                                                 hBakPpIntFwd, hGwakPpIntFwd,
                                                 outDir, 1);

      DrawComparisonCanvasWithIntegratedRightPad("compare_eff_pp_pt_mid_v7_vs_Gwak",
                                                 "Prompt J/#psi pp: |y| < 1.6",
                                                 "p_{T} (GeV/c)",
                                                 hBakPpPtMid, hGwakPpPtMid,
                                                 hBakPpIntMid, hGwakPpIntMid,
                                                 outDir, 1);

      fBakPp->Close();
      fGwakPp->Close();
      delete hBakPpIntFwd;
      delete hBakPpIntMid;
      delete hGwakPpIntFwdBuilt;
      delete hGwakPpIntMidBuilt;
      delete hGwakPpNumFwd;
      delete hGwakPpDenFwd;
      delete hGwakPpNumMid;
      delete hGwakPpDenMid;
      if (hGwakPpIntFwd != hGwakPpIntFwdBuilt)
        delete hGwakPpIntFwd;
      if (hGwakPpIntMid != hGwakPpIntMidBuilt)
        delete hGwakPpIntMid;
    }
  }
  else
  {
    std::cout << "[INFO] state=2 (nonprompt): skip pp comparison block.\n";
  }

  fBak->Close();
  fGwak->Close();
  std::cout << "[INFO] saved comparison plots to: " << outDir << "\n";
}
