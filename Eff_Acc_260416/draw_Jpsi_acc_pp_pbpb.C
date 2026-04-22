#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <iostream>

namespace
{
TString GetBaseDir()
{
  return gSystem->DirName(__FILE__);
}

TString ResolvePath(const TString &path)
{
  if (path.IsNull())
    return "";
  if (gSystem->IsAbsoluteFileName(path))
    return path;
  if (!gSystem->AccessPathName(path))
    return path;
  return Form("%s/%s", GetBaseDir().Data(), path.Data());
}

TString BuildDefaultInputPath(bool isPbPb,
                              bool isPrompt,
                              int weightOption,
                              int ptWeightSystematic,
                              const TString &remark)
{
  const TString flavor = isPrompt ? "PromptJpsi" : "BtoJpsi";
  const TString collLabel = isPbPb ? "PbPb" : "pp";
  return Form("%s/roots/acceptance_%s_GenOnly_wgt%d_%s_SysUp%d_%s.root",
              GetBaseDir().Data(),
              flavor.Data(),
              weightOption,
              collLabel.Data(),
              ptWeightSystematic,
              remark.Data());
}

TH1 *LoadHistClone(TFile *file, const char *key, const char *cloneName)
{
  if (!file)
    return nullptr;

  TH1 *hist = dynamic_cast<TH1 *>(file->Get(key));
  if (!hist)
  {
    std::cout << "[ERROR] missing histogram: " << key << "\n";
    return nullptr;
  }

  TH1 *clone = static_cast<TH1 *>(hist->Clone(cloneName));
  clone->SetDirectory(nullptr);
  return clone;
}

void StyleHist(TH1 *hist, int color, int markerStyle)
{
  if (!hist)
    return;

  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(1.1);
  hist->SetLineWidth(2);
}

double GetSafeMaximum(TH1 *a, TH1 *b, TH1 *c, TH1 *d)
{
  double ymax = 0.0;
  if (a)
    ymax = std::max(ymax, a->GetMaximum());
  if (b)
    ymax = std::max(ymax, b->GetMaximum());
  if (c)
    ymax = std::max(ymax, c->GetMaximum());
  if (d)
    ymax = std::max(ymax, d->GetMaximum());
  return std::max(0.2, ymax);
}

TH1D *BuildIntegratedViewHist(const TH1 *histSrc, const char *name)
{
  if (!histSrc)
    return nullptr;

  TH1D *histOut = new TH1D(name, "", 1, 0.5, 1.5);
  histOut->SetDirectory(nullptr);
  histOut->SetBinContent(1, histSrc->GetBinContent(1));
  histOut->SetBinError(1, histSrc->GetBinError(1));
  return histOut;
}

void StyleGraph(TGraphErrors *graph, int color, int markerStyle)
{
  if (!graph)
    return;

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(1.3);
  graph->SetLineWidth(2);
}

TGraphErrors *BuildIntegratedGraph(const TH1 *histSrc,
                                   const char *name,
                                   double xValue,
                                   double xError,
                                   int color,
                                   int markerStyle)
{
  if (!histSrc)
    return nullptr;

  TGraphErrors *graph = new TGraphErrors(1);
  graph->SetName(name);
  graph->SetPoint(0, xValue, histSrc->GetBinContent(1));
  graph->SetPointError(0, xError, histSrc->GetBinError(1));
  StyleGraph(graph, color, markerStyle);
  return graph;
}

void SetFrameStyle(TH1D *frame,
                   const char *xTitle,
                   const char *yTitle,
                   double yMax)
{
  if (!frame)
    return;

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleSize(0.052);
  frame->GetYaxis()->SetTitleSize(0.052);
  frame->GetXaxis()->SetLabelSize(0.044);
  frame->GetYaxis()->SetLabelSize(0.044);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetRangeUser(0.0, yMax);
  frame->SetStats(0);
}

void DrawRegionCanvasWithIntegratedPad(const char *canvasName,
                                       const char *regionLabel,
                                       TH1 *hPpPt,
                                       TH1 *hPbPbPt,
                                       TH1 *hPpIntSrc,
                                       TH1 *hPbPbIntSrc,
                                       double xMin,
                                       double xMax,
                                       bool isPrompt,
                                       const TString &outStem)
{
  if (!hPpPt || !hPbPbPt || !hPpIntSrc || !hPbPbIntSrc)
    return;

  TH1D *hPpInt = BuildIntegratedViewHist(hPpIntSrc, Form("%s_ppIntHist", canvasName));
  TH1D *hPbPbInt = BuildIntegratedViewHist(hPbPbIntSrc, Form("%s_pbpbIntHist", canvasName));
  if (!hPpInt || !hPbPbInt)
  {
    delete hPpInt;
    delete hPbPbInt;
    return;
  }

  const double integratedHalfWidth = 0.018;
  TGraphErrors *gPpInt = BuildIntegratedGraph(hPpIntSrc, Form("%s_ppIntGraph", canvasName), 0.96, integratedHalfWidth, kBlue + 1, 20);
  TGraphErrors *gPbPbInt = BuildIntegratedGraph(hPbPbIntSrc, Form("%s_pbpbIntGraph", canvasName), 1.04, integratedHalfWidth, kRed + 1, 21);
  if (!gPpInt || !gPbPbInt)
  {
    delete hPpInt;
    delete hPbPbInt;
    delete gPpInt;
    delete gPbPbInt;
    return;
  }

  //const double yMax = std::min(1.05, 1.25 * GetSafeMaximum(hPpPt, hPbPbPt, hPpInt, hPbPbInt));
  double yMax = 0.83;

  TCanvas *canvas = new TCanvas(canvasName, "", 980, 560);
  TPad *padLeft = new TPad(Form("%s_left", canvasName), "", 0.00, 0.00, 0.82, 0.90);
  TPad *padRight = new TPad(Form("%s_right", canvasName), "", 0.82, 0.00, 1.00, 0.90);
  padLeft->SetLeftMargin(0.15);
  padLeft->SetRightMargin(0.02);
  padLeft->SetBottomMargin(0.14);
  padLeft->SetTopMargin(0.08);
  padRight->SetLeftMargin(0.02);
  padRight->SetRightMargin(0.10);
  padRight->SetBottomMargin(0.14);
  padRight->SetTopMargin(0.08);
  padLeft->SetTicks(1, 1);
  padRight->SetTicks(1, 1);
  padLeft->Draw();
  padRight->Draw();

  padLeft->cd();
  TH1D *frameLeft = new TH1D(Form("%s_frameLeft", canvasName), "", 100, xMin, xMax);
  SetFrameStyle(frameLeft, "p_{T} (GeV/c)", "Acceptance", yMax);
  frameLeft->Draw();

  hPpPt->Draw("PE SAME");
  hPbPbPt->Draw("PE SAME");

  TLatex latexLeft;
  latexLeft.SetNDC();
  latexLeft.SetTextSize(0.050);
  latexLeft.DrawLatex(0.16, 0.84, regionLabel);

  TLegend *legend = new TLegend(0.58, 0.69, 0.90, 0.84);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.040);
  legend->AddEntry(hPpPt, "pp", "lep");
  legend->AddEntry(hPbPbPt, "PbPb", "lep");
  legend->Draw();

  padRight->cd();
  TH1D *frameRight = new TH1D(Form("%s_frameRight", canvasName), "", 1, 0.5, 1.5);
  frameRight->SetDirectory(nullptr);
  frameRight->SetTitle("");
  frameRight->GetXaxis()->SetBinLabel(1, "Int.");
  frameRight->GetXaxis()->SetLabelSize(0.16);
  frameRight->GetXaxis()->SetTickLength(0.08);
  frameRight->GetYaxis()->SetRangeUser(0.0, yMax);
  frameRight->GetYaxis()->SetLabelSize(0.0);
  frameRight->GetYaxis()->SetTickLength(0.06);
  frameRight->GetYaxis()->SetTitleSize(0.0);
  frameRight->Draw();

  gPpInt->Draw("PE SAME");
  gPbPbInt->Draw("PE SAME");

  TLatex latexRight;
  latexRight.SetNDC();
  latexRight.SetTextSize(0.12);
  latexRight.DrawLatex(0.18, 0.84, "Integrated");

  canvas->cd();
  TLatex latexTop;
  latexTop.SetNDC();
  latexTop.SetTextSize(0.048);
  latexTop.DrawLatex(0.12, 0.95, "CMS Simulation");
  latexTop.SetTextSize(0.044);
  latexTop.DrawLatex(0.35, 0.95, Form("%s J/#psi acceptance", isPrompt ? "Prompt" : "Nonprompt"));
  latexTop.DrawLatex(0.80, 0.95, "pp vs PbPb");

  canvas->SaveAs(outStem + ".pdf");
  canvas->SaveAs(outStem + ".png");

  delete hPpInt;
  delete hPbPbInt;
  delete gPpInt;
  delete gPbPbInt;
}
} // namespace

void draw_Jpsi_acc_pp_pbpb(bool isPrompt = true,
                           const TString &remark = "20260420",
                           int weightOption = 1,
                           int ptWeightSystematic = 0,
                           const TString &ppFilePath = "",
                           const TString &pbpbFilePath = "",
                           const TString &outDirPath = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(8);

  const TString ppPath = ppFilePath.IsNull()
                             ? BuildDefaultInputPath(false, isPrompt, weightOption, ptWeightSystematic, remark)
                             : ResolvePath(ppFilePath);
  const TString pbpbPath = pbpbFilePath.IsNull()
                               ? BuildDefaultInputPath(true, isPrompt, weightOption, ptWeightSystematic, remark)
                               : ResolvePath(pbpbFilePath);

  TFile *fPp = TFile::Open(ppPath, "READ");
  TFile *fPbPb = TFile::Open(pbpbPath, "READ");
  if (!fPp || fPp->IsZombie())
  {
    std::cout << "[ERROR] failed to open pp file: " << ppPath << "\n";
    return;
  }
  if (!fPbPb || fPbPb->IsZombie())
  {
    std::cout << "[ERROR] failed to open PbPb file: " << pbpbPath << "\n";
    return;
  }

  TH1 *hPpMid = LoadHistClone(fPp, "hAccPt_2021_midy", "hPpMid");
  TH1 *hPpFwd = LoadHistClone(fPp, "hAccPt_2021_Fory", "hPpFwd");
  TH1 *hPpMidInt = LoadHistClone(fPp, "hAccPt_2021_midy_Int", "hPpMidInt");
  TH1 *hPpFwdInt = LoadHistClone(fPp, "hAccPt_2021_Fory_Int", "hPpFwdInt");
  TH1 *hPbPbMid = LoadHistClone(fPbPb, "hAccPt_2021_midy", "hPbPbMid");
  TH1 *hPbPbFwd = LoadHistClone(fPbPb, "hAccPt_2021_Fory", "hPbPbFwd");
  TH1 *hPbPbMidInt = LoadHistClone(fPbPb, "hAccPt_2021_midy_Int", "hPbPbMidInt");
  TH1 *hPbPbFwdInt = LoadHistClone(fPbPb, "hAccPt_2021_Fory_Int", "hPbPbFwdInt");
  fPp->Close();
  fPbPb->Close();

  if (!hPpMid || !hPpFwd || !hPpMidInt || !hPpFwdInt ||
      !hPbPbMid || !hPbPbFwd || !hPbPbMidInt || !hPbPbFwdInt)
    return;

  StyleHist(hPpMid, kBlue + 1, 20);
  StyleHist(hPpFwd, kBlue + 1, 20);
  StyleHist(hPpMidInt, kBlue + 1, 20);
  StyleHist(hPpFwdInt, kBlue + 1, 20);
  StyleHist(hPbPbMid, kRed + 1, 21);
  StyleHist(hPbPbFwd, kRed + 1, 21);
  StyleHist(hPbPbMidInt, kRed + 1, 21);
  StyleHist(hPbPbFwdInt, kRed + 1, 21);

  TString outDir = outDirPath.IsNull()
                       ? Form("%s/plot_outputs_acc_compare_jpsi", GetBaseDir().Data())
                       : ResolvePath(outDirPath);
  gSystem->mkdir(outDir, true);

  const TString stateTag = isPrompt ? "prompt" : "nonprompt";
  const TString outStemMid = Form("%s/compare_acc_jpsi_%s_pp_pbpb_mid_%s",
                                  outDir.Data(),
                                  stateTag.Data(),
                                  remark.Data());
  const TString outStemFwd = Form("%s/compare_acc_jpsi_%s_pp_pbpb_fwd_%s",
                                  outDir.Data(),
                                  stateTag.Data(),
                                  remark.Data());

  DrawRegionCanvasWithIntegratedPad(Form("c_Jpsi_acc_pp_pbpb_mid_%s", stateTag.Data()),
                                    "|y| < 1.6",
                                    hPpMid, hPbPbMid,
                                    hPpMidInt, hPbPbMidInt,
                                    6.5, 40.0,
                                    isPrompt,
                                    outStemMid);

  DrawRegionCanvasWithIntegratedPad(Form("c_Jpsi_acc_pp_pbpb_fwd_%s", stateTag.Data()),
                                    "1.6 < |y| < 2.4",
                                    hPpFwd, hPbPbFwd,
                                    hPpFwdInt, hPbPbFwdInt,
                                    3.5, 40.0,
                                    isPrompt,
                                    outStemFwd);

  std::cout << "[INFO] pp input   : " << ppPath << "\n";
  std::cout << "[INFO] PbPb input : " << pbpbPath << "\n";
  std::cout << "[INFO] output(mid): " << outStemMid << ".pdf\n";
  std::cout << "[INFO] output(fwd): " << outStemFwd << ".pdf\n";
}

void draw_Jpsi_acc_pp_pbpb_prompt(const TString &remark = "20260420")
{
  draw_Jpsi_acc_pp_pbpb(true, remark);
}

void draw_Jpsi_acc_pp_pbpb_nonprompt(const TString &remark = "20260420")
{
  draw_Jpsi_acc_pp_pbpb(false, remark);
}
