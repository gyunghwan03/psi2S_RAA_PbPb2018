#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
const int kNSample = 500;

struct CompareSpec
{
  TString tag;
  TString title;
  TString rapLabel;
  double ptMin;
  double ptMax;
  TString currentPath;
  TString oldPath;
};

TString BaseDir()
{
  return gSystem->DirName(__FILE__);
}

TString OutDir()
{
  return Form("%s/plot_outputs_ptweight_compare_20260411", BaseDir().Data());
}

TF1 *LoadPtWeight(const TString &path, const TString &cloneName)
{
  TFile *file = TFile::Open(path, "READ");
  if (!file || file->IsZombie())
  {
    std::cout << "[ERROR] failed to open pT weight file: " << path << "\n";
    if (file)
      file->Close();
    return nullptr;
  }

  TF1 *func = dynamic_cast<TF1 *>(file->Get("dataMC_Ratio1"));
  if (!func)
    func = dynamic_cast<TF1 *>(file->Get("fitRatio1"));

  if (!func)
  {
    std::cout << "[ERROR] missing dataMC_Ratio1/fitRatio1 in: " << path << "\n";
    file->Close();
    return nullptr;
  }

  TF1 *clone = static_cast<TF1 *>(func->Clone(cloneName));
  file->Close();
  return clone;
}

TGraph *MakeFunctionGraph(TF1 *func, const CompareSpec &spec, const TString &name)
{
  if (!func)
    return nullptr;

  TGraph *graph = new TGraph(kNSample);
  graph->SetName(name);
  for (int i = 0; i < kNSample; ++i)
  {
    const double x = spec.ptMin + (spec.ptMax - spec.ptMin) * i / (kNSample - 1);
    graph->SetPoint(i, x, func->Eval(x));
  }
  return graph;
}

TGraph *MakeRatioGraph(TF1 *num, TF1 *den, const CompareSpec &spec, const TString &name)
{
  if (!num || !den)
    return nullptr;

  TGraph *graph = new TGraph(kNSample);
  graph->SetName(name);
  for (int i = 0; i < kNSample; ++i)
  {
    const double x = spec.ptMin + (spec.ptMax - spec.ptMin) * i / (kNSample - 1);
    const double denVal = den->Eval(x);
    const double ratio = (denVal != 0.0) ? num->Eval(x) / denVal : 0.0;
    graph->SetPoint(i, x, ratio);
  }
  return graph;
}

void StyleGraph(TGraph *graph, int color, int lineStyle)
{
  if (!graph)
    return;
  graph->SetLineColor(color);
  graph->SetLineStyle(lineStyle);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(color);
}

std::pair<double, double> GraphRange(const std::vector<TGraph *> &graphs)
{
  double minVal = std::numeric_limits<double>::max();
  double maxVal = -std::numeric_limits<double>::max();
  for (TGraph *graph : graphs)
  {
    if (!graph)
      continue;
    for (int i = 0; i < graph->GetN(); ++i)
    {
      double x = 0.0;
      double y = 0.0;
      graph->GetPoint(i, x, y);
      if (!std::isfinite(y))
        continue;
      minVal = std::min(minVal, y);
      maxVal = std::max(maxVal, y);
    }
  }

  if (minVal == std::numeric_limits<double>::max())
    return {0.0, 1.0};
  return {minVal, maxVal};
}

void SetPaddedRange(TH1D *frame, double minVal, double maxVal, bool includeOne)
{
  if (includeOne)
  {
    minVal = std::min(minVal, 1.0);
    maxVal = std::max(maxVal, 1.0);
  }

  const double span = std::max(0.03, maxVal - minVal);
  frame->GetYaxis()->SetRangeUser(std::max(0.0, minVal - 0.18 * span), maxVal + 0.22 * span);
}

void DrawLatexLabel(const TString &title, const TString &rapLabel)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.043);
  latex.DrawLatex(0.16, 0.86, title);
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.16, 0.805, Form("PbPb efficiency p_{T} weight, %s", rapLabel.Data()));
}

void DrawOne(const CompareSpec &spec)
{
  TF1 *current = LoadPtWeight(spec.currentPath, spec.tag + "_current");
  TF1 *old = LoadPtWeight(spec.oldPath, spec.tag + "_old");
  if (!current || !old)
  {
    delete current;
    delete old;
    return;
  }

  TGraph *gCurrent = MakeFunctionGraph(current, spec, spec.tag + "_gCurrent");
  TGraph *gOld = MakeFunctionGraph(old, spec, spec.tag + "_gOld");
  TGraph *gRatio = MakeRatioGraph(current, old, spec, spec.tag + "_gRatio");
  StyleGraph(gCurrent, kRed + 1, 1);
  StyleGraph(gOld, kBlue + 1, 2);
  StyleGraph(gRatio, kBlack, 1);

  TCanvas *canvas = new TCanvas("c_" + spec.tag, "", 820, 820);
  TPad *top = new TPad(spec.tag + "_top", "", 0.0, 0.30, 1.0, 1.0);
  TPad *bot = new TPad(spec.tag + "_bot", "", 0.0, 0.0, 1.0, 0.30);
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
  TH1D *frame = new TH1D(spec.tag + "_frame", "", 1, spec.ptMin, spec.ptMax);
  frame->SetDirectory(nullptr);
  frame->GetYaxis()->SetTitle("p_{T} weight");
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.05);
  const auto weightRange = GraphRange({gCurrent, gOld});
  SetPaddedRange(frame, weightRange.first, weightRange.second, true);
  frame->Draw();
  gOld->Draw("L SAME");
  gCurrent->Draw("L SAME");
  TLine *oneTop = new TLine(spec.ptMin, 1.0, spec.ptMax, 1.0);
  oneTop->SetLineStyle(7);
  oneTop->SetLineColor(kGray + 2);
  oneTop->Draw("SAME");

  TLegend *leg = new TLegend(0.53, 0.68, 0.91, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.034);
  leg->AddEntry(gCurrent, "Current eff (260410)", "l");
  leg->AddEntry(gOld, "Old eff (260310)", "l");
  leg->Draw();
  DrawLatexLabel(spec.title, spec.rapLabel);

  bot->cd();
  TH1D *ratioFrame = new TH1D(spec.tag + "_ratioFrame", "", 1, spec.ptMin, spec.ptMax);
  ratioFrame->SetDirectory(nullptr);
  ratioFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ratioFrame->GetYaxis()->SetTitle("Current/Old");
  ratioFrame->GetXaxis()->SetTitleSize(0.12);
  ratioFrame->GetXaxis()->SetLabelSize(0.10);
  ratioFrame->GetYaxis()->SetTitleSize(0.10);
  ratioFrame->GetYaxis()->SetLabelSize(0.09);
  ratioFrame->GetYaxis()->SetTitleOffset(0.55);
  ratioFrame->GetYaxis()->SetNdivisions(505);
  const auto ratioRange = GraphRange({gRatio});
  SetPaddedRange(ratioFrame, ratioRange.first, ratioRange.second, true);
  ratioFrame->Draw();
  gRatio->Draw("L SAME");
  TLine *oneBot = new TLine(spec.ptMin, 1.0, spec.ptMax, 1.0);
  oneBot->SetLineStyle(2);
  oneBot->Draw("SAME");

  const TString outDir = OutDir();
  gSystem->mkdir(outDir, true);
  canvas->SaveAs(Form("%s/%s.pdf", outDir.Data(), spec.tag.Data()));

  delete current;
  delete old;
}
} // namespace

void plot_compare_ptweight_jpsi_pbpb_eff_20260411()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const std::vector<CompareSpec> specs = {
      {"ptweight_jpsi_pbpb_mid_eff_prompt_current_vs_20260411",
       "Prompt J/#psi",
       "|y| < 1.6",
       6.5,
       40.0,
       "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260410.root",
       "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310.root"},
      {"ptweight_jpsi_pbpb_mid_eff_nonprompt_current_vs_20260411",
       "Nonprompt J/#psi",
       "|y| < 1.6",
       6.5,
       40.0,
       "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260410.root",
       "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310.root"},
      {"ptweight_jpsi_pbpb_fwd_eff_prompt_current_vs_20260411",
       "Prompt J/#psi",
       "1.6 < |y| < 2.4",
       3.5,
       40.0,
       "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_260410.root",
       "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310.root"},
      {"ptweight_jpsi_pbpb_fwd_eff_nonprompt_current_vs_20260411",
       "Nonprompt J/#psi",
       "1.6 < |y| < 2.4",
       3.5,
       40.0,
       "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_260410.root",
       "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310.root"}};

  for (const auto &spec : specs)
    DrawOne(spec);
}
