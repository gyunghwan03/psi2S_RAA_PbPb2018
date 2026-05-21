#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../commonUtility.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

namespace DrawDndptHNP260521
{
constexpr const char *kDir =
    "/data/hwan/psi2S_RAA_PbPb2018/compareDataToMC/ptw_candidates_260520_v3/H_NP_rational";

struct Entry {
  TString file;
};

TString BaseName(const TString &path)
{
  TString base = gSystem->BaseName(path);
  base.ReplaceAll(".root", "");
  return base;
}

TString ComponentLabel(const TString &base)
{
  if (base.Contains("BtoJpsi")) return "Nonprompt J/#psi  (B#rightarrowJ/#psi)";
  return "Prompt J/#psi";
}

TString SystemLabel(const TString &base)
{
  return base.Contains("ratioDataMC_AA_") ? "PbPb  (5.02 TeV)" : "pp  (5.02 TeV)";
}

TString RapidityLabel(const TString &base)
{
  if (base.Contains("_fwd") || base.Contains("y1p6_2p4")) return "1.6 < |y| < 2.4";
  return "|y| < 1.6";
}

void SaveMissing(const TString &base, const TString &path)
{
  TCanvas c(Form("c_missing_%s", base.Data()), "missing dN/dpT input", 600, 700);
  TPad p1("p1", "p1", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.30);
  p1.Draw();
  p2.Draw();

  p1.cd();
  p1.SetTicks(1, 1);
  TPaveText box(0.12, 0.34, 0.88, 0.68, "NDC");
  box.SetBorderSize(0);
  box.SetFillColor(0);
  box.SetTextAlign(22);
  box.SetTextFont(42);
  box.SetTextSize(0.038);
  box.AddText("Missing ROOT input");
  box.AddText(base);
  box.AddText(path);
  box.Draw();

  p2.cd();
  p2.SetTicks(1, 1);
  TPaveText low(0.12, 0.30, 0.88, 0.70, "NDC");
  low.SetBorderSize(0);
  low.SetFillColor(0);
  low.SetTextAlign(22);
  low.SetTextFont(42);
  low.SetTextSize(0.060);
  low.AddText("No Data/MC ratio");
  low.Draw();

  TString out = Form("%s/dndpt_%s", kDir, base.Data());
  c.SaveAs(out + ".pdf");
  c.SaveAs(out + ".png");
}

void DrawOne(const TString &path)
{
  const TString base = BaseName(path);
  if (gSystem->AccessPathName(path)) {
    std::cerr << "[missing] " << path << "\n";
    SaveMissing(base, path);
    return;
  }

  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "[open fail] " << path << "\n";
    SaveMissing(base, path);
    return;
  }

  TH1D *hDataIn = dynamic_cast<TH1D *>(f->Get("hptData_norm_width"));
  TH1D *hMcIn = dynamic_cast<TH1D *>(f->Get("hptMC_norm_width"));
  TH1D *hRatioIn = dynamic_cast<TH1D *>(f->Get("WeightFactor"));
  if (!hDataIn || !hMcIn) {
    std::cerr << "[missing hist] " << path << "\n";
    SaveMissing(base, path);
    return;
  }

  std::unique_ptr<TH1D> hData(dynamic_cast<TH1D *>(hDataIn->Clone("hData_dndpt")));
  std::unique_ptr<TH1D> hMc(dynamic_cast<TH1D *>(hMcIn->Clone("hMc_dndpt")));
  std::unique_ptr<TH1D> hRatio;
  if (hRatioIn) {
    hRatio.reset(dynamic_cast<TH1D *>(hRatioIn->Clone("hRatio_dndpt")));
  } else {
    hRatio.reset(dynamic_cast<TH1D *>(hDataIn->Clone("hRatio_dndpt")));
    hRatio->Divide(hMcIn);
  }
  TF1 *fitIn = dynamic_cast<TF1 *>(f->Get("fitRatio1"));
  if (!fitIn) fitIn = dynamic_cast<TF1 *>(f->Get("dataMC_Ratio1"));

  hData->SetDirectory(nullptr);
  hMc->SetDirectory(nullptr);
  hRatio->SetDirectory(nullptr);

  TCanvas c(Form("c_dndpt_%s", base.Data()), "dN/dpT", 600, 700);
  TPad p1("p1", "p1", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.30);
  p1.Draw();
  p2.Draw();

  p1.cd();
  p1.SetTicks(1, 1);
  p1.SetLeftMargin(0.13);
  p1.SetBottomMargin(0.02);
  p1.SetTopMargin(0.07);
  p1.SetRightMargin(0.04);

  hData->SetTitle("");
  hData->GetXaxis()->SetLabelSize(0);
  hData->GetXaxis()->SetTitle("");
  hData->GetYaxis()->SetTitle("Normalized dN/dp_{T}");
  hData->GetYaxis()->CenterTitle();
  hData->GetYaxis()->SetTitleOffset(1.25);
  hData->GetYaxis()->SetTitleSize(0.045);
  hData->GetYaxis()->SetLabelSize(0.040);
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(1.0);
  hData->SetLineColor(kBlack);
  hData->SetMarkerColor(kBlack);
  hMc->SetLineColor(kRed + 1);
  hMc->SetLineWidth(2);
  hMc->SetMarkerColor(kRed + 1);

  const double ymax = std::max(hData->GetMaximum(), hMc->GetMaximum()) * 1.35;
  hData->SetMinimum(0.0);
  hData->SetMaximum(ymax > 0.0 ? ymax : 1.0);
  hData->Draw("e1");
  hMc->Draw("hist same");
  hData->Draw("e1 same");

  drawText(ComponentLabel(base).Data(), 0.17, 0.87, 1, 18);
  drawText(SystemLabel(base).Data(), 0.17, 0.82, 1, 18);
  drawText(RapidityLabel(base).Data(), 0.17, 0.77, 1, 18);

  TLegend leg(0.62, 0.73, 0.91, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.038);
  leg.AddEntry(hData.get(), "Data", "pe");
  leg.AddEntry(hMc.get(), "MC", "l");
  leg.Draw();

  p2.cd();
  p2.SetTicks(1, 1);
  p2.SetLeftMargin(0.13);
  p2.SetBottomMargin(0.30);
  p2.SetTopMargin(0.02);
  p2.SetRightMargin(0.04);

  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatio->GetYaxis()->SetTitle("Data / MC");
  hRatio->GetXaxis()->CenterTitle();
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetXaxis()->SetTitleSize(0.105);
  hRatio->GetXaxis()->SetLabelSize(0.085);
  hRatio->GetYaxis()->SetTitleSize(0.090);
  hRatio->GetYaxis()->SetLabelSize(0.075);
  hRatio->GetYaxis()->SetTitleOffset(0.60);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.0);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMinimum(0.0);
  hRatio->SetMaximum(4.0);
  hRatio->Draw("e1");
  if (fitIn) {
    fitIn->SetLineColor(kBlue + 1);
    fitIn->SetLineWidth(2);
    fitIn->Draw("same");
  }
  const double xmin = hRatio->GetXaxis()->GetXmin();
  const double xmax = hRatio->GetXaxis()->GetXmax();
  TLine unity(xmin, 1.0, xmax, 1.0);
  unity.SetLineStyle(2);
  unity.Draw("same");

  TString out = Form("%s/dndpt_%s", kDir, base.Data());
  c.SaveAs(out + ".pdf");
  c.SaveAs(out + ".png");
  std::cout << "[wrote] " << out << ".{pdf,png}\n";
}

}  // namespace DrawDndptHNP260521

void draw_dndpt_H_NP_rational_260521()
{
  using namespace DrawDndptHNP260521;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  std::vector<Entry> entries = {
      {"ratioDataMC_AA_BtoJpsi_DATA_fwd.root"},
      {"ratioDataMC_AA_BtoJpsi_DATA_mid.root"},
      {"ratioDataMC_AA_BtoJpsi_DATA_y0_1p6_260515.root"},
      {"ratioDataMC_AA_BtoJpsi_DATA_y1p6_2p4_260515.root"},
      {"ratioDataMC_AA_Jpsi_DATA_fwd.root"},
      {"ratioDataMC_AA_Jpsi_DATA_mid.root"},
      {"ratioDataMC_AA_Jpsi_DATA_y0_1p6_260515.root"},
      {"ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_260515.root"},
      {"ratioDataMC_pp_BtoJpsi_DATA_fwd.root"},
      {"ratioDataMC_pp_BtoJpsi_DATA_mid.root"},
      {"ratioDataMC_pp_BtoJpsi_DATA_noCtau_y0_1p6_260515.root"},
      {"ratioDataMC_pp_Jpsi_DATA_fwd.root"},
      {"ratioDataMC_pp_Jpsi_DATA_mid.root"},
      {"ratioDataMC_pp_Jpsi_DATA_noCtau_y0_1p6_260515.root"},
      {"ratioDataMC_pp_Jpsi_DATA_noCtau_y1p6_2p4_260515.root"},
  };

  for (const auto &entry : entries) {
    DrawOne(Form("%s/%s", kDir, entry.file.Data()));
  }
}
