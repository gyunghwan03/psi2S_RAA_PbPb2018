#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <vector>

using namespace std;

namespace {

TString firstExisting(const vector<TString>& candidates) {
  for (const auto& p : candidates) {
    if (!gSystem->AccessPathName(p)) return p;
  }
  return "";
}

TF1* getRatioFunction(const TString& fileName, const TString& newName) {
  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    cout << "[ERROR] Cannot open file: " << fileName << endl;
    return nullptr;
  }

  TF1* fin = dynamic_cast<TF1*>(f->Get("dataMC_Ratio1"));
  if (!fin) {
    cout << "[ERROR] dataMC_Ratio1 not found in: " << fileName << endl;
    f->Close();
    return nullptr;
  }

  TF1* fout = dynamic_cast<TF1*>(fin->Clone(newName));
  if (!fout) {
    cout << "[ERROR] Failed to clone TF1 from: " << fileName << endl;
    f->Close();
    return nullptr;
  }
  f->Close();
  return fout;
}

void sampleMinMax(TF1* fA, TF1* fB, double xMin, double xMax, double& yMin, double& yMax) {
  yMin = 1e9;
  yMax = -1e9;

  const int n = 300;
  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double yA = fA->Eval(x);
    const double yB = fB->Eval(x);
    if (TMath::Finite(yA)) { if (yA < yMin) yMin = yA; if (yA > yMax) yMax = yA; }
    if (TMath::Finite(yB)) { if (yB < yMin) yMin = yB; if (yB > yMax) yMax = yB; }
  }

  if (yMin >= yMax) {
    yMin = 0.0;
    yMax = 2.0;
  } else {
    const double dy = yMax - yMin;
    yMin -= 0.15 * dy;
    yMax += 0.20 * dy;
  }
}

TGraph* buildDiffGraph(TF1* fA, TF1* fB, double xMin, double xMax) {
  vector<double> vx;
  vector<double> vy;
  const int n = 300;
  vx.reserve(n);
  vy.reserve(n);

  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double yA = fA->Eval(x);
    const double yB = fB->Eval(x);
    if (!TMath::Finite(yA) || !TMath::Finite(yB)) continue;
    vx.push_back(x);
    vy.push_back(yA - yB);
  }

  if (vx.empty()) return nullptr;
  TGraph* g = new TGraph(static_cast<int>(vx.size()), &vx[0], &vy[0]);
  g->SetName("g_diff");
  return g;
}

void sampleDiffAbsMax(TF1* fA, TF1* fB, double xMin, double xMax, double& yAbsMax) {
  yAbsMax = 0.0;
  const int n = 300;
  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double yA = fA->Eval(x);
    const double yB = fB->Eval(x);
    if (!TMath::Finite(yA) || !TMath::Finite(yB)) continue;
    const double d = TMath::Abs(yA - yB);
    if (d > yAbsMax) yAbsMax = d;
  }
  if (yAbsMax <= 0.) yAbsMax = 0.1;
}

void drawCompare(TF1* fA, TF1* fB,
                 const TString& labelA, const TString& labelB,
                 const TString& titleText,
                 double xMin, double xMax,
                 const TString& outPdf, const TString& outPng) {
  static int idx = 0;
  ++idx;

  double yMin = 0.0, yMax = 2.0;
  sampleMinMax(fA, fB, xMin, xMax, yMin, yMax);
  double yAbsDiff = 0.1;
  sampleDiffAbsMax(fA, fB, xMin, xMax, yAbsDiff);

  TCanvas* c = new TCanvas(Form("c_cmp_%d", idx), "", 850, 800);
  TPad* pTop = new TPad(Form("pTop_cmp_%d", idx), "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad(Form("pBot_cmp_%d", idx), "", 0.0, 0.0, 1.0, 0.30);

  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.04);
  pTop->SetTicks(1, 1);

  pBot->SetTopMargin(0.03);
  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.04);
  pBot->SetTicks(1, 1);

  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  TH1D* hFrame = new TH1D(Form("hFrame_%d", idx), ";p_{T} (GeV/c);Data/MC fit", 100, xMin, xMax);
  hFrame->SetMinimum(yMin);
  hFrame->SetMaximum(yMax);
  hFrame->GetXaxis()->SetLabelSize(0.0);
  hFrame->Draw("axis");

  fA->SetLineWidth(3);
  fA->SetLineStyle(1);
  fB->SetLineWidth(3);
  fB->SetLineStyle(2);
  fA->Draw("same");
  fB->Draw("same");

  TLegend* leg = new TLegend(0.52, 0.72, 0.90, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fA, labelA, "l");
  leg->AddEntry(fB, labelB, "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.14, 0.92, titleText);

  pBot->cd();
  TH1D* hDiffFrame = new TH1D(Form("hDiffFrame_%d", idx), ";p_{T} (GeV/c);A-B", 100, xMin, xMax);
  hDiffFrame->SetMinimum(-1.3 * yAbsDiff);
  hDiffFrame->SetMaximum( 1.3 * yAbsDiff);
  hDiffFrame->GetXaxis()->SetTitleSize(0.12);
  hDiffFrame->GetXaxis()->SetLabelSize(0.10);
  hDiffFrame->GetYaxis()->SetTitleSize(0.10);
  hDiffFrame->GetYaxis()->SetLabelSize(0.09);
  hDiffFrame->GetYaxis()->SetTitleOffset(0.50);
  hDiffFrame->GetYaxis()->SetNdivisions(505);
  hDiffFrame->Draw("axis");

  TGraph* gDiff = buildDiffGraph(fA, fB, xMin, xMax);
  if (gDiff) {
    gDiff->SetLineColor(kBlack);
    gDiff->SetLineWidth(2);
    gDiff->Draw("L same");
  }

  TLine* l0 = new TLine(xMin, 0.0, xMax, 0.0);
  l0->SetLineColor(kGray + 2);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  l0->Draw("same");

  gSystem->mkdir("./figs", kTRUE);
  c->SaveAs(outPdf);
  c->SaveAs(outPng);
}

} // namespace

// PR = 0: Prompt, PR = 1: NonPrompt(B->psi2S)
// system = "pp" or "AA" (PbPb)
void compare_psi2S_260323_2exp_vs_260310_2exp(int PR, TString system) {
  gStyle->SetOptStat(0);

  if (system != "pp" && system != "AA") {
    cout << "[ERROR] system must be \"pp\" or \"AA\"" << endl;
    return;
  }

  TString tag = (PR == 0) ? "psi2S" : "Btopsi2S";
  TString prLabel = (PR == 0) ? "Prompt" : "NonPrompt";
  TString sysLabel = (system == "AA") ? "PbPb" : "pp";

  vector<TString> tags;
  if (PR == 0) {
    tags = {"psi2S"};
  } else {
    tags = {"Btopsi2S", "BtoPsi2S"};
  }

  vector<TString> c310y0, c323y0, c323y1;
  for (const auto& t : tags) {
    c310y0.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y0_2p4_260310_2exp.root", system.Data(), t.Data()));
    c310y0.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y0_2p4_260310.root", system.Data(), t.Data()));

    c323y0.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y0_2p4_260323_2exp.root", system.Data(), t.Data()));
    c323y0.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y0_2p4_260323.root", system.Data(), t.Data()));

    c323y1.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y1p6_2p4_260323_2exp.root", system.Data(), t.Data()));
    c323y1.push_back(Form("./ratioDataMC_%s_%s_DATA_ctauCut_y1p6_2p4_260323.root", system.Data(), t.Data()));
  }

  TString f310_y0 = firstExisting(c310y0);
  TString f323_y0 = firstExisting(c323y0);
  TString f323_y1 = firstExisting(c323y1);

  if (f310_y0.IsNull()) {
    cout << "[ERROR] y0_2p4 260310(2exp) input not found for " << system << ", tag: " << tag << endl;
    return;
  }
  if (f323_y0.IsNull()) {
    cout << "[ERROR] y0_2p4 260323(2exp) input not found for " << system << ", tag: " << tag << endl;
    return;
  }
  if (f323_y1.IsNull()) {
    cout << "[ERROR] y1p6_2p4 260323(2exp) input not found for " << system << ", tag: " << tag << endl;
    return;
  }

  TF1* f260310_y0 = getRatioFunction(f310_y0, "f260310_y0");
  TF1* f260323_y0 = getRatioFunction(f323_y0, "f260323_y0");
  TF1* f260323_y1 = getRatioFunction(f323_y1, "f260323_y1");
  if (!f260310_y0 || !f260323_y0 || !f260323_y1) return;

  f260323_y1->SetLineColor(kRed + 1);
  f260310_y0->SetLineColor(kBlue + 1);
  f260323_y0->SetLineColor(kGreen + 2);

  // 1) y1p6_2p4(260323_2exp) vs y0_2p4(260310_2exp), pt: 3-12
  drawCompare(
    f260323_y1, f260310_y0,
    Form("260323 y1.6-2.4 (%s)", prLabel.Data()),
    Form("260310 y0-2.4 (%s)", prLabel.Data()),
    Form("psi(2S) %s fit comparison (%s): 260323 y1.6-2.4 vs 260310 y0-2.4", sysLabel.Data(), prLabel.Data()),
    3.0, 12.0,
    Form("./figs/compare_psi2S_fit_%s_%s_260323y1p6_2p4_vs_260310y0_2p4_pt3_12.pdf", system.Data(), tag.Data()),
    Form("./figs/compare_psi2S_fit_%s_%s_260323y1p6_2p4_vs_260310y0_2p4_pt3_12.png", system.Data(), tag.Data())
  );

  // 2) y0_2p4(260323_2exp) vs y0_2p4(260310_2exp), pt: 6.5-40
  drawCompare(
    f260323_y0, f260310_y0,
    Form("260323 y0-2.4 (%s)", prLabel.Data()),
    Form("260310 y0-2.4 (%s)", prLabel.Data()),
    Form("psi(2S) %s fit comparison (%s): 260323 vs 260310 at y0-2.4", sysLabel.Data(), prLabel.Data()),
    6.5, 40.0,
    Form("./figs/compare_psi2S_fit_%s_%s_260323y0_2p4_vs_260310y0_2p4_pt6p5_40.pdf", system.Data(), tag.Data()),
    Form("./figs/compare_psi2S_fit_%s_%s_260323y0_2p4_vs_260310y0_2p4_pt6p5_40.png", system.Data(), tag.Data())
  );

  cout << "[DONE] Saved TF1 comparison plots in ./figs/ for system=" << system << endl;
}

// Backward-compatible default: pp
void compare_psi2S_260323_2exp_vs_260310_2exp(int PR = 0) {
  compare_psi2S_260323_2exp_vs_260310_2exp(PR, "pp");
}

void compare_psi2S_260323_2exp_vs_260310_2exp_allSystems(int PR = 0) {
  compare_psi2S_260323_2exp_vs_260310_2exp(PR, "pp");
  compare_psi2S_260323_2exp_vs_260310_2exp(PR, "AA");
}
