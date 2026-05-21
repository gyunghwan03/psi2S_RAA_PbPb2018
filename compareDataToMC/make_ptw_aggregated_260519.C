// =============================================================
//  make_ptw_aggregated_260519.C
//
//  Option 2 (rebin): high-pT data points are noisy and pull the
//  2-exp fit to negative values around pT ~ 25-30 GeV (especially
//  for PbPb). To stabilize the fit we read the same native-binned
//  yields used by B_nominal_2exp, but BEFORE fitting we aggregate
//  the last two bins ({20-25, 25-40} mid; {15-20, 20-40} fwd) into
//  one effective high-pT data point. The fit then sees 5 (mid) or
//  7 (fwd) points and stays positive throughout 6.5-40 GeV.
//
//  Output: compareDataToMC/ptw_candidates_260519_v2/G_aggregateHighPt_2exp/
// =============================================================

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../commonUtility.h"
#include "../cutsAndBin.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace PtWAgg260519
{

constexpr const char *kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
constexpr const char *kOutBase = "ptw_candidates_260519_v2";
constexpr const char *kCandTag = "G_aggregateHighPt_2exp";

struct YieldErr { double val = 0.0, err = 0.0; bool ok = false; };

YieldErr GetPbPbYield(double pl, double ph, double yl, double yh)
{
  YieldErr r;
  TString lab = getKineLabel(pl, ph, yl, yh, 0.0, 0, 180);
  TString p = Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Mass/"
                   "Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                   kRepo, lab.Data());
  std::unique_ptr<TFile> f(TFile::Open(p, "READ"));
  if (!f || f->IsZombie()) return r;
  TH1D *h = dynamic_cast<TH1D *>(f->Get("fitResults"));
  if (!h) return r;
  r.val = h->GetBinContent(1);
  r.err = h->GetBinError(1);
  r.ok = r.val > 0.0;
  return r;
}

double GetPbPbFrac(double pl, double ph, double yl, double yh)
{
  TString lab = getKineLabel(pl, ph, yl, yh, 0.0, 0, 180);
  TString p = Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Final/"
                   "2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                   kRepo, lab.Data());
  std::unique_ptr<TFile> f(TFile::Open(p, "READ"));
  if (!f || f->IsZombie()) return 0.0;
  TH1D *h = dynamic_cast<TH1D *>(f->Get("2DfitResults"));
  return h ? h->GetBinContent(1) : 0.0;
}

YieldErr GetPpYield(bool prompt, double pl, double ph, double yl, double yh)
{
  YieldErr r;
  TString lab = getKineLabelpp(pl, ph, yl, yh, 0.0);
  TString massPath = Form("%s/Macros/pp_Jpsi/roots/2DFit_No_Weight/Mass/"
                          "Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                          kRepo, lab.Data());
  TString fracPath = Form("%s/Macros/pp_Jpsi/roots/2DFit_No_Weight/Final/"
                          "2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                          kRepo, lab.Data());
  if (gSystem->AccessPathName(massPath) || gSystem->AccessPathName(fracPath)) {
    std::cerr << "[GetPpYield] missing pp_Jpsi Mass/Final input for " << lab
              << " (legacy fallback is disabled)" << std::endl;
    return r;
  }

  std::unique_ptr<TFile> fMass(TFile::Open(massPath, "READ"));
  std::unique_ptr<TFile> fFrac(TFile::Open(fracPath, "READ"));
  if (!fMass || fMass->IsZombie() || !fFrac || fFrac->IsZombie()) return r;
  TH1D *hMass = dynamic_cast<TH1D *>(fMass->Get("fitResults"));
  TH1D *hFrac = dynamic_cast<TH1D *>(fFrac->Get("2DfitResults"));
  if (!hMass || !hFrac) return r;
  const double y = hMass->GetBinContent(1);
  const double yErr = hMass->GetBinError(1);
  const double frac = hFrac->GetBinContent(1);
  const double fracErr = hFrac->GetBinError(1);
  const double compFrac = prompt ? (1.0 - frac) : frac;
  const double compFracErr = fracErr;
  if (y <= 0.0 || compFrac <= 0.0) return r;
  r.val = y * compFrac;
  r.err = r.val * std::sqrt(std::pow(yErr / y, 2) + std::pow(compFracErr / compFrac, 2));
  r.ok = r.val > 0.0;
  return r;
}

struct McSample
{
  TString chainPath;
  double massLow, massHigh, yLow, yHigh;
  std::unique_ptr<TH1D> hFine;
  bool loaded = false;
};

void LoadMc(McSample &s)
{
  if (s.loaded) return;
  s.loaded = true;
  s.hFine.reset(new TH1D(Form("hFine_%p", (void *)&s), "; p_{T} (GeV/c); ", 2000, 0.0, 50.0));
  s.hFine->SetDirectory(nullptr);

  TChain tree("mmepevt");
  tree.Add(s.chainPath.Data());
  if (tree.GetEntries() <= 0) { std::cerr << "[mc] no entries " << s.chainPath << "\n"; return; }

  constexpr int N = 1000;
  float mass[N], pt[N], y[N];
  Int_t nD = 0;
  tree.SetBranchAddress("nDimu", &nD);
  tree.SetBranchAddress("mass", mass);
  tree.SetBranchAddress("y", y);
  tree.SetBranchAddress("pt", pt);
  Long64_t n = tree.GetEntries();
  std::cout << "[mc] " << s.chainPath << " entries=" << n << "\n";
  for (Long64_t i = 0; i < n; ++i) {
    tree.GetEntry(i);
    for (int j = 0; j < nD; ++j) {
      double ay = std::fabs(y[j]);
      if (mass[j] <= s.massLow || mass[j] >= s.massHigh) continue;
      if (pt[j] <= 0 || pt[j] >= 50) continue;
      if (ay <= s.yLow || ay >= s.yHigh) continue;
      s.hFine->Fill(pt[j]);
    }
  }
}

void ProcessComp(bool isPbPb, bool prompt, bool fwd, McSample &mc, const TString &outDir)
{
  LoadMc(mc);
  if (!mc.hFine || mc.hFine->Integral() <= 0) return;

  const std::vector<double> midNative = {6.5, 9, 12, 15, 20, 25, 40};
  const std::vector<double> fwdNative = {3.0, 4.0, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 40.0};
  const std::vector<double> &nat = fwd ? fwdNative : midNative;
  const double yLow = fwd ? 1.6 : 0.0;
  const double yHigh = fwd ? 2.4 : 1.6;
  const TString rapLab = fwd ? "y1p6_2p4" : "y0_1p6";
  const TString sysLab = isPbPb ? "AA" : "pp";
  const TString compLab = prompt ? "Jpsi" : "BtoJpsi";

  // ----- aggregated edges: merge last 2 native bins into one -----
  std::vector<double> agg(nat.begin(), nat.end() - 2);
  agg.push_back(nat[nat.size() - 1]);  // outer edge of last native bin
  const int nAgg = static_cast<int>(agg.size()) - 1;

  // ----- read native-bin data, then sum into aggregated bins -----
  // The first (nAgg - 1) aggregated bins map 1:1 to native bins; the LAST
  // aggregated bin is the sum of the last two native bins.
  std::vector<double> aggY(nAgg, 0.0), aggE(nAgg, 0.0);
  for (int i = 0; i < static_cast<int>(nat.size()) - 1; ++i) {
    YieldErr ye;
    double w = 1.0;
    if (isPbPb) {
      ye = GetPbPbYield(nat[i], nat[i + 1], yLow, yHigh);
      const double frac = GetPbPbFrac(nat[i], nat[i + 1], yLow, yHigh);
      w = prompt ? (1.0 - frac) : frac;
    } else {
      ye = GetPpYield(prompt, nat[i], nat[i + 1], yLow, yHigh);
    }
    if (!ye.ok) continue;
    const int aggBin = std::min(i, nAgg - 1);
    aggY[aggBin] += ye.val * w;
    aggE[aggBin] = std::sqrt(aggE[aggBin] * aggE[aggBin] + std::pow(ye.err * std::max(w, 1e-3), 2));
  }

  // ----- build histograms at aggregated binning -----
  std::unique_ptr<TH1D> hData(new TH1D(Form("hData_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                       "; p_{T} (GeV/c); ", nAgg, agg.data()));
  std::unique_ptr<TH1D> hDataR(new TH1D(Form("hDataR_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                        "; p_{T} (GeV/c); ", nAgg, agg.data()));
  hData->SetDirectory(nullptr);
  hDataR->SetDirectory(nullptr);
  hData->Sumw2();
  hDataR->Sumw2();
  for (int b = 1; b <= nAgg; ++b) {
    hData->SetBinContent(b, aggY[b - 1]);
    hData->SetBinError(b, aggE[b - 1]);
    hDataR->SetBinContent(b, aggY[b - 1]);
    hDataR->SetBinError(b, aggE[b - 1]);
  }
  if (hData->Integral() <= 0) { std::cerr << "[skip] zero data " << sysLab << compLab << rapLab << "\n"; return; }

  std::unique_ptr<TH1D> hMc(new TH1D(Form("hMC_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                     "; p_{T} (GeV/c); ", nAgg, agg.data()));
  std::unique_ptr<TH1D> hMcR(new TH1D(Form("hMCR_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                      "; p_{T} (GeV/c); ", nAgg, agg.data()));
  hMc->SetDirectory(nullptr);
  hMcR->SetDirectory(nullptr);
  hMc->Sumw2();
  hMcR->Sumw2();
  // rebin fine MC into aggregated bins
  for (int ix = 1; ix <= mc.hFine->GetNbinsX(); ++ix) {
    double x = mc.hFine->GetXaxis()->GetBinCenter(ix);
    double w = mc.hFine->GetBinContent(ix);
    if (w <= 0) continue;
    if (x < agg.front() || x >= agg.back()) continue;
    hMc->Fill(x, w);
    hMcR->Fill(x, w);
  }
  if (hMc->Integral() <= 0) { std::cerr << "[skip] zero MC\n"; return; }

  hData->Scale(1.0 / hData->Integral());
  hDataR->Scale(1.0 / hDataR->Integral());
  hMc->Scale(1.0 / hMc->Integral());
  hMcR->Scale(1.0 / hMcR->Integral());
  TH1ScaleByWidth(hData.get());
  TH1ScaleByWidth(hMc.get());
  hDataR->Divide(hMcR.get());

  // 2-exp fit (same form as B_nominal_2exp) on aggregated points
  std::unique_ptr<TF1> fit(new TF1("fitRatio1",
                                   "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]",
                                   fwd ? 3.0 : 3.0, 40.0));
  fit->SetNpx(1000);
  fit->SetParameters(0.5, 0.5, 0.5, 0.05, 0.5);
  // Constrain to keep PtW positive: e (offset) must be >= 0, and shape
  // coefficients must be >= 0 so the function cannot dive through zero.
  fit->SetParLimits(0, 0.0, 1e4);
  fit->SetParLimits(1, 1e-4, 5.0);
  fit->SetParLimits(2, 0.0, 1e4);
  fit->SetParLimits(3, 1e-4, 5.0);
  fit->SetParLimits(4, 0.0, 10.0);
  hDataR->Fit(fit.get(), "IEQ", "", 3.0, 40.0);
  TFitResultPtr res = hDataR->Fit(fit.get(), "SQ", "", 3.0, 40.0);
  if (res.Get())
    std::cout << "[fit] " << sysLab << " " << compLab << " " << rapLab
              << " chi2/ndf = " << res->Chi2() << " / " << res->Ndf() << "\n";

  // ----- write output ROOT in 260515 naming convention -----
  TString outName;
  if (isPbPb)
    outName = Form("%s/ratioDataMC_AA_%s_DATA_%s_260515.root",
                   outDir.Data(), compLab.Data(), rapLab.Data());
  else
    outName = Form("%s/ratioDataMC_pp_%s_DATA_noCtau_%s_260515.root",
                   outDir.Data(), compLab.Data(), rapLab.Data());
  std::unique_ptr<TFile> out(TFile::Open(outName, "RECREATE"));
  if (!out || out->IsZombie()) { std::cerr << "[write fail] " << outName << "\n"; return; }
  out->cd();
  hDataR->Write("WeightFactor", TObject::kOverwrite);
  fit->Write("dataMC_Ratio1", TObject::kOverwrite);
  fit->Write("fitRatio1", TObject::kOverwrite);
  hData->Write("hptData_norm_width");
  hMc->Write("hptMC_norm_width");

  TCanvas cnv("c_ptw", kCandTag, 600, 700);
  TPad p1("p1", "p1", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.30);
  p1.Draw(); p2.Draw();
  p1.cd();
  hData->SetMarkerStyle(20); hData->SetLineColor(kBlack); hData->SetMarkerColor(kBlack);
  hMc->SetLineColor(kRed + 1); hMc->SetLineWidth(2);
  double ymax = std::max(hData->GetMaximum(), hMc->GetMaximum()) * 1.3;
  hData->SetMaximum(ymax > 0 ? ymax : 1.0);
  hData->Draw("e1");
  hMc->Draw("hist same");
  // Identification labels
  {
    const TString cLbl = prompt ? "Prompt J/#psi" : "Nonprompt J/#psi  (B#rightarrowJ/#psi)";
    const TString sLbl = isPbPb ? "PbPb  (5.02 TeV)" : "pp  (5.02 TeV)";
    const TString rLbl = fwd ? "1.6 < |y| < 2.4" : "|y| < 1.6";
    drawText(cLbl.Data(), 0.17, 0.87, 1, 18);
    drawText(sLbl.Data(), 0.17, 0.82, 1, 18);
    drawText(rLbl.Data(), 0.17, 0.77, 1, 18);
  }
  p2.cd();
  hDataR->SetMarkerStyle(20); hDataR->SetLineColor(kBlack); hDataR->SetMarkerColor(kBlack);
  hDataR->SetMinimum(0); hDataR->SetMaximum(4);
  hDataR->Draw("e1");
  fit->SetLineColor(kBlue + 1); fit->SetLineWidth(2); fit->Draw("same");
  TLine unity(3, 1, 40, 1); unity.SetLineStyle(2); unity.Draw("same");
  cnv.Write("canvas_A");

  TNamed note("candidate_label", "G_aggregateHighPt_2exp (option 2): last 2 native pT bins merged into one for fit");
  note.Write();
  TString pdfP = outName; pdfP.ReplaceAll(".root", ".pdf");
  cnv.SaveAs(pdfP);
  TString pngP = outName; pngP.ReplaceAll(".root", ".png");
  cnv.SaveAs(pngP);
  out->Close();
  std::cout << "[wrote] " << outName << "\n";
}

}  // namespace

void make_ptw_aggregated_260519()
{
  using namespace PtWAgg260519;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  McSample mcPpPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcPpNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_miniAOD_251103.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcAaPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcAaNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_NonPrompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};

  TString outDir = Form("%s/compareDataToMC/%s/%s", kRepo, kOutBase, kCandTag);
  if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir))
  {
    std::cerr << "[err] cannot mkdir " << outDir << "\n";
    return;
  }

  std::cout << "\n##### G_aggregateHighPt_2exp (option 2: high-pT bin aggregation) #####\n";
  ProcessComp(false, true, false, mcPpPr, outDir);
  ProcessComp(false, true, true, mcPpPr, outDir);
  ProcessComp(false, false, false, mcPpNp, outDir);
  ProcessComp(false, false, true, mcPpNp, outDir);
  ProcessComp(true, true, false, mcAaPr, outDir);
  ProcessComp(true, true, true, mcAaPr, outDir);
  ProcessComp(true, false, false, mcAaNp, outDir);
  ProcessComp(true, false, true, mcAaNp, outDir);
}
