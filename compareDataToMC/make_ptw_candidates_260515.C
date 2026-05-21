// =============================================================
//  make_ptw_candidates_260515.C
//
//  Generates 6 pT reweighting candidates for the J/psi PbPb 2018
//  RAA closure study against HIN-16-025.
//
//  For each candidate we produce 8 ROOT files:
//      ratioDataMC_AA_Jpsi_DATA_y0_1p6_260515.root      (PR PbPb mid)
//      ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_260515.root    (PR PbPb fwd)
//      ratioDataMC_AA_BtoJpsi_DATA_y0_1p6_260515.root   (NP PbPb mid)
//      ratioDataMC_AA_BtoJpsi_DATA_y1p6_2p4_260515.root (NP PbPb fwd)
//      ratioDataMC_pp_Jpsi_DATA_noCtau_y0_1p6_260515.root   (PR pp  mid)
//      ratioDataMC_pp_Jpsi_DATA_noCtau_y1p6_2p4_260515.root (PR pp  fwd)
//      ratioDataMC_pp_BtoJpsi_DATA_noCtau_y0_1p6_260515.root   (NP pp  mid)
//      ratioDataMC_pp_BtoJpsi_DATA_noCtau_y1p6_2p4_260515.root (NP pp  fwd)
//
//  Each ROOT exposes the standard "WeightFactor" TH1D and
//  "dataMC_Ratio1" / "fitRatio1" TF1 picked up by the
//  Eff_Acc_260515/efficiency_1d*.C and acceptance_1d.C producers.
//
//  No ctau cut is applied; b-fraction splits prompt/nonprompt from
//  the same yield extraction used by dndpt_Jpsi_pp_noCtau_260513.C
//  and dndpt_y0_1p6_Jpsi.C (the legacy ctau-free PbPb script).
//
//  Candidates vary along:
//     binning (mid coarse / mid nominal / mid HIN-like)
//     fwd binning (nominal vs coarse)
//     fit function (rational, double-exp, modified-Tsallis-like
//                   power, simple power).
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
#include <string>
#include <vector>

namespace PtWCandidates260515
{

constexpr const char *kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
constexpr const char *kOutBase = "ptw_candidates_260515";

// ---- candidate definition ------------------------------------------------

struct Candidate
{
  TString tag;            // dir name under ptw_candidates_260515/
  TString label;          // free-form label
  std::vector<double> ptMid;
  std::vector<double> ptFwd;
  // ROOT TFormula string. {a,b,c,d,e} are par[0..4].
  TString fitFormula;
  double fitMin;
  double fitMax;
  // optional pre-seed parameters (size <= 5). If empty, ROOT defaults.
  std::vector<double> seedPars;
};

std::vector<Candidate> BuildCandidates()
{
  std::vector<Candidate> cs;

  // A: nominal binning + rational
  cs.push_back(Candidate{
      "A_nominal_rational",
      "nominal 13/8-bin, rational (a+bx+cx^2+ex^3)/(x-d)^3",
      {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 40.0},
      {3.0, 4.0, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 40.0},
      "( [0] + [1]*x + [2]*x*x + [4]*x*x*x ) / ( (x-[3])*(x-[3])*(x-[3]) )",
      3.0, 40.0,
      {1.0, 0.0, 0.0, -1.0, 0.0}});

  // B: same nominal binning + double exponential
  cs.push_back(Candidate{
      "B_nominal_2exp",
      "nominal 13/8-bin, double-exponential a*e^-bx + c*e^-dx + e",
      {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 40.0},
      {3.0, 4.0, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 40.0},
      "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]",
      3.0, 40.0,
      {0.5, 0.5, 0.5, 0.05, 0.5}});

  // C: nominal binning + simple power-law: a*x^b + c
  cs.push_back(Candidate{
      "C_nominal_powerLaw",
      "nominal 13/8-bin, power-law a*x^b + c",
      {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 40.0},
      {3.0, 4.0, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 40.0},
      "[0]*TMath::Power(x,[1]) + [2]",
      3.0, 40.0,
      {1.0, 0.0, 0.0}});

  // D: coarse binning + rational
  cs.push_back(Candidate{
      "D_coarse_rational",
      "coarse 6/6-bin, rational",
      {6.5, 9.0, 12.0, 15.0, 20.0, 25.0, 40.0},
      {3.0, 4.0, 5.0, 6.5, 9.0, 12.0, 40.0},
      "( [0] + [1]*x + [2]*x*x + [4]*x*x*x ) / ( (x-[3])*(x-[3])*(x-[3]) )",
      3.0, 40.0,
      {1.0, 0.0, 0.0, -1.0, 0.0}});

  // E: coarse binning + double exponential
  cs.push_back(Candidate{
      "E_coarse_2exp",
      "coarse 6/6-bin, double-exponential",
      {6.5, 9.0, 12.0, 15.0, 20.0, 25.0, 40.0},
      {3.0, 4.0, 5.0, 6.5, 9.0, 12.0, 40.0},
      "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]",
      3.0, 40.0,
      {0.5, 0.5, 0.5, 0.05, 0.5}});

  // F: HIN-16-025-style ultra-coarse fwd, coarse mid + rational
  // HIN-16-025 NP fwd is a single 3-30 GeV bin; mid uses 6.5/9.5/16/30 style.
  cs.push_back(Candidate{
      "F_HINlike_rational",
      "HIN-16-025-style ultra-coarse, rational (mid 6-bin, fwd 3-bin)",
      {6.5, 9.0, 12.0, 15.0, 20.0, 25.0, 40.0},
      {3.5, 6.5, 12.0, 40.0},
      "( [0] + [1]*x + [2]*x*x + [4]*x*x*x ) / ( (x-[3])*(x-[3])*(x-[3]) )",
      3.5, 40.0,
      {1.0, 0.0, 0.0, -1.0, 0.0}});

  return cs;
}

// ---- yield / fraction lookup helpers ------------------------------------

struct YieldErr
{
  double val = 0.0;
  double err = 0.0;
  bool ok = false;
};

// PbPb yield: Macros/Jpsi_250423/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_*_PRw_Effw0_Accw0_PtW0_TnP0.root
YieldErr GetPbPbYield(double ptLow, double ptHigh, double yLow, double yHigh)
{
  YieldErr r;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, 0, 180);
  TString path = Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Mass/"
                      "Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                      kRepo, kineLabel.Data());
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie())
    return r;
  TH1D *h = dynamic_cast<TH1D *>(f->Get("fitResults"));
  if (!h)
    return r;
  r.val = h->GetBinContent(1);
  r.err = h->GetBinError(1);
  r.ok = (r.val > 0.0);
  return r;
}

double GetPbPbFrac(double ptLow, double ptHigh, double yLow, double yHigh, double *errOut = nullptr)
{
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, 0, 180);
  TString path = Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Final/"
                      "2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                      kRepo, kineLabel.Data());
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie())
    return 0.0;
  TH1D *h = dynamic_cast<TH1D *>(f->Get("2DfitResults"));
  if (!h)
    return 0.0;
  if (errOut)
    *errOut = h->GetBinError(1);
  return h->GetBinContent(1);
}

// pp yield: Macros/Jpsi_L_cut/roots_1S_pp/{PRMC,NPMC}/Mass_FixedFitResult_*.root
YieldErr GetPpYield(bool prompt, double ptLow, double ptHigh, double yLow, double yHigh)
{
  YieldErr r;
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);
  const char *sub = prompt ? "PRMC" : "NPMC";
  TString path = Form("%s/Macros/Jpsi_L_cut/roots_1S_pp/%s/"
                      "Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                      kRepo, sub, kineLabel.Data());
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie())
    return r;
  TH1D *h = dynamic_cast<TH1D *>(f->Get("fitResults"));
  if (!h)
    return r;
  r.val = h->GetBinContent(1);
  r.err = h->GetBinError(1);
  r.ok = (r.val > 0.0);
  return r;
}

double GetPpFrac(bool prompt, double ptLow, double ptHigh, double yLow, double yHigh,
                 double *errOut = nullptr)
{
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);
  const char *sub = prompt ? "PRMC" : "NPMC";
  TString path = Form("%s/Macros/Jpsi_L_cut/roots_1S_pp/%s/"
                      "2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                      kRepo, sub, kineLabel.Data());
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie())
  {
    // fallback to PRMC dir for fraction when NPMC dir lacks 2D
    path = Form("%s/Macros/Jpsi_L_cut/roots_1S_pp/PRMC/"
                "2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                kRepo, kineLabel.Data());
    f.reset(TFile::Open(path, "READ"));
    if (!f || f->IsZombie())
      return 0.0;
  }
  TH1D *h = dynamic_cast<TH1D *>(f->Get("2DfitResults"));
  if (!h)
    return 0.0;
  if (errOut)
    *errOut = h->GetBinError(1);
  return h->GetBinContent(1);
}

// ---- MC dN/dpT cache ----------------------------------------------------

struct McSample
{
  TString chainPath;
  double massLow;
  double massHigh;
  double yLow;
  double yHigh;
  // 0..50 GeV, 2000 bins -> 0.025 GeV resolution. Plenty for rebinning.
  std::unique_ptr<TH1D> hFine;
  bool loaded = false;
};

void LoadMc(McSample &s)
{
  if (s.loaded)
    return;
  s.loaded = true;
  s.hFine.reset(new TH1D(Form("hFine_%p", (void *)&s), "; p_{T} (GeV/c); ", 2000, 0.0, 50.0));
  s.hFine->SetDirectory(nullptr);

  TChain tree("mmepevt");
  tree.Add(s.chainPath.Data());
  if (tree.GetEntries() <= 0)
  {
    std::cerr << "[mc] no entries in " << s.chainPath << "\n";
    return;
  }

  constexpr int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu];
  Int_t nDimu = 0;

  tree.SetBranchAddress("nDimu", &nDimu);
  tree.SetBranchAddress("mass", mass);
  tree.SetBranchAddress("y", y);
  tree.SetBranchAddress("pt", pt);

  const Long64_t n = tree.GetEntries();
  std::cout << "[mc] " << s.chainPath << " entries=" << n << "\n";
  for (Long64_t i = 0; i < n; ++i)
  {
    tree.GetEntry(i);
    for (int j = 0; j < nDimu; ++j)
    {
      const double absY = std::fabs(y[j]);
      if (mass[j] <= s.massLow || mass[j] >= s.massHigh)
        continue;
      if (pt[j] <= 0.0 || pt[j] >= 50.0)
        continue;
      if (absY <= s.yLow || absY >= s.yHigh)
        continue;
      s.hFine->Fill(pt[j]);
    }
  }
  std::cout << "[mc] integral = " << s.hFine->Integral() << "\n";
}

// Rebin the fine MC histogram into the candidate's bin edges.
TH1D *RebinMc(const TH1D *hFine, const std::vector<double> &edges, const char *name)
{
  TH1D *h = new TH1D(name, "; p_{T} (GeV/c); ", static_cast<int>(edges.size()) - 1, edges.data());
  h->Sumw2();
  // sum fine-bin counts within each [edges[i], edges[i+1])
  for (int ix = 1; ix <= hFine->GetNbinsX(); ++ix)
  {
    const double x = hFine->GetXaxis()->GetBinCenter(ix);
    const double w = hFine->GetBinContent(ix);
    if (w <= 0)
      continue;
    if (x < edges.front() || x >= edges.back())
      continue;
    h->Fill(x, w);
  }
  return h;
}

// ---- one (system, component, rap) component for one candidate ----------

bool ProcessComponent(const Candidate &c, bool isPbPb, bool prompt, bool fwd,
                      McSample &mc, const TString &outDir)
{
  LoadMc(mc);
  if (!mc.hFine || mc.hFine->Integral() <= 0)
    return false;

  const std::vector<double> &edges = fwd ? c.ptFwd : c.ptMid;
  if (edges.size() < 2)
    return false;

  const double yLow = fwd ? 1.6 : 0.0;
  const double yHigh = fwd ? 2.4 : 1.6;
  const TString rapLab = fwd ? "y1p6_2p4" : "y0_1p6";
  const TString sysLab = isPbPb ? "AA" : "pp";
  const TString compLab = prompt ? "Jpsi" : "BtoJpsi";
  const TString tagSuffix = isPbPb ? "260515" : "noCtau_260515";

  const int nb = static_cast<int>(edges.size()) - 1;
  std::unique_ptr<TH1D> hData(new TH1D(Form("hptData_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                       "; p_{T} (GeV/c); ", nb, edges.data()));
  std::unique_ptr<TH1D> hDataRatio(new TH1D(Form("hptDataR_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data()),
                                            "; p_{T} (GeV/c); ", nb, edges.data()));
  hData->SetDirectory(nullptr);
  hDataRatio->SetDirectory(nullptr);
  hData->Sumw2();
  hDataRatio->Sumw2();

  for (int ib = 1; ib <= nb; ++ib)
  {
    const double ptLow = edges[ib - 1];
    const double ptHigh = edges[ib];

    YieldErr tot;
    double w = 1.0;
    if (isPbPb)
    {
      // PbPb fits are inclusive in prompt/nonprompt; split via b-fraction.
      tot = GetPbPbYield(ptLow, ptHigh, yLow, yHigh);
      const double frac = GetPbPbFrac(ptLow, ptHigh, yLow, yHigh);
      w = prompt ? (1.0 - frac) : frac;
    }
    else
    {
      // pp uses PRMC / NPMC component-specific yield templates -> no extra split.
      tot = GetPpYield(prompt, ptLow, ptHigh, yLow, yHigh);
    }
    if (!tot.ok)
    {
      std::cerr << "[skip] " << c.tag << " " << sysLab << " " << compLab << " " << rapLab
                << " pt=" << ptLow << "-" << ptHigh << " : no yield fit ROOT\n";
      continue;
    }
    hData->SetBinContent(ib, tot.val * w);
    hData->SetBinError(ib, tot.err * std::max(w, 1e-3));
    hDataRatio->SetBinContent(ib, tot.val * w);
    hDataRatio->SetBinError(ib, tot.err * std::max(w, 1e-3));
  }

  if (hData->Integral() <= 0.0)
  {
    std::cerr << "[skip] " << c.tag << " " << sysLab << " " << compLab << " " << rapLab
              << " : zero data integral\n";
    return false;
  }

  std::unique_ptr<TH1D> hMc(RebinMc(mc.hFine.get(), edges,
                                    Form("hptMC_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data())));
  std::unique_ptr<TH1D> hMcRatio(RebinMc(mc.hFine.get(), edges,
                                         Form("hptMCR_%s_%s_%s", sysLab.Data(), compLab.Data(), rapLab.Data())));
  hMc->SetDirectory(nullptr);
  hMcRatio->SetDirectory(nullptr);
  if (hMc->Integral() <= 0 || hMcRatio->Integral() <= 0)
  {
    std::cerr << "[skip] " << c.tag << " " << sysLab << " " << compLab << " " << rapLab
              << " : zero MC integral\n";
    return false;
  }

  hData->Scale(1.0 / hData->Integral());
  hDataRatio->Scale(1.0 / hDataRatio->Integral());
  hMc->Scale(1.0 / hMc->Integral());
  hMcRatio->Scale(1.0 / hMcRatio->Integral());
  TH1ScaleByWidth(hData.get());
  TH1ScaleByWidth(hMc.get());

  hDataRatio->Divide(hMcRatio.get());

  std::unique_ptr<TF1> fit(new TF1(Form("fitRatio1_%s", c.tag.Data()),
                                   c.fitFormula.Data(), c.fitMin, c.fitMax));
  fit->SetNpx(1000);
  for (std::size_t ip = 0; ip < c.seedPars.size(); ++ip)
    fit->SetParameter(static_cast<int>(ip), c.seedPars[ip]);
  hDataRatio->Fit(fit.get(), "IEQ", "", c.fitMin, c.fitMax);
  TFitResultPtr res = hDataRatio->Fit(fit.get(), "SQ", "", c.fitMin, c.fitMax);
  if (res.Get())
    std::cout << "[fit] " << c.tag << " " << sysLab << " " << compLab << " " << rapLab
              << " chi2/ndf = " << res->Chi2() << "/" << res->Ndf() << "\n";

  // ---- write output ROOT (in candidate dir) ----
  TString outName;
  if (isPbPb)
    outName = Form("%s/ratioDataMC_%s_%s_DATA_%s_%s.root", outDir.Data(),
                   sysLab.Data(), compLab.Data(), rapLab.Data(), "260515");
  else
    outName = Form("%s/ratioDataMC_%s_%s_DATA_noCtau_%s_%s.root", outDir.Data(),
                   sysLab.Data(), compLab.Data(), rapLab.Data(), "260515");

  std::unique_ptr<TFile> out(TFile::Open(outName, "RECREATE"));
  if (!out || out->IsZombie())
  {
    std::cerr << "[write] cannot open " << outName << "\n";
    return false;
  }
  out->cd();
  hDataRatio->SetName("WeightFactor");
  hDataRatio->Write("WeightFactor", TObject::kOverwrite);
  fit->SetName("dataMC_Ratio1");
  fit->Write("dataMC_Ratio1", TObject::kOverwrite);
  fit->SetName("fitRatio1");
  fit->Write("fitRatio1", TObject::kOverwrite);
  hData->Write("hptData_norm_width");
  hMc->Write("hptMC_norm_width");

  // sidecar PDF / PNG so we can eyeball each candidate quickly
  TCanvas cnv("c_ptw", c.tag.Data(), 600, 700);
  TPad p1("p1", "p1", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.30);
  p1.Draw();
  p2.Draw();
  p1.cd();
  hData->SetMarkerStyle(20);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hMc->SetLineColor(kRed + 1);
  hMc->SetLineWidth(2);
  const double yMax = std::max(hData->GetMaximum(), hMc->GetMaximum()) * 1.3;
  hData->SetMaximum(yMax > 0 ? yMax : 1.0);
  hData->Draw("e1");
  hMc->Draw("hist same");
  p2.cd();
  hDataRatio->SetMarkerStyle(20);
  hDataRatio->SetMarkerColor(kBlack);
  hDataRatio->SetLineColor(kBlack);
  hDataRatio->SetMinimum(0.0);
  hDataRatio->SetMaximum(4.0);
  hDataRatio->Draw("e1");
  fit->SetLineColor(kBlue + 1);
  fit->SetLineWidth(2);
  fit->Draw("same");
  TLine unity(c.fitMin, 1.0, c.fitMax, 1.0);
  unity.SetLineStyle(2);
  unity.Draw("same");
  cnv.Write("canvas_A");

  TNamed note("candidate_label", c.label.Data());
  note.Write();
  TNamed binNote("binning", Form("mid=%d-bin / fwd=%d-bin",
                                  static_cast<int>(c.ptMid.size()) - 1,
                                  static_cast<int>(c.ptFwd.size()) - 1));
  binNote.Write();
  TNamed fitNote("fitFormula", c.fitFormula.Data());
  fitNote.Write();

  out->Close();

  TString pdfPath = outName;
  pdfPath.ReplaceAll(".root", ".pdf");
  cnv.SaveAs(pdfPath);
  TString pngPath = outName;
  pngPath.ReplaceAll(".root", ".png");
  cnv.SaveAs(pngPath);
  std::cout << "[wrote] " << outName << "\n";
  return true;
}

void RunOneCandidate(const Candidate &c,
                     McSample &mcPpPr, McSample &mcPpNp,
                     McSample &mcAaPr, McSample &mcAaNp,
                     bool midOnly = false)
{
  const TString outDir = Form("%s/compareDataToMC/%s/%s", kRepo, kOutBase, c.tag.Data());
  if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir))
  {
    std::cerr << "[err] cannot mkdir " << outDir << "\n";
    return;
  }
  std::cout << "\n========== candidate " << c.tag
            << (midOnly ? " (mid only)" : "")
            << " : " << c.label << " ==========\n";
  // pp prompt
  ProcessComponent(c, false, true, false, mcPpPr, outDir);
  if (!midOnly)
    ProcessComponent(c, false, true, true, mcPpPr, outDir);
  // pp nonprompt
  ProcessComponent(c, false, false, false, mcPpNp, outDir);
  if (!midOnly)
    ProcessComponent(c, false, false, true, mcPpNp, outDir);
  // PbPb prompt
  ProcessComponent(c, true, true, false, mcAaPr, outDir);
  if (!midOnly)
    ProcessComponent(c, true, true, true, mcAaPr, outDir);
  // PbPb nonprompt
  ProcessComponent(c, true, false, false, mcAaNp, outDir);
  if (!midOnly)
    ProcessComponent(c, true, false, true, mcAaNp, outDir);
}

}  // namespace PtWCandidates260515

void make_ptw_candidates_260515(const char *only = "", bool midOnly = false)
{
  using namespace PtWCandidates260515;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  // MC samples - load once and re-use across candidates
  McSample mcPpPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcPpNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_miniAOD_251103.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcAaPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcAaNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_NonPrompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};

  auto cs = BuildCandidates();
  for (const auto &c : cs)
  {
    if (only && std::strlen(only) > 0 && c.tag != only)
      continue;
    RunOneCandidate(c, mcPpPr, mcPpNp, mcAaPr, mcAaNp, midOnly);
  }
}
