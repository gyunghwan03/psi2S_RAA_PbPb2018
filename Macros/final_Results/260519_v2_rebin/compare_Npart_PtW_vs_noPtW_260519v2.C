// =============================================================
//  compare_Npart_PtW_vs_noPtW_260519v2.C
//
//  Overlays Npart-axis J/psi RAA for two Eff/Acc treatments:
//    (a) PtW ON  : G_aggregateHighPt_2exp full-stat (option 2)
//                  Eff_Acc_260519_v2/roots/G_aggregate_full
//    (b) PtW OFF : full-stat without pT reweighting
//                  Eff_Acc_260515/roots/noPtW_full
//
//  Both treatments use the same pT-integrated yield + b-fraction
//  inputs and the standard "single integrated eff*acc" formula.
//  The HIN-16-025 reference is drawn for context.
//
//  4 panels (PR mid, PR fwd, NP mid, NP fwd) on one canvas.
//  Output: figs/compare_Npart_PtW_vs_noPtW_260519v2.pdf / .png
// =============================================================

#include "../../../commonUtility.h"
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "../../../Style.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace {
constexpr const char *kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
constexpr double Nmb = 11968044281.;
constexpr double lumi_pp = 3.002;
constexpr double lumi_pp_scale = 1e-9;

// Nominal PtW-on Eff/Acc: H_NP_rational (option 2 aggregated bins + 2-exp PR
// + rational NP fit). PtW-off baseline unchanged.
constexpr const char *kPtWDir = "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/H_NP_rational_full";
constexpr const char *kPtWTag = "260520_H_NP_rational";
constexpr const char *kNoPtWDir = "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/noPtW_full";
constexpr const char *kNoPtWTag = "260518full_noPtW";

struct V { double v = 0, e = 0; };

V readBin1(const TString &p, const char *hn) {
  V o;
  std::unique_ptr<TFile> f(TFile::Open(p, "READ"));
  if (!f || f->IsZombie()) return o;
  TH1D *h = (TH1D*)f->Get(hn);
  if (!h) return o;
  o.v = h->GetBinContent(1); o.e = h->GetBinError(1);
  return o;
}

V GetYieldPbPb(double pl, double ph, double yl, double yh, int cL, int cH) {
  TString lab = getKineLabel(pl, ph, yl, yh, 0.0, cL, cH);
  return readBin1(Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kRepo, lab.Data()), "fitResults");
}
V GetFracPbPb(double pl, double ph, double yl, double yh, int cL, int cH) {
  TString lab = getKineLabel(pl, ph, yl, yh, 0.0, cL, cH);
  return readBin1(Form("%s/Macros/Jpsi_250423/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kRepo, lab.Data()), "2DfitResults");
}
V GetYieldPP(double pl, double ph, double yl, double yh) {
  TString lab = getKineLabelpp(pl, ph, yl, yh, 0.0);
  return readBin1(Form("%s/Macros/pp_Jpsi/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kRepo, lab.Data()), "fitResults");
}
V GetFracPP(double pl, double ph, double yl, double yh) {
  TString lab = getKineLabelpp(pl, ph, yl, yh, 0.0);
  return readBin1(Form("%s/Macros/pp_Jpsi/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kRepo, lab.Data()), "2DfitResults");
}

double rel(double e, double v) { return v != 0 ? e/v : 0.0; }

double integratedPpEff(TFile *f, const char *numN, const char *denN) {
  if (!f) return 1.0;
  TH1 *num = (TH1*)f->Get(numN);
  TH1 *den = (TH1*)f->Get(denN);
  if (!num || !den) return 1.0;
  double sN = 0, sD = 0;
  for (int i = 1; i <= num->GetNbinsX(); ++i) sN += num->GetBinContent(i);
  for (int i = 1; i <= den->GetNbinsX(); ++i) sD += den->GetBinContent(i);
  return sD > 0 ? sN/sD : 1.0;
}

struct RAAOut {
  std::vector<double> raaPR, raaPR_err, raaNP, raaNP_err;
  std::vector<double> npart;  // central -> peripheral order
};

// Compute one (rap, PtW-state) Npart RAA.  Returns 6 mid or 4 fwd points.
RAAOut compute(bool isMid, bool isPtWOn)
{
  RAAOut out;

  std::vector<int> centHi;
  std::vector<double> centBin, taa, npart_central_first;
  double ptLow, ptHigh, yLow, yHigh;
  const char *ptTag, *absY, *accRap;
  int nCent;
  if (isMid) {
    nCent = 6;
    centBin = {0, 10, 20, 30, 40, 50, 90};
    centHi = {0, 20, 40, 60, 80, 100, 180};
    taa = {23.05, 14.39, 8.798, 5.124, 2.777, 0.5803};
    npart_central_first = {356.9, 262.3, 188.2, 131.0, 87.19, 27.12};
    ptLow = 6.5; ptHigh = 40.0; yLow = 0.0; yHigh = 1.6;
    ptTag = "pt_6p5_to_40"; absY = "absy0_1p6"; accRap = "midy";
  } else {
    nCent = 4;
    centBin = {0, 10, 30, 50, 90};
    centHi = {0, 20, 60, 100, 180};
    taa = {23.05, 11.60, 3.950, 0.5803};
    npart_central_first = {356.9, 225.2, 109.1, 27.12};
    ptLow = 3.5; ptHigh = 40.0; yLow = 1.6; yHigh = 2.4;
    ptTag = "pt_3_to_40"; absY = "absy1p6_2p4"; accRap = "Fory";
  }

  const TString dir = isPtWOn ? kPtWDir : kNoPtWDir;
  const TString tag = isPtWOn ? kPtWTag : kNoPtWTag;
  const TString suffix = isPtWOn ? ".root" : "_noPtW.root";
  const TString ptwLabel = isPtWOn ? "PtWnomi" : "PtW1";
  const char *accWgt = isPtWOn ? "wgt1" : "wgt0";

  TFile *fEff_AAPR = TFile::Open(Form("%s/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_%s_tnp1_%s%s",
                                       dir.Data(), ptwLabel.Data(), tag.Data(), suffix.Data()));
  TFile *fEff_AANP = TFile::Open(Form("%s/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_%s_tnp1_%s%s",
                                       dir.Data(), ptwLabel.Data(), tag.Data(), suffix.Data()));
  TFile *fEff_ppPR = TFile::Open(Form("%s/mc_eff_vs_pt_rap_prompt_pp_Jpsi_%s_tnp1_%s%s",
                                       dir.Data(), ptwLabel.Data(), tag.Data(), suffix.Data()));
  TFile *fEff_ppNP = TFile::Open(Form("%s/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_%s_tnp1_%s%s",
                                       dir.Data(), ptwLabel.Data(), tag.Data(), suffix.Data()));
  TFile *fAcc_AAPR = TFile::Open(Form("%s/acceptance_PromptJpsi_GenOnly_%s_PbPb_SysUp0_%s.root",
                                       dir.Data(), accWgt, tag.Data()));
  TFile *fAcc_AANP = TFile::Open(Form("%s/acceptance_BtoJpsi_GenOnly_%s_PbPb_SysUp0_%s.root",
                                       dir.Data(), accWgt, tag.Data()));
  TFile *fAcc_ppPR = TFile::Open(Form("%s/acceptance_PromptJpsi_GenOnly_%s_pp_SysUp0_%s.root",
                                       dir.Data(), accWgt, tag.Data()));
  TFile *fAcc_ppNP = TFile::Open(Form("%s/acceptance_BtoJpsi_GenOnly_%s_pp_SysUp0_%s.root",
                                       dir.Data(), accWgt, tag.Data()));
  if (!fEff_AAPR || !fEff_AANP || !fEff_ppPR || !fEff_ppNP ||
      !fAcc_AAPR || !fAcc_AANP || !fAcc_ppPR || !fAcc_ppNP) {
    std::cerr << "[compute] missing Eff/Acc file(s)\n";
    return out;
  }

  TH1D *hEff_AAPR = (TH1D*)fEff_AAPR->Get(Form("mc_eff_vs_cent_TnP1_PtW1_%s_%s", ptTag, absY));
  TH1D *hEff_AANP = (TH1D*)fEff_AANP->Get(Form("mc_eff_vs_cent_TnP1_PtW1_%s_%s", ptTag, absY));
  const char *numKey = isMid ? "hist_eff_num_mid" : "hist_eff_num_fwd";
  const char *denKey = isMid ? "hist_eff_den_mid" : "hist_eff_den_fwd";
  double eff_pp_PR_int = integratedPpEff(fEff_ppPR, numKey, denKey);
  double eff_pp_NP_int = integratedPpEff(fEff_ppNP, numKey, denKey);
  TH1D *hAcc_AAPR = (TH1D*)fAcc_AAPR->Get(Form("hAccPt_2021_%s_Int", accRap));
  TH1D *hAcc_AANP = (TH1D*)fAcc_AANP->Get(Form("hAccPt_2021_%s_Int", accRap));
  TH1D *hAcc_ppPR = (TH1D*)fAcc_ppPR->Get(Form("hAccPt_2021_%s_Int", accRap));
  TH1D *hAcc_ppNP = (TH1D*)fAcc_ppNP->Get(Form("hAccPt_2021_%s_Int", accRap));

  V yPp = GetYieldPP(ptLow, ptHigh, yLow, yHigh);
  V frPp = GetFracPP(ptLow, ptHigh, yLow, yHigh);
  double yPpPR = yPp.v * (1 - frPp.v);
  double yPpNP = yPp.v * frPp.v;
  double yPpPR_err = yPpPR * std::sqrt(std::pow(rel(yPp.e, yPp.v), 2) + std::pow(rel(frPp.e, 1 - frPp.v), 2));
  double yPpNP_err = yPpNP * std::sqrt(std::pow(rel(yPp.e, yPp.v), 2) + std::pow(rel(frPp.e, frPp.v), 2));
  double acc_pp_PR = hAcc_ppPR ? hAcc_ppPR->GetBinContent(1) : 1.0;
  double acc_pp_NP = hAcc_ppNP ? hAcc_ppNP->GetBinContent(1) : 1.0;
  double wPP_PR = eff_pp_PR_int * acc_pp_PR;
  double wPP_NP = eff_pp_NP_int * acc_pp_NP;
  double X_pp_PR = lumi_pp_scale * (yPpPR / wPP_PR) /
                   (lumi_pp * 1e2 * (ptHigh - ptLow) * 2 * (yHigh - yLow));
  double X_pp_NP = lumi_pp_scale * (yPpNP / wPP_NP) /
                   (lumi_pp * 1e2 * (ptHigh - ptLow) * 2 * (yHigh - yLow));

  for (int c = 0; c < nCent; ++c) {
    V yA = GetYieldPbPb(ptLow, ptHigh, yLow, yHigh, centHi[c], centHi[c + 1]);
    V frA = GetFracPbPb(ptLow, ptHigh, yLow, yHigh, centHi[c], centHi[c + 1]);
    double yAAPR = yA.v * (1 - frA.v);
    double yAANP = yA.v * frA.v;
    double yAAPR_err = yAAPR * std::sqrt(std::pow(rel(yA.e, yA.v), 2) + std::pow(rel(frA.e, 1 - frA.v), 2));
    double yAANP_err = yAANP * std::sqrt(std::pow(rel(yA.e, yA.v), 2) + std::pow(rel(frA.e, frA.v), 2));

    double effAA_PR = hEff_AAPR ? hEff_AAPR->GetBinContent(c + 1) : 0;
    double effAA_NP = hEff_AANP ? hEff_AANP->GetBinContent(c + 1) : 0;
    double accAA_PR = hAcc_AAPR ? hAcc_AAPR->GetBinContent(1) : 0;
    double accAA_NP = hAcc_AANP ? hAcc_AANP->GetBinContent(1) : 0;
    if (effAA_PR <= 0 || effAA_NP <= 0 || accAA_PR <= 0 || accAA_NP <= 0) {
      out.raaPR.push_back(0); out.raaPR_err.push_back(0);
      out.raaNP.push_back(0); out.raaNP_err.push_back(0);
      out.npart.push_back(npart_central_first[c]);
      continue;
    }
    double wAA_PR = effAA_PR * accAA_PR;
    double wAA_NP = effAA_NP * accAA_NP;
    double centFrac = (centBin[c + 1] - centBin[c]) / 90.0;
    double X_AA_PR = (yAAPR / wAA_PR) /
                     (Nmb * taa[c] * (ptHigh - ptLow) * 2 * (yHigh - yLow) * centFrac);
    double X_AA_NP = (yAANP / wAA_NP) /
                     (Nmb * taa[c] * (ptHigh - ptLow) * 2 * (yHigh - yLow) * centFrac);
    double raaPR = X_pp_PR > 0 ? X_AA_PR / X_pp_PR : 0;
    double raaNP = X_pp_NP > 0 ? X_AA_NP / X_pp_NP : 0;
    double raaPR_err = raaPR * std::sqrt(std::pow(rel(yAAPR_err, yAAPR), 2) +
                                          std::pow(rel(yPpPR_err, yPpPR), 2));
    double raaNP_err = raaNP * std::sqrt(std::pow(rel(yAANP_err, yAANP), 2) +
                                          std::pow(rel(yPpNP_err, yPpNP), 2));
    out.raaPR.push_back(raaPR);
    out.raaPR_err.push_back(raaPR_err);
    out.raaNP.push_back(raaNP);
    out.raaNP_err.push_back(raaNP_err);
    out.npart.push_back(npart_central_first[c]);
  }
  return out;
}

// HIN-16-025 references (central-first index → we'll plot at given Npart).
struct HINNpart { std::vector<double> npart, raa, err; };
HINNpart LoadHIN(const char *file, const char *histName, const char *errName,
                  const std::vector<double> &npartCentralFirst) {
  HINNpart o;
  std::unique_ptr<TFile> f(TFile::Open(Form("%s/Macros/final_Results/roots/%s", kRepo, file), "READ"));
  if (!f || f->IsZombie()) return o;
  TH1 *h = (TH1*)f->Get(histName);
  TH1 *he = (TH1*)f->Get(errName);
  if (!h || !he) return o;
  int n = std::min(h->GetNbinsX(), (int)npartCentralFirst.size());
  for (int i = 0; i < n; ++i) {
    o.npart.push_back(npartCentralFirst[i]);
    o.raa.push_back(h->GetBinContent(i + 1));
    o.err.push_back(he->GetBinContent(i + 1));
  }
  return o;
}

TGraphErrors *makeGraph(const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &ey,
                         int color, int marker, double size = 1.4)
{
  std::vector<double> ex(x.size(), 0.0);
  TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(),
                                     ex.data(), ey.data());
  g->SetMarkerColor(color); g->SetLineColor(color);
  g->SetMarkerStyle(marker); g->SetMarkerSize(size); g->SetLineWidth(2);
  return g;
}

TGraphErrors *makeHinGraph(const HINNpart &h) {
  return makeGraph(h.npart, h.raa, h.err, kGray + 3, 28, 1.2);
}
}  // namespace

void compare_Npart_PtW_vs_noPtW_260519v2()
{
  gStyle->SetOptStat(0);
  setTDRStyle();
  int iPeriod = 101, iPos = 33;

  std::cout << "\n== Computing Npart RAA: mid + PtW on ==\n";
  RAAOut mid_on  = compute(true, true);
  std::cout << "== Computing Npart RAA: mid + PtW off ==\n";
  RAAOut mid_off = compute(true, false);
  std::cout << "== Computing Npart RAA: fwd + PtW on ==\n";
  RAAOut fwd_on  = compute(false, true);
  std::cout << "== Computing Npart RAA: fwd + PtW off ==\n";
  RAAOut fwd_off = compute(false, false);

  const std::vector<double> np_mid_HIN = {358.8, 264.2, 189.2, 131.4, 86.9, 21.9};
  const std::vector<double> np_fwdPR_HIN = {311.5, 160.3, 32.7};
  // Nonprompt fwd in HIN-16-025 is published in two pT slices, both at the
  // mid-rapidity Npart binning (6 points): Table 30/y2 = 6.5<pT<50 GeV/c,
  // Table 32/y1 = 3<pT<6.5 GeV/c.  Both at 1.8<|y|<2.4.
  HINNpart hin_prMid = LoadHIN("RAA_PR_JPsi_HIN_16_025_mid_Npart.root",
                                "Table 15/Hist1D_y1", "Table 15/Hist1D_y1_e1", np_mid_HIN);
  HINNpart hin_prFwd = LoadHIN("RAA_PR_JPsi_HIN_16_025_fwd_Npart.root",
                                "Table 19/Hist1D_y1", "Table 19/Hist1D_y1_e1", np_fwdPR_HIN);
  HINNpart hin_npMid = LoadHIN("RAA_NP_JPsi_HIN_16_025_Npart.root",
                                "Table 30/Hist1D_y1", "Table 30/Hist1D_y1_e1", np_mid_HIN);
  HINNpart hin_npFwd_hipt = LoadHIN("RAA_NP_JPsi_HIN_16_025_Npart.root",
                                     "Table 30/Hist1D_y2", "Table 30/Hist1D_y2_e1", np_mid_HIN);
  HINNpart hin_npFwd_lopt = LoadHIN("RAA_NP_JPsi_HIN_16_025_fwd_Npart.root",
                                     "Table 32/Hist1D_y1", "Table 32/Hist1D_y1_e1", np_mid_HIN);

  // ----- print a quick comparison table -----
  std::cout << "\n=== Npart RAA: PtW on vs PtW off vs HIN-16-025 ===\n";
  std::cout << "             |       PR mid       |       NP mid       |\n";
  std::cout << "  Npart      | PtW_on  PtW_off HIN| PtW_on  PtW_off HIN|\n";
  for (size_t i = 0; i < mid_on.npart.size(); ++i) {
    double hinPR = 0, hinNP = 0;
    for (size_t j = 0; j < hin_prMid.npart.size(); ++j)
      if (std::abs(hin_prMid.npart[j] - mid_on.npart[i]) < 30) hinPR = hin_prMid.raa[j];
    for (size_t j = 0; j < hin_npMid.npart.size(); ++j)
      if (std::abs(hin_npMid.npart[j] - mid_on.npart[i]) < 30) hinNP = hin_npMid.raa[j];
    printf("  Npart=%6.1f | %.3f  %.3f  %.3f| %.3f  %.3f  %.3f|\n",
           mid_on.npart[i], mid_on.raaPR[i], mid_off.raaPR[i], hinPR,
                            mid_on.raaNP[i], mid_off.raaNP[i], hinNP);
  }
  std::cout << "\n             |  PR fwd               |       NP fwd (hipt / lopt)        |\n";
  for (size_t i = 0; i < fwd_on.npart.size(); ++i) {
    double hinPR = 0, hinNPhi = 0, hinNPlo = 0;
    for (size_t j = 0; j < hin_prFwd.npart.size(); ++j)
      if (std::abs(hin_prFwd.npart[j] - fwd_on.npart[i]) < 50) hinPR = hin_prFwd.raa[j];
    for (size_t j = 0; j < hin_npFwd_hipt.npart.size(); ++j)
      if (std::abs(hin_npFwd_hipt.npart[j] - fwd_on.npart[i]) < 30) hinNPhi = hin_npFwd_hipt.raa[j];
    for (size_t j = 0; j < hin_npFwd_lopt.npart.size(); ++j)
      if (std::abs(hin_npFwd_lopt.npart[j] - fwd_on.npart[i]) < 30) hinNPlo = hin_npFwd_lopt.raa[j];
    printf("  Npart=%6.1f | on=%.3f off=%.3f HIN=%.3f | on=%.3f off=%.3f HIN_hi=%.3f HIN_lo=%.3f\n",
           fwd_on.npart[i], fwd_on.raaPR[i], fwd_off.raaPR[i], hinPR,
                            fwd_on.raaNP[i], fwd_off.raaNP[i], hinNPhi, hinNPlo);
  }

  // ----- one canvas per (component, rapidity) -----
  // Accepts up to two HIN-16-025 series. For NP fwd we draw both the high-pT
  // (Table 30/y2, 6.5<pT<50) and the low-pT (Table 32/y1, 3<pT<6.5) HIN
  // slices, following the convention of compare_NP_fwd_Npart_Jpsi_260519v2.
  auto drawOne = [&](const char *outBase,
                     const char *compTxt, const char *rapTxt, const char *ptTxt,
                     const RAAOut &on, const RAAOut &off,
                     const HINNpart &hin1, const char *hin1Lbl,
                     const HINNpart *hin2, const char *hin2Lbl,
                     bool isPR)
  {
    TCanvas *c = new TCanvas(outBase, outBase, 700, 700);
    c->cd();
    c->SetTicks(1, 1);
    c->DrawFrame(0, 0, 400, 1.5, ";#LTN_{part}#GT;R_{AA}");
    TLine *one = new TLine(0, 1, 400, 1);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw("same");

    TGraphErrors *gOn  = makeGraph(on.npart,
                                    isPR ? on.raaPR : on.raaNP,
                                    isPR ? on.raaPR_err : on.raaNP_err,
                                    kRed + 1, 20, 1.4);
    TGraphErrors *gOff = makeGraph(off.npart,
                                    isPR ? off.raaPR : off.raaNP,
                                    isPR ? off.raaPR_err : off.raaNP_err,
                                    kBlue + 2, 21, 1.4);
    TGraphErrors *gHin1 = makeGraph(hin1.npart, hin1.raa, hin1.err, kGray + 3, 28, 1.2);
    TGraphErrors *gHin2 = nullptr;
    if (hin2)
      gHin2 = makeGraph(hin2->npart, hin2->raa, hin2->err, kGray + 1, 25, 1.2);
    if (gHin2) gHin2->Draw("p same");
    gHin1->Draw("p same");
    gOff->Draw("p same");
    gOn->Draw("p same");

    TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.78);
    leg->SetTextSize(0.030); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(gOn,  "PtW on", "pe");
    leg->AddEntry(gOff, "noPtW", "pe");
    leg->AddEntry(gHin1, hin1Lbl, "pe");
    if (gHin2) leg->AddEntry(gHin2, hin2Lbl, "pe");
    leg->Draw("same");

    // Rapidity / pT / component text on the upper-left
    const float pos_x = 0.21;
    const float pos_y = 0.86;
    const float dy = 0.055;
    const int color = 1;
    const float size = 22;
    drawText(compTxt, pos_x, pos_y,              color, size);
    drawText(rapTxt,  pos_x, pos_y - dy,         color, size);
    drawText(ptTxt,   pos_x, pos_y - 2 * dy,     color, size);

    CMS_lumi_v2mass(c, iPeriod, iPos);

    c->SaveAs(Form("figs/%s.pdf", outBase));
    c->SaveAs(Form("figs/%s.png", outBase));
    std::cout << "  [wrote] figs/" << outBase << ".pdf / .png\n";
  };

  std::cout << "\n";
  drawOne("compare_Npart_PtW_vs_noPtW_PR_mid_260519v2",
          "Prompt J/#psi", "|y| < 1.6", "6.5 < p_{T} < 40 GeV/c",
          mid_on, mid_off,
          hin_prMid, "HIN-16-025, 6.5 < p_{T} < 30 GeV/c",
          nullptr, nullptr, true);
  drawOne("compare_Npart_PtW_vs_noPtW_PR_fwd_260519v2",
          "Prompt J/#psi", "1.6 < |y| < 2.4", "3.5 < p_{T} < 40 GeV/c",
          fwd_on, fwd_off,
          hin_prFwd, "HIN-16-025, 3 < p_{T} < 30 GeV/c, 1.8<|y|<2.4",
          nullptr, nullptr, true);
  drawOne("compare_Npart_PtW_vs_noPtW_NP_mid_260519v2",
          "Nonprompt J/#psi", "|y| < 1.6", "6.5 < p_{T} < 40 GeV/c",
          mid_on, mid_off,
          hin_npMid, "HIN-16-025, 6.5 < p_{T} < 50 GeV/c, |y|<0.6",
          nullptr, nullptr, false);
  drawOne("compare_Npart_PtW_vs_noPtW_NP_fwd_260519v2",
          "Nonprompt J/#psi", "1.6 < |y| < 2.4", "3.5 < p_{T} < 40 GeV/c",
          fwd_on, fwd_off,
          hin_npFwd_hipt, "HIN-16-025, 6.5 < p_{T} < 50, 1.8<|y|<2.4",
          &hin_npFwd_lopt, "HIN-16-025, 3 < p_{T} < 6.5, 1.8<|y|<2.4", false);
}
