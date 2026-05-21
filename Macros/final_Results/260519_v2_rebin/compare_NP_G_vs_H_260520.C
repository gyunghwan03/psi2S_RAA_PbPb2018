// =============================================================
//  compare_NP_G_vs_H_260520.C
//
//  Overlay J/psi RAA from two Eff/Acc treatments that differ ONLY
//  in the nonprompt (B->J/psi) pT-weight shape:
//    G : option-2 aggregated bins + 2-exp fit (current baseline)
//        Eff_Acc_260519_v2/roots/G_aggregate_full
//    H : option-2 aggregated bins + rational fit (NP only)
//        Eff_Acc_260519_v2/roots/H_NP_rational_full
//  PR (Jpsi) PtW is identical between G and H, so PR panels are
//  a sanity check (curves should overlap).
//
//  4 panels (Npart axis): PR mid, PR fwd, NP mid, NP fwd
//  4 panels (pT axis):    PR mid, PR fwd, NP mid, NP fwd
//
//  Outputs: figs/compare_NP_G_vs_H_260520_{PR,NP}_{mid,fwd}_{pT,Npart}.pdf
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

constexpr const char *kGDir = "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/G_aggregate_full";
constexpr const char *kGTag = "260519v2_G_aggregate";
constexpr const char *kHDir = "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/H_NP_rational_full";
constexpr const char *kHTag = "260520_H_NP_rational";

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
  std::vector<double> npart;
};

// Cent-bin Npart RAA for one Eff/Acc treatment.
RAAOut computeNpart(bool isMid, const char *dir, const char *tag)
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

  TFile *fEff_AAPR = TFile::Open(Form("%s/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_%s.root",  dir, tag));
  TFile *fEff_AANP = TFile::Open(Form("%s/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_%s.root", dir, tag));
  TFile *fEff_ppPR = TFile::Open(Form("%s/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_%s.root", dir, tag));
  TFile *fEff_ppNP = TFile::Open(Form("%s/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_%s.root", dir, tag));
  TFile *fAcc_AAPR = TFile::Open(Form("%s/acceptance_PromptJpsi_GenOnly_wgt1_PbPb_SysUp0_%s.root", dir, tag));
  TFile *fAcc_AANP = TFile::Open(Form("%s/acceptance_BtoJpsi_GenOnly_wgt1_PbPb_SysUp0_%s.root", dir, tag));
  TFile *fAcc_ppPR = TFile::Open(Form("%s/acceptance_PromptJpsi_GenOnly_wgt1_pp_SysUp0_%s.root", dir, tag));
  TFile *fAcc_ppNP = TFile::Open(Form("%s/acceptance_BtoJpsi_GenOnly_wgt1_pp_SysUp0_%s.root", dir, tag));

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
  double X_pp_PR = lumi_pp_scale * (yPpPR / (eff_pp_PR_int * acc_pp_PR)) /
                   (lumi_pp * 1e2 * (ptHigh - ptLow) * 2 * (yHigh - yLow));
  double X_pp_NP = lumi_pp_scale * (yPpNP / (eff_pp_NP_int * acc_pp_NP)) /
                   (lumi_pp * 1e2 * (ptHigh - ptLow) * 2 * (yHigh - yLow));

  for (int c = 0; c < nCent; ++c) {
    V yA = GetYieldPbPb(ptLow, ptHigh, yLow, yHigh, centHi[c], centHi[c+1]);
    V frA = GetFracPbPb(ptLow, ptHigh, yLow, yHigh, centHi[c], centHi[c+1]);
    double yAAPR = yA.v * (1 - frA.v);
    double yAANP = yA.v * frA.v;
    double yAAPR_err = yAAPR * std::sqrt(std::pow(rel(yA.e, yA.v), 2) + std::pow(rel(frA.e, 1 - frA.v), 2));
    double yAANP_err = yAANP * std::sqrt(std::pow(rel(yA.e, yA.v), 2) + std::pow(rel(frA.e, frA.v), 2));
    double effAA_PR = hEff_AAPR ? hEff_AAPR->GetBinContent(c+1) : 0;
    double effAA_NP = hEff_AANP ? hEff_AANP->GetBinContent(c+1) : 0;
    double accAA_PR = hAcc_AAPR ? hAcc_AAPR->GetBinContent(1) : 0;
    double accAA_NP = hAcc_AANP ? hAcc_AANP->GetBinContent(1) : 0;
    if (effAA_PR <= 0 || effAA_NP <= 0 || accAA_PR <= 0 || accAA_NP <= 0) {
      out.raaPR.push_back(0); out.raaPR_err.push_back(0);
      out.raaNP.push_back(0); out.raaNP_err.push_back(0);
      out.npart.push_back(npart_central_first[c]); continue;
    }
    double wAA_PR = effAA_PR * accAA_PR;
    double wAA_NP = effAA_NP * accAA_NP;
    double centFrac = (centBin[c+1] - centBin[c]) / 90.0;
    double X_AA_PR = (yAAPR / wAA_PR) / (Nmb * taa[c] * (ptHigh - ptLow) * 2 * (yHigh - yLow) * centFrac);
    double X_AA_NP = (yAANP / wAA_NP) / (Nmb * taa[c] * (ptHigh - ptLow) * 2 * (yHigh - yLow) * centFrac);
    double raaPR = X_pp_PR > 0 ? X_AA_PR / X_pp_PR : 0;
    double raaNP = X_pp_NP > 0 ? X_AA_NP / X_pp_NP : 0;
    double raaPR_err = raaPR * std::sqrt(std::pow(rel(yAAPR_err, yAAPR), 2) + std::pow(rel(yPpPR_err, yPpPR), 2));
    double raaNP_err = raaNP * std::sqrt(std::pow(rel(yAANP_err, yAANP), 2) + std::pow(rel(yPpNP_err, yPpNP), 2));
    out.raaPR.push_back(raaPR); out.raaPR_err.push_back(raaPR_err);
    out.raaNP.push_back(raaNP); out.raaNP_err.push_back(raaNP_err);
    out.npart.push_back(npart_central_first[c]);
  }
  return out;
}

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
    o.raa.push_back(h->GetBinContent(i+1)); o.err.push_back(he->GetBinContent(i+1));
  }
  return o;
}

TGraphErrors *makeGraph(const std::vector<double> &x, const std::vector<double> &y,
                        const std::vector<double> &ey, int color, int marker, double size = 1.4) {
  std::vector<double> ex(x.size(), 0.0);
  TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetMarkerColor(color); g->SetLineColor(color);
  g->SetMarkerStyle(marker); g->SetMarkerSize(size); g->SetLineWidth(2);
  return g;
}
}  // namespace

void compare_NP_G_vs_H_260520()
{
  gStyle->SetOptStat(0);
  setTDRStyle();
  int iPeriod = 101, iPos = 33;

  std::cout << "\n== Computing G_aggregate (2-exp NP) ==\n";
  RAAOut G_mid = computeNpart(true, kGDir, kGTag);
  RAAOut G_fwd = computeNpart(false, kGDir, kGTag);
  std::cout << "== Computing H_NP_rational (rational NP) ==\n";
  RAAOut H_mid = computeNpart(true, kHDir, kHTag);
  RAAOut H_fwd = computeNpart(false, kHDir, kHTag);

  const std::vector<double> np_mid_HIN = {358.8, 264.2, 189.2, 131.4, 86.9, 21.9};
  const std::vector<double> np_fwdPR_HIN = {311.5, 160.3, 32.7};
  HINNpart hin_prMid = LoadHIN("RAA_PR_JPsi_HIN_16_025_mid_Npart.root", "Table 15/Hist1D_y1", "Table 15/Hist1D_y1_e1", np_mid_HIN);
  HINNpart hin_prFwd = LoadHIN("RAA_PR_JPsi_HIN_16_025_fwd_Npart.root", "Table 19/Hist1D_y1", "Table 19/Hist1D_y1_e1", np_fwdPR_HIN);
  HINNpart hin_npMid = LoadHIN("RAA_NP_JPsi_HIN_16_025_Npart.root",     "Table 30/Hist1D_y1", "Table 30/Hist1D_y1_e1", np_mid_HIN);
  HINNpart hin_npFwd_hi = LoadHIN("RAA_NP_JPsi_HIN_16_025_Npart.root",  "Table 30/Hist1D_y2", "Table 30/Hist1D_y2_e1", np_mid_HIN);
  HINNpart hin_npFwd_lo = LoadHIN("RAA_NP_JPsi_HIN_16_025_fwd_Npart.root", "Table 32/Hist1D_y1", "Table 32/Hist1D_y1_e1", np_mid_HIN);

  // ----- Quick tabular comparison -----
  std::cout << "\n=== Mid Npart RAA: G vs H vs HIN-16-025 ===\n";
  std::cout << " <Npart>    |    PR (G/H/HIN)         |    NP (G/H/HIN)\n";
  for (size_t i = 0; i < G_mid.npart.size(); ++i) {
    double hPR = 0, hNP = 0;
    for (size_t j = 0; j < hin_prMid.npart.size(); ++j)
      if (std::abs(hin_prMid.npart[j] - G_mid.npart[i]) < 30) hPR = hin_prMid.raa[j];
    for (size_t j = 0; j < hin_npMid.npart.size(); ++j)
      if (std::abs(hin_npMid.npart[j] - G_mid.npart[i]) < 30) hNP = hin_npMid.raa[j];
    printf("  Np=%6.1f | %.3f / %.3f / %.3f | %.3f / %.3f / %.3f\n",
           G_mid.npart[i], G_mid.raaPR[i], H_mid.raaPR[i], hPR,
                           G_mid.raaNP[i], H_mid.raaNP[i], hNP);
  }
  std::cout << "\n=== Fwd Npart RAA: G vs H vs HIN-16-025 (NP fwd shown vs lopt/hipt) ===\n";
  for (size_t i = 0; i < G_fwd.npart.size(); ++i) {
    double hPR = 0, hNP_hi = 0, hNP_lo = 0;
    for (size_t j = 0; j < hin_prFwd.npart.size(); ++j)
      if (std::abs(hin_prFwd.npart[j] - G_fwd.npart[i]) < 50) hPR = hin_prFwd.raa[j];
    for (size_t j = 0; j < hin_npFwd_hi.npart.size(); ++j)
      if (std::abs(hin_npFwd_hi.npart[j] - G_fwd.npart[i]) < 30) hNP_hi = hin_npFwd_hi.raa[j];
    for (size_t j = 0; j < hin_npFwd_lo.npart.size(); ++j)
      if (std::abs(hin_npFwd_lo.npart[j] - G_fwd.npart[i]) < 30) hNP_lo = hin_npFwd_lo.raa[j];
    printf("  Np=%6.1f | PR G=%.3f H=%.3f HIN=%.3f | NP G=%.3f H=%.3f HIN_hi=%.3f HIN_lo=%.3f\n",
           G_fwd.npart[i], G_fwd.raaPR[i], H_fwd.raaPR[i], hPR,
                           G_fwd.raaNP[i], H_fwd.raaNP[i], hNP_hi, hNP_lo);
  }

  // ----- Draw one PDF per (component, rapidity) -----
  auto drawOne = [&](const char *outBase, const char *compTxt, const char *rapTxt,
                     const char *ptTxt, const RAAOut &G, const RAAOut &H,
                     const HINNpart &hin1, const char *hin1Lbl,
                     const HINNpart *hin2, const char *hin2Lbl, bool isPR)
  {
    TCanvas *c = new TCanvas(outBase, outBase, 700, 700);
    c->cd(); c->SetTicks(1, 1);
    c->DrawFrame(0, 0, 400, 1.5, ";#LTN_{part}#GT;R_{AA}");
    TLine *one = new TLine(0, 1, 400, 1);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw("same");

    TGraphErrors *gG = makeGraph(G.npart,
                                  isPR ? G.raaPR : G.raaNP,
                                  isPR ? G.raaPR_err : G.raaNP_err,
                                  kBlue + 2, 20, 1.4);
    TGraphErrors *gH = makeGraph(H.npart,
                                  isPR ? H.raaPR : H.raaNP,
                                  isPR ? H.raaPR_err : H.raaNP_err,
                                  kRed + 1, 21, 1.4);
    TGraphErrors *gHin1 = makeGraph(hin1.npart, hin1.raa, hin1.err, kGray+3, 28, 1.2);
    TGraphErrors *gHin2 = hin2 ? makeGraph(hin2->npart, hin2->raa, hin2->err, kGray+1, 25, 1.2) : nullptr;
    if (gHin2) gHin2->Draw("p same");
    gHin1->Draw("p same");
    gG->Draw("p same");
    gH->Draw("p same");

    TLegend *leg = new TLegend(0.40, 0.62, 0.88, 0.78);
    leg->SetTextSize(0.028); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(gG, "G (2-exp NP)", "pe");
    leg->AddEntry(gH, "H (rational NP)", "pe");
    leg->AddEntry(gHin1, hin1Lbl, "pe");
    if (gHin2) leg->AddEntry(gHin2, hin2Lbl, "pe");
    leg->Draw("same");

    drawText(compTxt, 0.21, 0.86, 1, 22);
    drawText(rapTxt,  0.21, 0.86 - 0.055, 1, 22);
    drawText(ptTxt,   0.21, 0.86 - 0.110, 1, 22);
    CMS_lumi_v2mass(c, iPeriod, iPos);

    c->SaveAs(Form("figs/%s.pdf", outBase));
    c->SaveAs(Form("figs/%s.png", outBase));
    std::cout << "  [wrote] figs/" << outBase << ".pdf\n";
  };

  std::cout << "\n";
  drawOne("compare_NP_G_vs_H_260520_PR_mid", "Prompt J/#psi", "|y| < 1.6", "6.5 < p_{T} < 40 GeV/c",
          G_mid, H_mid, hin_prMid, "HIN-16-025, 6.5 < p_{T} < 30", nullptr, nullptr, true);
  drawOne("compare_NP_G_vs_H_260520_PR_fwd", "Prompt J/#psi", "1.6 < |y| < 2.4", "3.5 < p_{T} < 40 GeV/c",
          G_fwd, H_fwd, hin_prFwd, "HIN-16-025, 3 < p_{T} < 30, 1.8<|y|<2.4", nullptr, nullptr, true);
  drawOne("compare_NP_G_vs_H_260520_NP_mid", "Nonprompt J/#psi", "|y| < 1.6", "6.5 < p_{T} < 40 GeV/c",
          G_mid, H_mid, hin_npMid, "HIN-16-025, 6.5<p_{T}<50, |y|<0.6", nullptr, nullptr, false);
  drawOne("compare_NP_G_vs_H_260520_NP_fwd", "Nonprompt J/#psi", "1.6 < |y| < 2.4", "3.5 < p_{T} < 40 GeV/c",
          G_fwd, H_fwd, hin_npFwd_hi, "HIN-16-025, 6.5<p_{T}<50, 1.8<|y|<2.4",
          &hin_npFwd_lo, "HIN-16-025, 3<p_{T}<6.5, 1.8<|y|<2.4", false);
}
