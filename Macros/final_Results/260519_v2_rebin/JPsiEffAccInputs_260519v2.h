// Helper for 260519_v2 (G_aggregateHighPt_2exp full-stat Eff/Acc).
#ifndef JPSI_EFF_ACC_INPUTS_260519V2_H
#define JPSI_EFF_ACC_INPUTS_260519V2_H

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

#include <cmath>
#include <iostream>

namespace JPsiEffAcc260519v2
{
// Switched to H_NP_rational nominal on 260520 to lift NP mid RAA toward
// HIN-16-025. PR PtW is unchanged (rational fit is NP-only; PR files in
// H_NP_rational_full are symlinks back to G_aggregate's PR ROOTs).
inline TString Dir() { return "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/H_NP_rational_full"; }
inline TString Tag() { return "260520_H_NP_rational"; }

inline TFile *Open(const char *fileName)
{
  TFile *f = TFile::Open(Form("%s/%s", Dir().Data(), fileName), "READ");
  if (!f || f->IsZombie()) std::cerr << "[260519v2] cannot open " << fileName << "\n";
  return f;
}

inline TH1D *BuildIntegratedRatioHist(TFile *file, const char *numName,
                                       const char *denName, const char *outName)
{
  TH1D *out = new TH1D(outName, "", 1, 0.5, 1.5);
  out->SetDirectory(0);
  if (!file || file->IsZombie()) return out;
  TH1 *num = dynamic_cast<TH1 *>(file->Get(numName));
  TH1 *den = dynamic_cast<TH1 *>(file->Get(denName));
  if (!num || !den) return out;
  double sN = 0.0, eN2 = 0.0, sD = 0.0, eD2 = 0.0;
  for (int ib = 1; ib <= num->GetNbinsX(); ++ib) { sN += num->GetBinContent(ib); eN2 += std::pow(num->GetBinError(ib), 2); }
  for (int ib = 1; ib <= den->GetNbinsX(); ++ib) { sD += den->GetBinContent(ib); eD2 += std::pow(den->GetBinError(ib), 2); }
  if (sD > 0.0) {
    double r = sN / sD;
    out->SetBinContent(1, r);
    out->SetBinError(1, r * std::sqrt(eN2/(sN*sN) + eD2/(sD*sD)));
  }
  return out;
}
}  // namespace
#endif
