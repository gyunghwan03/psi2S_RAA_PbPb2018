#ifndef JPSI_EFF_ACC_INPUTS_260429_H
#define JPSI_EFF_ACC_INPUTS_260429_H

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace JPsiEffAcc260429
{
inline TString RootDir()
{
  return "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260429/roots";
}

inline TFile *Open(const char *fileName)
{
  TFile *file = TFile::Open(Form("%s/%s", RootDir().Data(), fileName), "READ");
  if (!file || file->IsZombie())
    std::cerr << "[JPsiEffAcc260429] failed to open " << fileName << std::endl;
  return file;
}

inline TH1D *BuildIntegratedRatioHist(TFile *file, const char *numName, const char *denName, const char *outName)
{
  TH1D *out = new TH1D(outName, "", 1, 0.5, 1.5);
  out->SetDirectory(0);

  if (!file || file->IsZombie())
    return out;

  TH1 *num = dynamic_cast<TH1 *>(file->Get(numName));
  TH1 *den = dynamic_cast<TH1 *>(file->Get(denName));
  if (!num || !den)
  {
    std::cerr << "[JPsiEffAcc260429] missing " << numName << " or " << denName
              << " in " << file->GetName() << std::endl;
    return out;
  }

  double sumNum = 0.0;
  double err2Num = 0.0;
  double sumDen = 0.0;
  double err2Den = 0.0;
  for (int ib = 1; ib <= num->GetNbinsX(); ++ib)
  {
    sumNum += num->GetBinContent(ib);
    err2Num += std::pow(num->GetBinError(ib), 2);
  }
  for (int ib = 1; ib <= den->GetNbinsX(); ++ib)
  {
    sumDen += den->GetBinContent(ib);
    err2Den += std::pow(den->GetBinError(ib), 2);
  }

  if (sumDen <= 0.0)
    return out;

  const double value = sumNum / sumDen;
  double error = 0.0;
  if (sumNum > 0.0)
  {
    const double relNum = std::sqrt(err2Num) / sumNum;
    const double relDen = std::sqrt(err2Den) / sumDen;
    error = value * std::sqrt(relNum * relNum + relDen * relDen);
  }
  out->SetBinContent(1, value);
  out->SetBinError(1, error);
  return out;
}

inline TH1D *MapCentralityHistToPercentBins(const TH1 *source, int nBins, const double *centBinsPercent, const char *outName)
{
  TH1D *out = new TH1D(outName, "", nBins, centBinsPercent);
  out->SetDirectory(0);

  if (!source)
    return out;

  for (int ib = 1; ib <= nBins; ++ib)
  {
    const double targetLow = 2.0 * centBinsPercent[ib - 1];
    const double targetHigh = 2.0 * centBinsPercent[ib];
    double weight = 0.0;
    double value = 0.0;
    double err2 = 0.0;

    for (int jb = 1; jb <= source->GetNbinsX(); ++jb)
    {
      const double sourceLow = source->GetBinLowEdge(jb);
      const double sourceHigh = source->GetBinLowEdge(jb + 1);
      const double overlap = std::min(targetHigh, sourceHigh) - std::max(targetLow, sourceLow);
      if (overlap <= 0.0)
        continue;

      weight += overlap;
      value += overlap * source->GetBinContent(jb);
      err2 += std::pow(overlap * source->GetBinError(jb), 2);
    }

    if (weight <= 0.0)
      continue;

    out->SetBinContent(ib, value / weight);
    out->SetBinError(ib, std::sqrt(err2) / weight);
  }

  return out;
}
}

#endif
