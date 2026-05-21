#include "TFile.h"
#include "TF1.h"
#include "TKey.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <memory>

namespace JPsiNPFwdLowPtBoost260509
{
const TString kBaseDir = "/data/hwan/psi2S_RAA_PbPb2018/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear";

TString OutDir(double alpha)
{
  const int alphaTag = static_cast<int>(alpha * 100.0 + 0.5);
  return Form("/data/hwan/psi2S_RAA_PbPb2018/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear_npFwdLowBoostA%03d",
              alphaTag);
}

void CopyFile(const TString &srcPath, const TString &outPath)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie())
  {
    std::cerr << "[make_JPsi_npFwd_lowPtBoost_260509] failed to open " << srcPath << std::endl;
    return;
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie())
  {
    std::cerr << "[make_JPsi_npFwd_lowPtBoost_260509] failed to create " << outPath << std::endl;
    return;
  }

  TIter next(src->GetListOfKeys());
  TKey *key = nullptr;
  while ((key = static_cast<TKey *>(next())))
  {
    TObject *obj = key->ReadObj();
    if (!obj)
      continue;
    out->cd();
    obj->Write(obj->GetName(), TObject::kOverwrite);
  }
  std::cout << "[COPY] " << outPath << std::endl;
}

void MakeBoostedNPFwd(const TString &srcPath, const TString &outPath, double alpha)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie())
  {
    std::cerr << "[make_JPsi_npFwd_lowPtBoost_260509] failed to open " << srcPath << std::endl;
    return;
  }

  TF1 *base = dynamic_cast<TF1 *>(src->Get("dataMC_Ratio1"));
  if (!base)
  {
    std::cerr << "[make_JPsi_npFwd_lowPtBoost_260509] missing dataMC_Ratio1 in " << srcPath << std::endl;
    return;
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie())
  {
    std::cerr << "[make_JPsi_npFwd_lowPtBoost_260509] failed to create " << outPath << std::endl;
    return;
  }

  TIter next(src->GetListOfKeys());
  TKey *key = nullptr;
  while ((key = static_cast<TKey *>(next())))
  {
    TObject *obj = key->ReadObj();
    if (!obj || TString(obj->GetName()) == "dataMC_Ratio1")
      continue;
    out->cd();
    obj->Write(obj->GetName(), TObject::kOverwrite);
  }

  TF1 boosted("dataMC_Ratio1",
              "([0]+[1]*x)*((x<6.5)*exp([2]*(6.5-x))+(x>=6.5))",
              3.5, 40.0);
  boosted.SetParameter(0, base->GetParameter(0));
  boosted.SetParameter(1, base->GetParameter(1));
  boosted.SetParameter(2, alpha);
  boosted.Write("dataMC_Ratio1", TObject::kOverwrite);

  TF1 source("source_dataMC_Ratio1", "[0]+[1]*x", 3.5, 40.0);
  source.SetParameter(0, base->GetParameter(0));
  source.SetParameter(1, base->GetParameter(1));
  source.Write("source_dataMC_Ratio1", TObject::kOverwrite);

  std::cout << "[BOOST] " << outPath << " alpha=" << alpha
            << " base p0=" << base->GetParameter(0)
            << " p1=" << base->GetParameter(1) << std::endl;
}

void Run(double alpha)
{
  const TString outDir = OutDir(alpha);
  gSystem->mkdir(outDir, true);
  CopyFile(kBaseDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root",
           outDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root");
  CopyFile(kBaseDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root",
           outDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root");
  CopyFile(kBaseDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
           outDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root");
  MakeBoostedNPFwd(kBaseDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                   outDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root", alpha);
}
}

void make_JPsi_npFwd_lowPtBoost_260509(double alpha = 0.20)
{
  JPsiNPFwdLowPtBoost260509::Run(alpha);
}
