#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
struct SourceSet {
  TString tag;
  TString promptMid;
  TString promptFwd;
  TString nonpromptMid;
  TString nonpromptFwd;
};

struct ChannelInput {
  TString stateTag;
  TString rapTag;
  TString fileName;
};

TString BaseDir()
{
  return gSystem->DirName(__FILE__);
}

TString SrcPath(const TString &fileName)
{
  if (fileName.BeginsWith("/"))
    return fileName;
  return Form("%s/%s", BaseDir().Data(), fileName.Data());
}

TString OutFileName(const TString &outDir, const TString &stateTag, const TString &rapTag)
{
  const TString particle = (stateTag == "prompt") ? "Jpsi" : "BtoJpsi";
  return Form("%s/ratioDataMC_AA_%s_DATA_%s.root", outDir.Data(), particle.Data(), rapTag.Data());
}

TH1D *CloneWeightHist(TFile *file, const char *name)
{
  if (!file)
    return nullptr;
  TH1D *hist = dynamic_cast<TH1D *>(file->Get("WeightFactor"));
  if (!hist)
    return nullptr;
  TH1D *out = static_cast<TH1D *>(hist->Clone(name));
  out->SetDirectory(nullptr);
  return out;
}

TF1 *CloneSourceFunc(TFile *file, const char *name)
{
  if (!file)
    return nullptr;
  TF1 *func = dynamic_cast<TF1 *>(file->Get("dataMC_Ratio1"));
  if (!func)
    func = dynamic_cast<TF1 *>(file->Get("fitRatio1"));
  if (!func)
    return nullptr;
  return static_cast<TF1 *>(func->Clone(name));
}

TF1 *FitDoubleExp(TH1D *hist, const char *name)
{
  if (!hist)
    return nullptr;
  TF1 *func = new TF1(name, "[0]*TMath::Exp(-[1]*x)+[2]*TMath::Exp(-[3]*x)+[4]", 3.0, 40.0);
  func->SetParameters(1.0, 0.30, 0.5, 0.08, 0.0);
  func->SetParLimits(1, 0.0, 5.0);
  func->SetParLimits(3, 0.0, 5.0);
  hist->Fit(func, "Q0", "", 3.0, 40.0);
  return func;
}

TF1 *FitPoly3Pole(TH1D *hist, const char *name)
{
  if (!hist)
    return nullptr;
  TF1 *func = new TF1(name, "([0]+[1]*x+[2]*x*x+[4]*x*x*x)/TMath::Max(1e-3,(x-[3])*(x-[3])*(x-[3]))", 3.0, 40.0);
  func->SetParameters(1.0, 0.0, 0.0, 1.0, 0.0);
  func->SetParLimits(3, 0.5, 2.5);
  hist->Fit(func, "Q0", "", 3.0, 40.0);
  return func;
}

TF1 *BuildSelectedFunction(const TString &mode, TFile *sourceFile, TH1D *weightHist)
{
  if (mode == "sourceFunc")
  {
    TF1 *func = CloneSourceFunc(sourceFile, "dataMC_Ratio1");
    if (func)
      return func;
    return FitDoubleExp(weightHist, "dataMC_Ratio1");
  }
  if (mode == "doubleExp")
    return FitDoubleExp(weightHist, "dataMC_Ratio1");
  if (mode == "poly3Pole")
    return FitPoly3Pole(weightHist, "dataMC_Ratio1");

  std::cout << "[WARN] unknown mode " << mode << ", falling back to sourceFunc\n";
  return CloneSourceFunc(sourceFile, "dataMC_Ratio1");
}

bool WriteOneCandidate(const TString &mode, const TString &sourceTag, const ChannelInput &input,
                       std::ofstream &manifest)
{
  const TString srcPath = SrcPath(input.fileName);
  TFile *src = TFile::Open(srcPath, "READ");
  if (!src || src->IsZombie())
  {
    std::cout << "[WARN] missing source: " << srcPath << "\n";
    if (src)
      delete src;
    return false;
  }

  TH1D *weight = CloneWeightHist(src, "WeightFactor");
  if (!weight)
  {
    std::cout << "[WARN] missing WeightFactor in: " << srcPath << "\n";
    src->Close();
    delete src;
    return false;
  }

  const TString outDir = Form("%s/ptw_candidates_260508/%s_%s", BaseDir().Data(), sourceTag.Data(), mode.Data());
  gSystem->mkdir(outDir, true);
  const TString outPath = OutFileName(outDir, input.stateTag, input.rapTag);
  TFile *out = TFile::Open(outPath, "RECREATE");
  if (!out || out->IsZombie())
  {
    std::cout << "[ERROR] cannot write: " << outPath << "\n";
    src->Close();
    delete src;
    delete weight;
    return false;
  }

  TF1 *selected = BuildSelectedFunction(mode, src, weight);
  TF1 *sourceFunc = CloneSourceFunc(src, "source_dataMC_Ratio1");
  TF1 *doubleExp = FitDoubleExp(weight, "fit_doubleExp");
  TF1 *poly3Pole = FitPoly3Pole(weight, "fit_poly3Pole");

  out->cd();
  weight->Write("WeightFactor");
  if (selected)
    selected->Write("dataMC_Ratio1");
  if (sourceFunc)
    sourceFunc->Write("source_dataMC_Ratio1");
  if (doubleExp)
    doubleExp->Write("fit_doubleExp");
  if (poly3Pole)
    poly3Pole->Write("fit_poly3Pole");
  out->Close();

  manifest << sourceTag << "," << mode << "," << input.stateTag << "," << input.rapTag << ","
           << srcPath << "," << outPath << "\n";

  delete out;
  delete weight;
  delete selected;
  delete sourceFunc;
  delete doubleExp;
  delete poly3Pole;
  src->Close();
  delete src;
  return true;
}

std::vector<SourceSet> SourceSets()
{
  return {
      {"local260505_2exp",
       "ratioDataMC_AA_Jpsi_DATA_y0_1p6_260505_2exp.root",
       "ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_260505_2exp.root",
       "ratioDataMC_AA_BtoJpsi_DATA_y0_1p6_260505_2exp.root",
       "ratioDataMC_AA_BtoJpsi_DATA_y1p6_2p4_260505_2exp.root"},
      {"ctau260310_2exp_y024",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root"},
      {"split251118",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251118.root",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251118.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root"},
      {"split251103",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251103.root",
       "ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root",
       "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root"}};
}
}

void make_JPsi_pt_weight_candidates_260508()
{
  const TString outBase = Form("%s/ptw_candidates_260508", BaseDir().Data());
  gSystem->mkdir(outBase, true);

  std::ofstream manifest(Form("%s/manifest.csv", outBase.Data()));
  manifest << "source_tag,mode,state,rapidity,source_file,output_file\n";

  const std::vector<TString> modes = {"sourceFunc", "doubleExp", "poly3Pole"};
  int nWritten = 0;

  for (const auto &source : SourceSets())
  {
    const std::vector<ChannelInput> inputs = {
        {"prompt", "mid", source.promptMid},
        {"prompt", "fwd", source.promptFwd},
        {"nonprompt", "mid", source.nonpromptMid},
        {"nonprompt", "fwd", source.nonpromptFwd}};

    for (const auto &mode : modes)
    {
      for (const auto &input : inputs)
      {
        if (WriteOneCandidate(mode, source.tag, input, manifest))
          ++nWritten;
      }
    }
  }

  std::cout << "[DONE] wrote " << nWritten << " candidate pT-weight files under "
            << outBase << "\n";
  std::cout << "[DONE] manifest: " << outBase << "/manifest.csv\n";
}
