#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
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

struct BinningScheme {
  TString tag;
  std::vector<int> edgeBins;
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

TH1D *GetWeightHist(TFile *file, const char *name)
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

std::vector<int> UniqueSortedEdges(const std::vector<int> &edges, int nBins)
{
  std::set<int> keep;
  for (int edge : edges)
  {
    if (edge < 1)
      edge = 1;
    if (edge > nBins + 1)
      edge = nBins + 1;
    keep.insert(edge);
  }
  keep.insert(1);
  keep.insert(nBins + 1);
  return std::vector<int>(keep.begin(), keep.end());
}

std::vector<BinningScheme> BuildSchemes(int nBins)
{
  std::vector<int> nominal;
  for (int edge = 1; edge <= nBins + 1; ++edge)
    nominal.push_back(edge);

  std::vector<BinningScheme> out;
  out.push_back({"nominal", nominal});

  std::vector<int> mergeLow = nominal;
  if (nBins >= 3)
    mergeLow.erase(mergeLow.begin() + 1);
  out.push_back({"mergeLow", UniqueSortedEdges(mergeLow, nBins)});

  std::vector<int> mergeHigh = nominal;
  if (nBins >= 3)
    mergeHigh.erase(mergeHigh.end() - 2);
  out.push_back({"mergeHigh", UniqueSortedEdges(mergeHigh, nBins)});

  std::vector<int> mergeBoth = nominal;
  if (nBins >= 5)
  {
    mergeBoth.erase(mergeBoth.end() - 2);
    mergeBoth.erase(mergeBoth.begin() + 1);
  }
  out.push_back({"mergeLowHigh", UniqueSortedEdges(mergeBoth, nBins)});

  std::vector<int> everyOther;
  for (int edge = 1; edge <= nBins + 1; edge += 2)
    everyOther.push_back(edge);
  everyOther.push_back(nBins + 1);
  out.push_back({"everyOther", UniqueSortedEdges(everyOther, nBins)});

  std::vector<int> coarse3 = {1, std::max(2, nBins / 2), nBins + 1};
  out.push_back({"coarse3", UniqueSortedEdges(coarse3, nBins)});

  return out;
}

TH1D *MergeWeightHist(const TH1D *src, const BinningScheme &scheme, const char *name)
{
  if (!src)
    return nullptr;

  std::vector<int> edgeBins = UniqueSortedEdges(scheme.edgeBins, src->GetNbinsX());
  std::vector<double> edges;
  for (int edgeBin : edgeBins)
    edges.push_back(src->GetXaxis()->GetBinLowEdge(edgeBin));
  edges.back() = src->GetXaxis()->GetBinUpEdge(src->GetNbinsX());

  TH1D *out = new TH1D(name, src->GetTitle(), static_cast<int>(edges.size()) - 1, edges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();

  for (int ib = 1; ib <= out->GetNbinsX(); ++ib)
  {
    const int first = edgeBins[ib - 1];
    const int last = edgeBins[ib] - 1;
    double sumW = 0.0;
    double sumVW = 0.0;
    double sumWidth = 0.0;
    double sumVWidth = 0.0;

    for (int jb = first; jb <= last; ++jb)
    {
      const double value = src->GetBinContent(jb);
      const double err = src->GetBinError(jb);
      const double width = src->GetBinWidth(jb);
      sumWidth += width;
      sumVWidth += value * width;
      if (err > 0.0)
      {
        const double w = 1.0 / (err * err);
        sumW += w;
        sumVW += value * w;
      }
    }

    if (sumW > 0.0)
    {
      out->SetBinContent(ib, sumVW / sumW);
      out->SetBinError(ib, std::sqrt(1.0 / sumW));
    }
    else if (sumWidth > 0.0)
    {
      out->SetBinContent(ib, sumVWidth / sumWidth);
      out->SetBinError(ib, 0.0);
    }
  }
  return out;
}

TF1 *FitFunction(TH1D *hist, const TString &mode, const char *name)
{
  if (!hist)
    return nullptr;

  const double xMin = hist->GetXaxis()->GetXmin();
  const double xMax = hist->GetXaxis()->GetXmax();
  TF1 *func = nullptr;

  if (mode == "linear")
  {
    func = new TF1(name, "[0]+[1]*x", xMin, xMax);
    func->SetParameters(1.0, 0.0);
  }
  else if (mode == "poly2")
  {
    func = new TF1(name, "[0]+[1]*x+[2]*x*x", xMin, xMax);
    func->SetParameters(1.0, 0.0, 0.0);
  }
  else if (mode == "expoConst")
  {
    func = new TF1(name, "[0]*TMath::Exp(-[1]*x)+[2]", xMin, xMax);
    func->SetParameters(1.0, 0.1, 0.2);
    func->SetParLimits(1, 0.0, 5.0);
  }
  else if (mode == "doubleExp")
  {
    func = new TF1(name, "[0]*TMath::Exp(-[1]*x)+[2]*TMath::Exp(-[3]*x)+[4]", xMin, xMax);
    func->SetParameters(1.0, 0.30, 0.5, 0.08, 0.0);
    func->SetParLimits(1, 0.0, 5.0);
    func->SetParLimits(3, 0.0, 5.0);
  }

  if (!func)
    return nullptr;

  if (hist->GetNbinsX() <= func->GetNpar())
    hist->Fit(func, "Q0W", "", xMin, xMax);
  else
    hist->Fit(func, "Q0", "", xMin, xMax);
  return func;
}

bool WriteOneCandidate(const TString &sourceTag, const TString &fitMode,
                       const BinningScheme &scheme, const ChannelInput &input,
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

  TH1D *sourceWeight = GetWeightHist(src, "source_WeightFactor");
  if (!sourceWeight)
  {
    std::cout << "[WARN] missing WeightFactor in: " << srcPath << "\n";
    src->Close();
    delete src;
    return false;
  }

  TH1D *weight = MergeWeightHist(sourceWeight, scheme, "WeightFactor");
  if (!weight)
  {
    src->Close();
    delete src;
    delete sourceWeight;
    return false;
  }

  const TString caseTag = Form("%s_%s_%s", sourceTag.Data(), scheme.tag.Data(), fitMode.Data());
  const TString outDir = Form("%s/ptw_binning_scan_260509/%s", BaseDir().Data(), caseTag.Data());
  gSystem->mkdir(outDir, true);
  const TString outPath = OutFileName(outDir, input.stateTag, input.rapTag);
  TFile *out = TFile::Open(outPath, "RECREATE");
  if (!out || out->IsZombie())
  {
    std::cout << "[ERROR] cannot write: " << outPath << "\n";
    src->Close();
    delete src;
    delete sourceWeight;
    delete weight;
    return false;
  }

  TF1 *selected = FitFunction(weight, fitMode, "dataMC_Ratio1");
  TF1 *sourceFunc = CloneSourceFunc(src, "source_dataMC_Ratio1");

  out->cd();
  weight->Write("WeightFactor");
  sourceWeight->Write("source_WeightFactor");
  if (selected)
    selected->Write("dataMC_Ratio1");
  if (sourceFunc)
    sourceFunc->Write("source_dataMC_Ratio1");
  out->Close();

  manifest << caseTag << "," << sourceTag << "," << scheme.tag << "," << fitMode << ","
           << input.stateTag << "," << input.rapTag << "," << srcPath << ","
           << outPath << "," << weight->GetNbinsX() << "\n";

  delete out;
  delete sourceWeight;
  delete weight;
  delete selected;
  delete sourceFunc;
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
      {"local260503",
       "ratioDataMC_AA_Jpsi_DATA_y0_1p6_260503.root",
       "ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_260503.root",
       "ratioDataMC_AA_BtoJpsi_DATA_y0_1p6_260503.root",
       "ratioDataMC_AA_BtoJpsi_DATA_y1p6_2p4_260503.root"}};
}
}

void make_JPsi_pt_weight_binning_scan_260509()
{
  const TString outBase = Form("%s/ptw_binning_scan_260509", BaseDir().Data());
  gSystem->mkdir(outBase, true);

  std::ofstream manifest(Form("%s/manifest.csv", outBase.Data()));
  manifest << "case_tag,source_tag,binning_tag,fit_mode,state,rapidity,source_file,output_file,n_weight_bins\n";

  const std::vector<TString> fitModes = {"linear", "poly2", "expoConst", "doubleExp"};
  int nWritten = 0;

  for (const auto &source : SourceSets())
  {
    const std::vector<ChannelInput> inputs = {
        {"prompt", "mid", source.promptMid},
        {"prompt", "fwd", source.promptFwd},
        {"nonprompt", "mid", source.nonpromptMid},
        {"nonprompt", "fwd", source.nonpromptFwd}};

    for (const auto &input : inputs)
    {
      TFile *probe = TFile::Open(SrcPath(input.fileName), "READ");
      TH1D *probeHist = GetWeightHist(probe, "probe");
      if (!probeHist)
      {
        std::cout << "[WARN] cannot build schemes for " << input.fileName << "\n";
        if (probe)
        {
          probe->Close();
          delete probe;
        }
        continue;
      }
      const std::vector<BinningScheme> schemes = BuildSchemes(probeHist->GetNbinsX());
      delete probeHist;
      probe->Close();
      delete probe;

      for (const auto &scheme : schemes)
      {
        for (const auto &fitMode : fitModes)
        {
          if (WriteOneCandidate(source.tag, fitMode, scheme, input, manifest))
            ++nWritten;
        }
      }
    }
  }

  std::cout << "[DONE] wrote " << nWritten << " binned pT-weight files under "
            << outBase << "\n";
  std::cout << "[DONE] manifest: " << outBase << "/manifest.csv\n";
}
