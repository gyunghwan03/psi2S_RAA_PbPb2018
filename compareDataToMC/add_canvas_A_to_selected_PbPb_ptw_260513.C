#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TNamed.h"
#include "TPad.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <vector>

namespace AddPbPbCanvasA260513
{
struct Spec {
  TString label;
  TString source;
  TString target;
};

const TString kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
const TString kSelectedDir = kRepo + "/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear_npFwdLowBoostA095";

TH1 *CloneWeightFactor(TFile *file)
{
  if (!file)
    return nullptr;
  TH1 *hist = dynamic_cast<TH1 *>(file->Get("WeightFactor"));
  if (!hist)
    return nullptr;
  TH1 *clone = dynamic_cast<TH1 *>(hist->Clone("WeightFactor"));
  if (clone)
    clone->SetDirectory(nullptr);
  return clone;
}

TF1 *CloneSelectedFunction(TFile *file)
{
  if (!file)
    return nullptr;
  TF1 *func = dynamic_cast<TF1 *>(file->Get("dataMC_Ratio1"));
  if (!func)
    func = dynamic_cast<TF1 *>(file->Get("fitRatio1"));
  if (!func)
    return nullptr;
  TF1 *clone = dynamic_cast<TF1 *>(func->Clone("fitRatio1_canvas"));
  if (clone)
  {
    clone->SetLineColor(kRed + 1);
    clone->SetLineWidth(2);
  }
  return clone;
}

void StyleBottomHistogram(TH1 *hist)
{
  if (!hist)
    return;
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleSize(0.15);
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetLabelOffset(0.04);
  hist->GetXaxis()->SetLabelSize(0.15);
  hist->GetXaxis()->SetTickSize(0.03);
  hist->GetYaxis()->SetTickSize(0.04);
  hist->GetYaxis()->SetNdivisions(404);
  hist->GetYaxis()->SetTitle("Data/MC");
  hist->GetYaxis()->SetTitleOffset(0.25);
  hist->GetYaxis()->SetTitleSize(0.15);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetLabelSize(0.15);
  hist->GetYaxis()->SetNdivisions(404);
}

TCanvas *BuildCanvas(const Spec &spec, TFile *targetFile)
{
  std::unique_ptr<TFile> sourceFile(TFile::Open(spec.source, "READ"));
  if (!sourceFile || sourceFile->IsZombie())
  {
    std::cerr << "[add_canvas_A_to_selected_PbPb_ptw_260513] cannot open source "
              << spec.source << "\n";
    return nullptr;
  }

  TCanvas *sourceCanvas = dynamic_cast<TCanvas *>(sourceFile->Get("canvas_A"));
  if (!sourceCanvas)
  {
    std::cerr << "[add_canvas_A_to_selected_PbPb_ptw_260513] source has no canvas_A "
              << spec.source << "\n";
    return nullptr;
  }

  std::unique_ptr<TH1> ratio(CloneWeightFactor(targetFile));
  std::unique_ptr<TF1> curve(CloneSelectedFunction(targetFile));
  if (!ratio || !curve)
  {
    std::cerr << "[add_canvas_A_to_selected_PbPb_ptw_260513] target missing WeightFactor/dataMC_Ratio1 "
              << spec.target << "\n";
    return nullptr;
  }

  std::unique_ptr<TCanvas> canvas(dynamic_cast<TCanvas *>(sourceCanvas->Clone("canvas_A")));
  if (!canvas)
    return nullptr;

  TPad *padBottom = dynamic_cast<TPad *>(canvas->FindObject("pad_A_2"));
  if (!padBottom)
    padBottom = dynamic_cast<TPad *>(canvas->GetPrimitive("pad_A_2"));
  if (!padBottom)
  {
    std::cerr << "[add_canvas_A_to_selected_PbPb_ptw_260513] source canvas_A has no pad_A_2 "
              << spec.source << "\n";
    return nullptr;
  }

  padBottom->cd();
  padBottom->Clear();
  padBottom->SetFillColor(0);
  padBottom->SetBorderMode(0);
  padBottom->SetBorderSize(2);
  padBottom->SetTicks(1, 1);
  padBottom->SetBottomMargin(0.4361001);

  ratio->SetName("WeightFactor");
  StyleBottomHistogram(ratio.get());
  ratio->Draw();
  curve->Draw("same");

  TLine *unity = new TLine(ratio->GetXaxis()->GetXmin(), 1.0,
                           ratio->GetXaxis()->GetXmax(), 1.0);
  unity->SetLineStyle(2);
  unity->Draw("same");

  canvas->cd();
  canvas->Modified();
  canvas->Update();
  ratio.release();
  curve.release();
  return canvas.release();
}

void WriteOne(const Spec &spec)
{
  std::unique_ptr<TFile> targetFile(TFile::Open(spec.target, "UPDATE"));
  if (!targetFile || targetFile->IsZombie())
  {
    std::cerr << "[add_canvas_A_to_selected_PbPb_ptw_260513] cannot open target "
              << spec.target << "\n";
    return;
  }

  std::unique_ptr<TCanvas> canvas(BuildCanvas(spec, targetFile.get()));
  if (!canvas)
  {
    TNamed status("canvas_A_status", "failed to build PbPb reference-style canvas_A");
    targetFile->cd();
    status.Write("canvas_A_status", TObject::kOverwrite);
    return;
  }

  targetFile->cd();
  canvas->Write("canvas_A", TObject::kOverwrite);
  TNamed source("canvas_A_source", spec.source.Data());
  TNamed status("canvas_A_status",
                "reference-style PbPb dN/dpT + selected Data/MC canvas cloned from source canvas_A");
  source.Write("canvas_A_source", TObject::kOverwrite);
  status.Write("canvas_A_status", TObject::kOverwrite);
  std::cout << "[wrote canvas_A] " << spec.target << "\n";
}

void Run()
{
  const std::vector<Spec> specs = {
      {"prompt mid",
       kRepo + "/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y0_1p6_260505_2exp.root",
       kSelectedDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root"},
      {"prompt fwd",
       kRepo + "/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_260505_2exp.root",
       kSelectedDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root"},
      {"nonprompt mid",
       kRepo + "/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_y0_1p6_260505_2exp.root",
       kSelectedDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root"},
      {"nonprompt fwd",
       kRepo + "/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_y1p6_2p4_260505_2exp.root",
       kSelectedDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root"}};

  for (const Spec &spec : specs)
    WriteOne(spec);
}
}

void add_canvas_A_to_selected_PbPb_ptw_260513()
{
  AddPbPbCanvasA260513::Run();
}
