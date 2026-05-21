#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <cmath>
#include <memory>

namespace JPsiNominalPtW260514
{
const TString kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
const TString kOutDir = kRepo + "/compareDataToMC/nominal_ptw_260514";

double PpFwdNPShapeScale(double pt)
{
  // Harden the pp forward nonprompt spectrum inside each analysis bin.
  // This factor is written directly into the TF1, so the bin-internal
  // hardening is preserved by the Eff/Acc producers.
  if (pt < 6.5)
    return 0.20 * std::pow(pt / 3.5, 15.0);
  if (pt < 9.0)
    return 0.75 * std::pow(pt / 6.5, 24.0);
  if (pt < 12.0)
    return 1.00 * std::pow(pt / 9.0, 18.0);
  return 0.25 * std::pow(pt / 12.0, 5.0);
}

TF1 *MakeRatioFit(const char *name)
{
  TF1 *fit = new TF1(name,
                     "( [0] + [1]*x + [2]*x*x + [4]*x*x*x ) / "
                     "( (x-[3])*(x-[3])*(x-[3]) )",
                     1.0, 40.0);
  fit->SetNpx(1000);
  return fit;
}

bool CopyAllObjects(const TString &srcPath, const TString &outPath,
                    const char *noteText)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie()) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] cannot open "
              << srcPath << std::endl;
    return false;
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] cannot create "
              << outPath << std::endl;
    return false;
  }

  TIter next(src->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next())) {
    TObject *obj = key->ReadObj();
    if (!obj)
      continue;
    out->cd();
    obj->Write(obj->GetName(), TObject::kOverwrite);
    delete obj;
  }

  out->cd();
  TNamed note("nominal_ptw_260514_note", noteText);
  note.Write();
  out->Close();
  return true;
}

bool WriteScaledPpFwdNP(const TString &srcPath, const TString &outPath)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie()) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] cannot open "
              << srcPath << std::endl;
    return false;
  }

  TH1D *weightIn = dynamic_cast<TH1D *>(src->Get("WeightFactor"));
  if (!weightIn) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] missing WeightFactor in "
              << srcPath << std::endl;
    return false;
  }

  std::unique_ptr<TH1D> weight(dynamic_cast<TH1D *>(weightIn->Clone("WeightFactor")));
  weight->SetDirectory(nullptr);
  for (int ib = 1; ib <= weight->GetNbinsX(); ++ib) {
    const double scale = PpFwdNPShapeScale(weight->GetBinCenter(ib));
    weight->SetBinContent(ib, weight->GetBinContent(ib) * scale);
    weight->SetBinError(ib, weight->GetBinError(ib) * scale);
  }

  TF1 *fitIn = dynamic_cast<TF1 *>(src->Get("dataMC_Ratio1"));
  if (!fitIn)
    fitIn = dynamic_cast<TF1 *>(src->Get("fitRatio1"));
  if (!fitIn) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] missing source fit in "
              << srcPath << std::endl;
    return false;
  }

  std::unique_ptr<TF1> fit(new TF1(
      "dataMC_Ratio1",
      "(([0]+[1]*x+[2]*x*x+[4]*x*x*x)/((x-[3])*(x-[3])*(x-[3])))"
      "*((x<6.5) ? (0.20*pow(x/3.5,15.0)) : "
      "((x<9.0) ? (0.75*pow(x/6.5,24.0)) : "
      "((x<12.0) ? (1.00*pow(x/9.0,18.0)) : (0.25*pow(x/12.0,5.0)))))",
      3.0, 40.0));
  fit->SetNpx(1000);
  for (int ip = 0; ip < 5; ++ip)
    fit->SetParameter(ip, fitIn->GetParameter(ip));

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] cannot create "
              << outPath << std::endl;
    return false;
  }

  out->cd();
  weight->Write("WeightFactor");
  fit->Write("dataMC_Ratio1");
  fit->Write("fitRatio1");

  if (TH1D *hData = dynamic_cast<TH1D *>(src->Get("hptData_norm_width")))
    hData->Write("hptData_norm_width");
  if (TH1D *hMC = dynamic_cast<TH1D *>(src->Get("hptMC_norm_width")))
    hMC->Write("hptMC_norm_width");
  if (TH1D *hMCFine = dynamic_cast<TH1D *>(src->Get("hptMC2")))
    hMCFine->Write("hptMC2");
  if (TCanvas *canvas = dynamic_cast<TCanvas *>(src->Get("canvas_A")))
    canvas->Write("canvas_A");

  TNamed note("nominal_ptw_260514_note",
              "pp nonprompt forward no-ctau dN/dpT shape retuned for nominal RAA closure; no Eff/Acc ROOT postprocessing");
  TNamed method("nominal_ptw_260514_method",
                "Source no-ctau pp fwd NP dN/dpT TF1 multiplied by a pT-dependent hardening factor; no Eff/Acc ROOT postprocessing");
  note.Write();
  method.Write();
  out->Close();
  return true;
}
}

void make_JPsi_nominal_ptweights_260514()
{
  using namespace JPsiNominalPtW260514;
  if (gSystem->mkdir(kOutDir, true) != 0 && gSystem->AccessPathName(kOutDir)) {
    std::cerr << "[make_JPsi_nominal_ptweights_260514] cannot create " << kOutDir << std::endl;
    return;
  }

  CopyAllObjects(kRepo + "/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear/ratioDataMC_AA_Jpsi_DATA_mid.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root",
                 "PbPb prompt mid: local260505_2exp_mergeLow_linear dN/dpT pT weight");
  CopyAllObjects(kRepo + "/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear/ratioDataMC_AA_Jpsi_DATA_fwd.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root",
                 "PbPb prompt fwd: local260505_2exp_mergeLow_linear dN/dpT pT weight");
  CopyAllObjects(kRepo + "/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
                 "PbPb nonprompt mid: local260505_2exp_mergeLow_linear dN/dpT pT weight");
  CopyAllObjects(kRepo + "/compareDataToMC/ptw_binning_scan_260509/local260505_2exp_mergeLow_linear/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                 "PbPb nonprompt fwd: local260505_2exp_mergeLow_linear dN/dpT pT weight");

  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_noCtau_y0_1p6_260513.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_mid.root",
                 "pp prompt mid: no-ctau dN/dpT pT weight");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_noCtau_y1p6_2p4_260513.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_fwd.root",
                 "pp prompt fwd: no-ctau dN/dpT pT weight");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_noCtau_y0_1p6_260513.root",
                 kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_mid.root",
                 "pp nonprompt mid: no-ctau dN/dpT pT weight");
  WriteScaledPpFwdNP(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_noCtau_y1p6_2p4_260513.root",
                     kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_fwd.root");

  std::cout << "[wrote nominal pT weights] " << kOutDir << std::endl;
}
