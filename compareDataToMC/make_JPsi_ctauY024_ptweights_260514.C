#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <memory>

namespace JPsiCtauY024PtW260514
{
const TString kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
const TString kOutDir = kRepo + "/compareDataToMC/ctau_y024_ptw_260514";

bool CopyAllObjects(const TString &srcPath, const TString &outPath,
                    const char *noteText)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie()) {
    std::cerr << "[make_JPsi_ctauY024_ptweights_260514] cannot open "
              << srcPath << std::endl;
    return false;
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "[make_JPsi_ctauY024_ptweights_260514] cannot create "
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
  TNamed note("ctau_y024_ptw_260514_note", noteText);
  TNamed method("ctau_y024_ptw_260514_method",
                "Single |y|<2.4 ctau-cut dN/dpT weight copied to both mid and fwd override names; no rapidity-split pT weights");
  note.Write();
  method.Write();
  out->Close();
  return true;
}
}

void make_JPsi_ctauY024_ptweights_260514()
{
  using namespace JPsiCtauY024PtW260514;
  if (gSystem->mkdir(kOutDir, true) != 0 && gSystem->AccessPathName(kOutDir)) {
    std::cerr << "[make_JPsi_ctauY024_ptweights_260514] cannot create "
              << kOutDir << std::endl;
    return;
  }

  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root",
                 "PbPb prompt: |y|<2.4 ctau-cut pT weight used for mid rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root",
                 "PbPb prompt: |y|<2.4 ctau-cut pT weight used for forward rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
                 "PbPb nonprompt: |y|<2.4 ctau-cut pT weight used for mid rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                 "PbPb nonprompt: |y|<2.4 ctau-cut pT weight used for forward rapidity");

  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_mid.root",
                 "pp prompt: |y|<2.4 ctau-cut pT weight used for mid rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_fwd.root",
                 "pp prompt: |y|<2.4 ctau-cut pT weight used for forward rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_mid.root",
                 "pp nonprompt: |y|<2.4 ctau-cut pT weight used for mid rapidity");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root",
                 kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_fwd.root",
                 "pp nonprompt: |y|<2.4 ctau-cut pT weight used for forward rapidity");

  std::cout << "[wrote ctau y0_2p4 pT weights] " << kOutDir << std::endl;
}
