#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <memory>

namespace JPsiNoCtauSplitPp260514PtW
{
const TString kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
const TString kOutDir = kRepo + "/compareDataToMC/noctau_split_pp260514_ptw_260514";

bool CopyAllObjects(const TString &srcPath, const TString &outPath,
                    const char *noteText)
{
  std::unique_ptr<TFile> src(TFile::Open(srcPath, "READ"));
  if (!src || src->IsZombie()) {
    std::cerr << "[make_JPsi_noCtauSplit_pp260514_ptweights_260514] cannot open "
              << srcPath << std::endl;
    return false;
  }

  std::unique_ptr<TFile> out(TFile::Open(outPath, "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "[make_JPsi_noCtauSplit_pp260514_ptweights_260514] cannot create "
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
  TNamed note("noctau_split_pp260514_ptw_note", noteText);
  TNamed method("noctau_split_pp260514_ptw_method",
                "Rapidity-split no-ctau pT weights. PbPb files are kept from the previous noCtauSplit set; pp files use the newly produced 260514 noCtau weights.");
  note.Write();
  method.Write();
  out->Close();
  return true;
}
}

void make_JPsi_noCtauSplit_pp260514_ptweights_260514()
{
  using namespace JPsiNoCtauSplitPp260514PtW;
  if (gSystem->mkdir(kOutDir, true) != 0 && gSystem->AccessPathName(kOutDir)) {
    std::cerr << "[make_JPsi_noCtauSplit_pp260514_ptweights_260514] cannot create "
              << kOutDir << std::endl;
    return;
  }

  CopyAllObjects(kRepo + "/compareDataToMC/noctau_split_ptw_260514/ratioDataMC_AA_Jpsi_DATA_mid.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_mid.root",
                 "PbPb prompt: unchanged no-ctau pT weight for |y|<1.6.");
  CopyAllObjects(kRepo + "/compareDataToMC/noctau_split_ptw_260514/ratioDataMC_AA_Jpsi_DATA_fwd.root",
                 kOutDir + "/ratioDataMC_AA_Jpsi_DATA_fwd.root",
                 "PbPb prompt: unchanged no-ctau pT weight for 1.6<|y|<2.4.");
  CopyAllObjects(kRepo + "/compareDataToMC/noctau_split_ptw_260514/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_mid.root",
                 "PbPb nonprompt: unchanged no-ctau pT weight for |y|<1.6.");
  CopyAllObjects(kRepo + "/compareDataToMC/noctau_split_ptw_260514/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                 kOutDir + "/ratioDataMC_AA_BtoJpsi_DATA_fwd.root",
                 "PbPb nonprompt: unchanged no-ctau pT weight for 1.6<|y|<2.4.");

  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_noCtau_y0_1p6_260514.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_mid.root",
                 "pp prompt: new 260514 no-ctau pT weight for |y|<1.6.");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_noCtau_y1p6_2p4_260514.root",
                 kOutDir + "/ratioDataMC_pp_Jpsi_DATA_fwd.root",
                 "pp prompt: new 260514 no-ctau pT weight for 1.6<|y|<2.4.");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_noCtau_y0_1p6_260514.root",
                 kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_mid.root",
                 "pp nonprompt: new 260514 no-ctau pT weight for |y|<1.6.");
  CopyAllObjects(kRepo + "/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_noCtau_y1p6_2p4_260514.root",
                 kOutDir + "/ratioDataMC_pp_BtoJpsi_DATA_fwd.root",
                 "pp nonprompt: new 260514 no-ctau pT weight for 1.6<|y|<2.4.");

  std::cout << "[wrote no-ctau split pp260514 pT weights] " << kOutDir << std::endl;
}
