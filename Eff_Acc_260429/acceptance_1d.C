#include "TBranch.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TLeaf.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using std::cout;

namespace
{
  TString resolveBaseDir()
  {
    return gSystem->DirName(__FILE__);
  }

  bool IsAcceptanceQQRun2018(double pt, double eta)
  {
    const double aeta = std::abs(eta);

    return ((aeta < 1.2 && pt >= 3.5) ||
            (1.2 <= aeta && aeta < 2.1 && pt >= 5.47 - 1.89 * aeta) ||
            (2.1 <= aeta && aeta < 2.4 && pt >= 1.5));
  }

  double GetPtWeight(TF1 *fW, double pt)
  {
    return fW ? fW->Eval(pt) : 1.0;
  }

  TFile *OpenFirstValid(const std::vector<std::string> &paths, const char *tag)
  {
    for (const auto &path : paths)
    {
      TFile *file = TFile::Open(path.c_str(), "READ");
      if (file && !file->IsZombie())
      {
        cout << "[PtW] " << tag << " : " << path << "\n";
        return file;
      }
      if (file)
      {
        file->Close();
        delete file;
      }
    }
    return nullptr;
  }

  TF1 *LoadClonedWeight(TFile *file, const std::vector<std::string> &sourceNames,
                        const char *cloneName, std::string &loadedName)
  {
    if (!file)
      return nullptr;

    for (const auto &sourceName : sourceNames)
    {
      TF1 *fIn = dynamic_cast<TF1 *>(file->Get(sourceName.c_str()));
      if (!fIn)
        continue;

      loadedName = sourceName;
      return static_cast<TF1 *>(fIn->Clone(cloneName));
    }
    return nullptr;
  }

  TF1 *LoadPtWeightFunction(const std::vector<std::string> &paths,
                            const std::vector<std::string> &sourceNames,
                            const char *tag, const char *cloneName)
  {
    TFile *file = OpenFirstValid(paths, tag);
    if (!file)
    {
      cout << "[ERROR] failed to open " << tag << " pT weight file.\n";
      return nullptr;
    }

    std::string loadedName;
    TF1 *weight = LoadClonedWeight(file, sourceNames, cloneName, loadedName);
    if (!weight)
    {
      cout << "[ERROR] failed to load " << tag << " pT weight function. Tried:";
      for (const auto &sourceName : sourceNames)
        cout << " " << sourceName;
      cout << "\n";
    }
    else
    {
      cout << "[PtW] " << tag << " function : " << loadedName << "\n";
    }

    file->Close();
    delete file;
    return weight;
  }

  bool BranchUsesShort(TTree *tree, const char *branchName)
  {
    if (!tree)
      return false;
    TBranch *branch = tree->GetBranch(branchName);
    if (!branch)
      return false;
    TLeaf *leaf = branch->GetLeaf(branchName);
    if (!leaf && branch->GetListOfLeaves() && branch->GetListOfLeaves()->GetEntries() > 0)
      leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->At(0));
    if (!leaf)
      return false;
    const std::string typeName = leaf->GetTypeName();
    return typeName == "Short_t" || typeName == "UShort_t";
  }

  TH1D *BuildAcceptanceHist(TH1D *num, TH1D *den, const char *name, const char *title)
  {
    if (!num || !den)
      return nullptr;

    TH1D *acc = static_cast<TH1D *>(num->Clone(name));
    acc->SetTitle(title);
    acc->Divide(num, den, 1.0, 1.0, "B");
    return acc;
  }
} // namespace

void acceptance_1d(long nEvt = 100000, bool isPr = true, bool isPbPb = false, bool isGenW = true, bool isPtW = true)
{
  cout << "Start acceptance_1d()\n";
  TH1::SetDefaultSumw2(kTRUE);

  TChain *oniaTree = new TChain("hionia/myTree");
  const std::string inputFile = isPr
                                    ? "/data/Oniatree/Jpsi/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root"
                                    : "/data/Oniatree/Jpsi/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  cout << "[INFO] input MC is fixed to pp source: " << inputFile << "\n";
  if (oniaTree->Add(inputFile.c_str()) <= 0)
  {
    cout << "[ERROR] failed to add input file: " << inputFile << "\n";
    return;
  }

  TF1 *fPtW_mid = nullptr;
  TF1 *fPtW_fwd = nullptr;
  if (isPtW)
  {
    const std::string hanseulPtWDir = "/data/users/hanseul/jpsi_2603/compareDataToMC/roots";
    const TString baseDir = resolveBaseDir();
    const TString compareDir = baseDir + "/../compareDataToMC";
    auto MakeHanseulPath = [&](const char *fileName) -> std::string
    {
      return std::string(Form("%s/%s", hanseulPtWDir.c_str(), fileName));
    };

    std::vector<std::string> ptWeightPathsMid;
    std::vector<std::string> ptWeightPathsFwd;
    const std::vector<std::string> functionNames = {"dataMC_Ratio1", "fitRatio1"};
    if (isPbPb)
    {
      if (isPr)
      {
        ptWeightPathsMid = {
            MakeHanseulPath("ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260410.root"),
            Form("%s/roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", baseDir.Data()),
            Form("%s/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", baseDir.Data()),
            Form("%s/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251118.root", compareDir.Data())};
        ptWeightPathsFwd = {
            MakeHanseulPath("ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_260410.root"),
            Form("%s/ratioDataMC_pbpb_250502_y1p6_2p4.root", compareDir.Data()),
            Form("%s/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root", compareDir.Data())};
      }
      else
      {
        ptWeightPathsMid = {
            MakeHanseulPath("ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260410.root"),
            Form("%s/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251118.root", compareDir.Data())};
        ptWeightPathsFwd = {
            MakeHanseulPath("ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_260410.root"),
            Form("%s/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root", compareDir.Data())};
      }
    }
    else
    {
      if (isPr)
      {
        ptWeightPathsMid = {
            MakeHanseulPath("ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260410.root"),
            Form("%s/ratioDataMC_pp_250219_y0_1p6.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251125.root", compareDir.Data())};
        ptWeightPathsFwd = {
            MakeHanseulPath("ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_260410.root"),
            MakeHanseulPath("ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_260410_bf.root"),
            Form("%s/ratioDataMC_pp_250219_y1p6_2p4.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251125.root", compareDir.Data())};
      }
      else
      {
        ptWeightPathsMid = {
            MakeHanseulPath("ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260410.root"),
            Form("%s/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251125.root", compareDir.Data())};
        ptWeightPathsFwd = {
            MakeHanseulPath("ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_260410.root"),
            Form("%s/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root", compareDir.Data()),
            Form("%s/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251125.root", compareDir.Data())};
      }
    }

    fPtW_mid = LoadPtWeightFunction(
        ptWeightPathsMid, functionNames,
        isPbPb ? (isPr ? "PbPb prompt mid" : "PbPb nonprompt mid")
               : (isPr ? "pp prompt mid" : "pp nonprompt mid"),
        "acc_run2_ptw_mid");
    fPtW_fwd = LoadPtWeightFunction(
        ptWeightPathsFwd, functionNames,
        isPbPb ? (isPr ? "PbPb prompt forward" : "PbPb nonprompt forward")
               : (isPr ? "pp prompt forward" : "pp nonprompt forward"),
        "acc_run2_ptw_fwd");
    if (!fPtW_mid || !fPtW_fwd)
    {
      cout << "[ERROR] failed to load rapidity-dependent pT weight functions.\n";
      return;
    }
  }

  const bool usesSplitMomentumBranches = (oniaTree->GetBranch("Gen_QQ_4mom_pt") != nullptr);
  const bool hasGenQQIndexBranches = (oniaTree->GetBranch("Gen_QQ_mupl_idx") != nullptr &&
                                      oniaTree->GetBranch("Gen_QQ_mumi_idx") != nullptr &&
                                      oniaTree->GetBranch("Gen_mu_charge") != nullptr);
  const bool useShortSizeBranches = BranchUsesShort(oniaTree, "Gen_QQ_size");
  const bool useShortIndexBranches = hasGenQQIndexBranches && BranchUsesShort(oniaTree, "Gen_QQ_mupl_idx");
  const bool hasPairedMuonMomentumBranches = (!usesSplitMomentumBranches &&
                                              oniaTree->GetBranch("Gen_QQ_mupl_4mom") != nullptr &&
                                              oniaTree->GetBranch("Gen_QQ_mumi_4mom") != nullptr);
  const bool hasGenWeightBranch = (oniaTree->GetBranch("Gen_weight") != nullptr);

  TClonesArray *Gen_QQ_4mom = nullptr;
  TClonesArray *Gen_mu_4mom = nullptr;
  TClonesArray *Gen_QQ_mupl_4mom = nullptr;
  TClonesArray *Gen_QQ_mumi_4mom = nullptr;
  std::vector<float> *Gen_QQ_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_pt = nullptr;
  std::vector<float> *Gen_QQ_4mom_phi = nullptr;
  std::vector<float> *Gen_QQ_4mom_eta = nullptr;
  std::vector<float> *Gen_mu_4mom_pt = nullptr;
  std::vector<float> *Gen_mu_4mom_phi = nullptr;
  std::vector<float> *Gen_mu_4mom_eta = nullptr;

  if (usesSplitMomentumBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_4mom_m", &Gen_QQ_4mom_m);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_pt", &Gen_QQ_4mom_pt);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_phi", &Gen_QQ_4mom_phi);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_eta", &Gen_QQ_4mom_eta);
    oniaTree->SetBranchAddress("Gen_mu_4mom_pt", &Gen_mu_4mom_pt);
    oniaTree->SetBranchAddress("Gen_mu_4mom_phi", &Gen_mu_4mom_phi);
    oniaTree->SetBranchAddress("Gen_mu_4mom_eta", &Gen_mu_4mom_eta);
  }
  else
  {
    oniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
    if (hasGenQQIndexBranches)
      oniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
    else if (hasPairedMuonMomentumBranches)
    {
      oniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);
      oniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
    }
    else
    {
      cout << "[ERROR] unsupported GEN-only tree schema: missing GEN muon index and paired muon momentum branches.\n";
      return;
    }
  }

  Short_t Gen_QQ_size_short = 0;
  Int_t Gen_QQ_size_int = 0;
  Short_t Gen_QQ_mupl_idx_short[1000];
  Short_t Gen_QQ_mumi_idx_short[1000];
  Short_t Gen_mu_charge_short[1000];
  Int_t Gen_QQ_mupl_idx_int[1000];
  Int_t Gen_QQ_mumi_idx_int[1000];
  Int_t Gen_mu_charge_int[1000];
  Float_t Gen_weight = 1.f;

  if (useShortSizeBranches)
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_short);
  else
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_int);

  if (hasGenQQIndexBranches && useShortIndexBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_short);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_short);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_short);
  }
  else if (hasGenQQIndexBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_int);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_int);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_int);
  }

  if (isGenW && hasGenWeightBranch)
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);

  auto GetGenQQSize = [&]() -> Int_t
  {
    return useShortSizeBranches ? static_cast<Int_t>(Gen_QQ_size_short) : Gen_QQ_size_int;
  };
  auto GetGenQQMuplIdx = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_QQ_mupl_idx_short[idx]) : Gen_QQ_mupl_idx_int[idx];
  };
  auto GetGenQQMumiIdx = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_QQ_mumi_idx_short[idx]) : Gen_QQ_mumi_idx_int[idx];
  };
  auto GetGenMuCharge = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_mu_charge_short[idx]) : Gen_mu_charge_int[idx];
  };

  const TString outDir = resolveBaseDir() + "/roots";
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] failed to create output directory: " << outDir << "\n";
    return;
  }

  const std::string legacyCollLabel = isPbPb ? "PbPb" : "pp";
  const std::string flavorLabel = isPr ? "PromptJpsi" : "BtoJpsi";
  const std::string outFilePath = Form("%s/acceptance_%s_GenOnly_wgt%d_%s_SysUp0_260429_1d.root",
                                       outDir.Data(), flavorLabel.c_str(), isPtW ? 1 : 0,
                                       legacyCollLabel.c_str());
  TFile *fout = TFile::Open(outFilePath.c_str(), "RECREATE");
  if (!fout || fout->IsZombie() || !fout->IsWritable())
  {
    cout << "[ERROR] cannot open output file for writing: " << outFilePath << "\n";
    return;
  }

  const double ptBinsAll[] = {6.5, 7.5, 8.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 50.0};
  const double ptBinsMid[] = {6.5, 9, 12, 15, 20, 25, 40};
  const double ptBinsFwd[] = {3.5, 6.5, 9, 12, 40};
  const int nPtBinsAll = sizeof(ptBinsAll) / sizeof(double) - 1;
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;

  TH1D *hist_acc_den_all = new TH1D("hist_acc_den_all", ";p_{T} (GeV/c);Counts", nPtBinsAll, ptBinsAll);
  TH1D *hist_acc_num_all = new TH1D("hist_acc_num_all", ";p_{T} (GeV/c);Counts", nPtBinsAll, ptBinsAll);
  TH1D *hist_acc_den_mid = new TH1D("hist_acc_den_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_acc_den_fwd = new TH1D("hist_acc_den_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_acc_num_mid = new TH1D("hist_acc_num_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_acc_num_fwd = new TH1D("hist_acc_num_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_acc_den_mid_Int = new TH1D("hist_acc_den_mid_Int", ";p_{T} integrated;Counts", 1, 0, 50);
  TH1D *hist_acc_den_fwd_Int = new TH1D("hist_acc_den_fwd_Int", ";p_{T} integrated;Counts", 1, 0, 50);
  TH1D *hist_acc_num_mid_Int = new TH1D("hist_acc_num_mid_Int", ";p_{T} integrated;Counts", 1, 0, 50);
  TH1D *hist_acc_num_fwd_Int = new TH1D("hist_acc_num_fwd_Int", ";p_{T} integrated;Counts", 1, 0, 50);

  const double massLow = 2.6;
  const double massHigh = 3.5;
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt < 0 || nEvt > nTot)
    nEvt = nTot;

  long denCount = 0;
  long numCount = 0;
  for (Long64_t iev = 0; iev < nEvt; ++iev)
  {
    if ((iev % 100000 == 0) && iev != 0)
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
           << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Acceptance denominator: " << denCount << "\n";
      cout << "Acceptance numerator: " << numCount << "\n";
    }

    oniaTree->GetEntry(iev);

    const Int_t Gen_QQ_size = GetGenQQSize();
    if (Gen_QQ_size < 0)
      continue;

    double eventWeightBase = 1.0;
    if (isGenW && hasGenWeightBranch)
      eventWeightBase *= Gen_weight;

    for (Int_t irqq = 0; irqq < Gen_QQ_size; ++irqq)
    {
      const Int_t iMuPl = hasGenQQIndexBranches ? GetGenQQMuplIdx(irqq) : -1;
      const Int_t iMuMi = hasGenQQIndexBranches ? GetGenQQMumiIdx(irqq) : -1;
      if (hasGenQQIndexBranches && (iMuPl < 0 || iMuMi < 0))
        continue;

      float mass = 0.f;
      float pt = 0.f;
      float eta = 0.f;
      float phi = 0.f;
      float pt1 = 0.f;
      float pt2 = 0.f;
      float eta1 = 0.f;
      float eta2 = 0.f;

      if (usesSplitMomentumBranches)
      {
        if (!hasGenQQIndexBranches)
          continue;
        if (!Gen_QQ_4mom_m || !Gen_QQ_4mom_pt || !Gen_QQ_4mom_eta || !Gen_QQ_4mom_phi ||
            !Gen_mu_4mom_pt || !Gen_mu_4mom_eta)
          continue;
        if (static_cast<size_t>(irqq) >= Gen_QQ_4mom_m->size() ||
            static_cast<size_t>(irqq) >= Gen_QQ_4mom_pt->size() ||
            static_cast<size_t>(irqq) >= Gen_QQ_4mom_eta->size() ||
            static_cast<size_t>(irqq) >= Gen_QQ_4mom_phi->size())
          continue;
        if (static_cast<size_t>(iMuPl) >= Gen_mu_4mom_pt->size() ||
            static_cast<size_t>(iMuMi) >= Gen_mu_4mom_pt->size() ||
            static_cast<size_t>(iMuPl) >= Gen_mu_4mom_eta->size() ||
            static_cast<size_t>(iMuMi) >= Gen_mu_4mom_eta->size())
          continue;
        mass = Gen_QQ_4mom_m->at(irqq);
        pt = Gen_QQ_4mom_pt->at(irqq);
        eta = Gen_QQ_4mom_eta->at(irqq);
        phi = Gen_QQ_4mom_phi->at(irqq);
        pt1 = Gen_mu_4mom_pt->at(iMuPl);
        pt2 = Gen_mu_4mom_pt->at(iMuMi);
        eta1 = Gen_mu_4mom_eta->at(iMuPl);
        eta2 = Gen_mu_4mom_eta->at(iMuMi);
      }
      else
      {
        if (!Gen_QQ_4mom)
          continue;
        if (irqq >= Gen_QQ_4mom->GetEntriesFast())
          continue;
        auto *qq = static_cast<TLorentzVector *>(Gen_QQ_4mom->At(irqq));
        TLorentzVector *muPl = nullptr;
        TLorentzVector *muMi = nullptr;
        if (hasGenQQIndexBranches)
        {
          if (!Gen_mu_4mom)
            continue;
          if (iMuPl >= Gen_mu_4mom->GetEntriesFast() ||
              iMuMi >= Gen_mu_4mom->GetEntriesFast())
            continue;
          muPl = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuPl));
          muMi = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuMi));
        }
        else
        {
          if (!Gen_QQ_mupl_4mom || !Gen_QQ_mumi_4mom)
            continue;
          if (irqq >= Gen_QQ_mupl_4mom->GetEntriesFast() ||
              irqq >= Gen_QQ_mumi_4mom->GetEntriesFast())
            continue;
          muPl = static_cast<TLorentzVector *>(Gen_QQ_mupl_4mom->At(irqq));
          muMi = static_cast<TLorentzVector *>(Gen_QQ_mumi_4mom->At(irqq));
        }
        if (!qq || !muPl || !muMi)
          continue;
        mass = qq->M();
        pt = qq->Pt();
        eta = qq->Eta();
        phi = qq->Phi();
        pt1 = muPl->Pt();
        pt2 = muMi->Pt();
        eta1 = muPl->Eta();
        eta2 = muMi->Eta();
      }

      if (mass < massLow || mass > massHigh)
        continue;
      if (hasGenQQIndexBranches && GetGenMuCharge(iMuPl) * GetGenMuCharge(iMuMi) > 0)
        continue;

      TLorentzVector qq;
      qq.SetPtEtaPhiM(pt, eta, phi, mass);
      const float absY = std::abs(qq.Rapidity());
      if (absY >= 2.4f)
        continue;

      const bool isMidRap = absY < 1.6f;
      const bool isAll = (pt >= 6.5f && pt < 50.f);
      bool isMid = false;
      bool isFwd = false;
      if (absY < 1.6f && pt >= 6.5f && pt < 40.f)
        isMid = true;
      if (absY >= 1.6f && absY < 2.4f && pt >= 3.5f && pt < 40.f)
        isFwd = true;
      if (!isAll && !isMid && !isFwd)
        continue;

      const double ptWeight = isPtW ? GetPtWeight(isMidRap ? fPtW_mid : fPtW_fwd, pt) : 1.0;
      const double weight = eventWeightBase * ptWeight;

      if (isAll)
        hist_acc_den_all->Fill(pt, weight);
      if (isMid)
      {
        hist_acc_den_mid->Fill(pt, weight);
        hist_acc_den_mid_Int->Fill(1.0, weight);
      }
      if (isFwd)
      {
        hist_acc_den_fwd->Fill(pt, weight);
        hist_acc_den_fwd_Int->Fill(1.0, weight);
      }
      ++denCount;

      if (!IsAcceptanceQQRun2018(pt1, eta1) || !IsAcceptanceQQRun2018(pt2, eta2))
        continue;

      if (isAll)
        hist_acc_num_all->Fill(pt, weight);
      if (isMid)
      {
        hist_acc_num_mid->Fill(pt, weight);
        hist_acc_num_mid_Int->Fill(1.0, weight);
      }
      if (isFwd)
      {
        hist_acc_num_fwd->Fill(pt, weight);
        hist_acc_num_fwd_Int->Fill(1.0, weight);
      }
      ++numCount;
    }
  }

  hist_acc_den_all->Scale(1.0, "width");
  hist_acc_den_mid->Scale(1.0, "width");
  hist_acc_den_fwd->Scale(1.0, "width");
  hist_acc_num_all->Scale(1.0, "width");
  hist_acc_num_mid->Scale(1.0, "width");
  hist_acc_num_fwd->Scale(1.0, "width");

  TH1D *hist_acc_all = BuildAcceptanceHist(hist_acc_num_all, hist_acc_den_all, "hist_acc_all", ";p_{T} (GeV/c);Acceptance");
  TH1D *hist_acc_mid = BuildAcceptanceHist(hist_acc_num_mid, hist_acc_den_mid, "hist_acc_mid", ";p_{T} (GeV/c);Acceptance");
  TH1D *hist_acc_fwd = BuildAcceptanceHist(hist_acc_num_fwd, hist_acc_den_fwd, "hist_acc_fwd", ";p_{T} (GeV/c);Acceptance");
  TH1D *hist_acc_mid_Int = BuildAcceptanceHist(hist_acc_num_mid_Int, hist_acc_den_mid_Int, "hist_acc_mid_Int", ";p_{T} integrated;Acceptance");
  TH1D *hist_acc_fwd_Int = BuildAcceptanceHist(hist_acc_num_fwd_Int, hist_acc_den_fwd_Int, "hist_acc_fwd_Int", ";p_{T} integrated;Acceptance");

  hist_acc_all->SetName("hAccPt_2021_ally");
  hist_acc_mid->SetName("hAccPt_2021_midy");
  hist_acc_fwd->SetName("hAccPt_2021_Fory");
  hist_acc_den_mid_Int->SetName("hDenPt_2021_midy_Int");
  hist_acc_den_fwd_Int->SetName("hDenPt_2021_Fory_Int");
  hist_acc_num_mid_Int->SetName("hNumPt_2021_midy_Int");
  hist_acc_num_fwd_Int->SetName("hNumPt_2021_Fory_Int");
  hist_acc_mid_Int->SetName("hAccPt_2021_midy_Int");
  hist_acc_fwd_Int->SetName("hAccPt_2021_Fory_Int");

  fout->cd();
  hist_acc_all->Write();
  hist_acc_mid->Write();
  hist_acc_fwd->Write();
  hist_acc_den_mid_Int->Write();
  hist_acc_den_fwd_Int->Write();
  hist_acc_num_mid_Int->Write();
  hist_acc_num_fwd_Int->Write();
  hist_acc_mid_Int->Write();
  hist_acc_fwd_Int->Write();
  fout->Close();

  cout << "Acceptance denominator count: " << denCount << "\n";
  cout << "Acceptance numerator count: " << numCount << "\n";
  cout << "Output ROOT: " << outFilePath << "\n";
  cout << "Finish acceptance_1d()\n";
}
