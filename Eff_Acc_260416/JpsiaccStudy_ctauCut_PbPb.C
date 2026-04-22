#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "../Style.h"
#include "../tdrstyle.C"

using namespace std;

namespace
{
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

  const TString typeName = leaf->GetTypeName();
  return typeName == "Short_t" || typeName == "UShort_t";
}

bool IsAcceptable(double pt, double eta)
{
  if ((fabs(eta) < 1.2 && pt >= 3.5) ||
      ((fabs(eta) >= 1.2 && fabs(eta) < 2.1) && pt >= (5.47 - 1.89 * fabs(eta))) ||
      ((fabs(eta) >= 2.1 && fabs(eta) < 2.4) && pt >= 1.5))
    return true;
  return false;
}

void PrintAcc(double bins[], TH1F *h)
{
  for (int i = 0; i < h->GetNbinsX(); i++)
  {
    cout << bins[i] << " - " << bins[i + 1] << " : "
         << h->GetBinContent(i + 1) << " +- " << h->GetBinError(i + 1) << endl;
  }
}

TFile *OpenFirstValid(const char *const *paths, int nPaths, const char *tag)
{
  for (int i = 0; i < nPaths; ++i)
  {
    TFile *f = TFile::Open(paths[i], "READ");
    if (f && !f->IsZombie())
    {
      cout << "[PtW] " << tag << " : " << paths[i] << endl;
      return f;
    }
    if (f)
    {
      f->Close();
      delete f;
    }
  }
  return nullptr;
}

TF1 *LoadClonedWeight(TFile *file, const char *sourceName, const char *cloneName)
{
  if (!file)
    return nullptr;

  TF1 *fIn = dynamic_cast<TF1 *>(file->Get(sourceName));
  if (!fIn)
    return nullptr;

  return static_cast<TF1 *>(fIn->Clone(cloneName));
}
} // namespace

void JpsiaccStudy_ctauCut_PbPb(int state = 1, int wtopt = 1, int isPtWgtUp = 0, TString rmk = "20260420")
{
  if (state != 1 && state != 2)
  {
    cout << "[ERROR] state must be 1(prompt) or 2(nonprompt). state=" << state << endl;
    return;
  }

  if (isPtWgtUp != 0)
  {
    cout << "[WARN] pt-weight systematic shift is disabled in this macro. isPtWgtUp is ignored." << endl;
    isPtWgtUp = 0;
  }

  TStopwatch timer;
  timer.Start();

  gStyle->SetOptStat(0);
  setTDRStyle();

  gSystem->mkdir("./figs", true);
  gSystem->mkdir("./roots", true);

  const char *inputPrompt = "/data/Oniatree/Jpsi/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  const char *inputNonPrompt = "/data/Oniatree/Jpsi/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  const char *inputPath = (state == 1 ? inputPrompt : inputNonPrompt);

  TFile *rf = TFile::Open(inputPath, "READ");
  if (!rf || rf->IsZombie())
  {
    cout << "[ERROR] cannot open input file: " << inputPath << endl;
    return;
  }
  TTree *tree = (TTree *)rf->Get("hionia/myTree");
  if (!tree)
  {
    cout << "[ERROR] cannot find tree hionia/myTree in " << inputPath << endl;
    return;
  }

  TF1 *fptwMid = nullptr;
  TF1 *fptwFor = nullptr;
  TFile *fPtWMid = nullptr;
  TFile *fPtWFor = nullptr;

  if (wtopt == 1)
  {
    const char *pbPRMid[] = {
        "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260410.root",
    };
    const char *pbPRFor[] = {
        "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_260410.root",
    };
    const char *pbNPMid[] = {
        "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260410.root",
    };
    const char *pbNPFor[] = {
        "/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_260410.root",
    };

    if (state == 1)
    {
      fPtWMid = OpenFirstValid(pbPRMid, 1, "PbPb prompt mid");
      fPtWFor = OpenFirstValid(pbPRFor, 1, "PbPb prompt forward");
    }
    else
    {
      fPtWMid = OpenFirstValid(pbNPMid, 1, "PbPb nonprompt mid");
      fPtWFor = OpenFirstValid(pbNPFor, 1, "PbPb nonprompt forward");
    }

    if (!fPtWMid || !fPtWFor)
    {
      cout << "[ERROR] failed to open rapidity-dependent pT weight files for PbPb." << endl;
      return;
    }

    fptwMid = LoadClonedWeight(fPtWMid, "dataMC_Ratio1", "fptwMid_pbpb");
    fptwFor = LoadClonedWeight(fPtWFor, "dataMC_Ratio1", "fptwFor_pbpb");
    if (!fptwMid || !fptwFor)
    {
      cout << "[ERROR] dataMC_Ratio1 is missing in one of pt-weight files." << endl;
      return;
    }
    fPtWMid->Close();
    fPtWFor->Close();
    delete fPtWMid;
    delete fPtWFor;
  }

  Int_t Gen_QQ_size;
  TClonesArray *Gen_QQ_4mom = nullptr;
  TClonesArray *Gen_QQ_mupl_4mom = nullptr;
  TClonesArray *Gen_QQ_mumi_4mom = nullptr;
  const bool hasGenQQIndexBranches = (tree->GetBranch("Gen_QQ_mupl_idx") != nullptr &&
                                      tree->GetBranch("Gen_QQ_mumi_idx") != nullptr &&
                                      tree->GetBranch("Gen_mu_charge") != nullptr);
  const bool useShortIndexBranches = hasGenQQIndexBranches && BranchUsesShort(tree, "Gen_QQ_mupl_idx");
  const int maxBranchSize = 1000;
  Short_t Gen_QQ_mupl_idx_short[maxBranchSize];
  Short_t Gen_QQ_mumi_idx_short[maxBranchSize];
  Short_t Gen_mu_charge_short[maxBranchSize];
  Int_t Gen_QQ_mupl_idx_int[maxBranchSize];
  Int_t Gen_QQ_mumi_idx_int[maxBranchSize];
  Int_t Gen_mu_charge_int[maxBranchSize];

  tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
  tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
  tree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);
  tree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
  if (hasGenQQIndexBranches && useShortIndexBranches)
  {
    tree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_short);
    tree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_short);
    tree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_short);
  }
  else if (hasGenQQIndexBranches)
  {
    tree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_int);
    tree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_int);
    tree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_int);
  }

  TLorentzVector *JP = nullptr;
  TLorentzVector *Mu1 = nullptr;
  TLorentzVector *Mu2 = nullptr;
  const double massLow = 2.6;
  const double massHigh = 3.5;

  double aDenPt_2021_allybin[] = {6.5, 7.5, 8.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 50.0};
  double aDenPt_2021_midybin[] = {6.5, 9, 12, 15, 20, 25, 40};
  double aDenPt_2021_Forybin[] = {3.5, 6.5, 9, 12, 40};

  TH1F *hDenPt_2021_ally = new TH1F("hDenPt_2021_ally", ";p_{T} (GeV/c};", 10, aDenPt_2021_allybin);
  TH1F *hNumPt_2021_ally = new TH1F("hNumPt_2021_ally", ";p_{T} (GeV/c};", 10, aDenPt_2021_allybin);
  TH1F *hAccPt_2021_ally = new TH1F("hAccPt_2021_ally", ";p_{T} (GeV/c};", 10, aDenPt_2021_allybin);
  TH1F *hDenPt_2021_midy = new TH1F("hDenPt_2021_midy", ";p_{T} (GeV/c};", 6, aDenPt_2021_midybin);
  TH1F *hNumPt_2021_midy = new TH1F("hNumPt_2021_midy", ";p_{T} (GeV/c};", 6, aDenPt_2021_midybin);
  TH1F *hAccPt_2021_midy = new TH1F("hAccPt_2021_midy", ";p_{T} (GeV/c};", 6, aDenPt_2021_midybin);
  TH1F *hDenPt_2021_Fory = new TH1F("hDenPt_2021_Fory", ";p_{T} (GeV/c};", 4, aDenPt_2021_Forybin);
  TH1F *hNumPt_2021_Fory = new TH1F("hNumPt_2021_Fory", ";p_{T} (GeV/c};", 4, aDenPt_2021_Forybin);
  TH1F *hAccPt_2021_Fory = new TH1F("hAccPt_2021_Fory", ";p_{T} (GeV/c};", 4, aDenPt_2021_Forybin);

  TH1F *hDenPt_2021_Fory_Int = new TH1F("hDenPt_2021_Fory_Int", ";p_{T} (GeV/c};", 1, 0, 50);
  TH1F *hNumPt_2021_Fory_Int = new TH1F("hNumPt_2021_Fory_Int", ";p_{T} (GeV/c};", 1, 0, 50);
  TH1F *hAccPt_2021_Fory_Int = new TH1F("hAccPt_2021_Fory_Int", ";p_{T} (GeV/c};", 1, 0, 50);
  TH1F *hDenPt_2021_midy_Int = new TH1F("hDenPt_2021_midy_Int", ";p_{T} (GeV/c};", 1, 0, 50);
  TH1F *hNumPt_2021_midy_Int = new TH1F("hNumPt_2021_midy_Int", ";p_{T} (GeV/c};", 1, 0, 50);
  TH1F *hAccPt_2021_midy_Int = new TH1F("hAccPt_2021_midy_Int", ";p_{T} (GeV/c};", 1, 0, 50);

  hDenPt_2021_ally->Sumw2();
  hNumPt_2021_ally->Sumw2();
  hDenPt_2021_midy->Sumw2();
  hNumPt_2021_midy->Sumw2();
  hDenPt_2021_Fory->Sumw2();
  hNumPt_2021_Fory->Sumw2();
  hDenPt_2021_Fory_Int->Sumw2();
  hNumPt_2021_Fory_Int->Sumw2();
  hDenPt_2021_midy_Int->Sumw2();
  hNumPt_2021_midy_Int->Sumw2();

  const Long64_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  for (Long64_t i = 0; i < nEvt; ++i)
  {
    tree->GetEntry(i);
    if (i % 100000 == 0)
      cout << ">>>>> EVENT " << i << " / " << nEvt << endl;

    for (int j = 0; j < Gen_QQ_size; ++j)
    {
      JP = (TLorentzVector *)Gen_QQ_4mom->At(j);
      Mu1 = (TLorentzVector *)Gen_QQ_mupl_4mom->At(j);
      Mu2 = (TLorentzVector *)Gen_QQ_mumi_4mom->At(j);
      if (!JP || !Mu1 || !Mu2)
        continue;

      const double pt = JP->Pt();
      const double mass = JP->M();
      if (mass < massLow || mass > massHigh)
        continue;

      if (hasGenQQIndexBranches)
      {
        const int muplIdx = useShortIndexBranches ? static_cast<int>(Gen_QQ_mupl_idx_short[j]) : Gen_QQ_mupl_idx_int[j];
        const int mumiIdx = useShortIndexBranches ? static_cast<int>(Gen_QQ_mumi_idx_short[j]) : Gen_QQ_mumi_idx_int[j];
        if (muplIdx < 0 || mumiIdx < 0 || muplIdx >= maxBranchSize || mumiIdx >= maxBranchSize)
          continue;

        const int muplCharge = useShortIndexBranches ? static_cast<int>(Gen_mu_charge_short[muplIdx]) : Gen_mu_charge_int[muplIdx];
        const int mumiCharge = useShortIndexBranches ? static_cast<int>(Gen_mu_charge_short[mumiIdx]) : Gen_mu_charge_int[mumiIdx];
        if (muplCharge * mumiCharge > 0)
          continue;
      }

      const double absY = fabs(JP->Rapidity());
      if (absY >= 2.4)
        continue;

      const bool isMid = (absY < 1.6 && pt >= 6.5 && pt < 40.0);
      const bool isFor = (absY >= 1.6 && absY < 2.4 && pt >= 3.5 && pt < 40.0);
      if (!isMid && !isFor)
        continue;

      double wtMid = 1.0;
      double wtFor = 1.0;
      if (wtopt == 1)
      {
        wtMid = (fptwMid ? fptwMid->Eval(pt) : 1.0);
        wtFor = (fptwFor ? fptwFor->Eval(pt) : 1.0);
      }
      const double wtEvt = (isMid ? wtMid : wtFor);

      hDenPt_2021_ally->Fill(pt, wtEvt);
      if (isFor)
      {
        hDenPt_2021_Fory->Fill(pt, wtFor);
        if (pt > 3.5 && pt < 40)
          hDenPt_2021_Fory_Int->Fill(1, wtFor);
      }
      if (isMid)
      {
        hDenPt_2021_midy->Fill(pt, wtMid);
        if (pt > 6.5 && pt < 40)
          hDenPt_2021_midy_Int->Fill(1, wtMid);
      }

      const bool mu1pass = IsAcceptable(Mu1->Pt(), Mu1->Eta());
      const bool mu2pass = IsAcceptable(Mu2->Pt(), Mu2->Eta());
      if (!mu1pass || !mu2pass)
        continue;

      hNumPt_2021_ally->Fill(pt, wtEvt);
      if (isFor)
      {
        hNumPt_2021_Fory->Fill(pt, wtFor);
        if (pt > 3.5 && pt < 40)
          hNumPt_2021_Fory_Int->Fill(1, wtFor);
      }
      if (isMid)
      {
        hNumPt_2021_midy->Fill(pt, wtMid);
        if (pt > 6.5 && pt < 40)
          hNumPt_2021_midy_Int->Fill(1, wtMid);
      }
    }
  }

  hAccPt_2021_ally->Divide(hNumPt_2021_ally, hDenPt_2021_ally, 1, 1, "B");
  hAccPt_2021_midy->Divide(hNumPt_2021_midy, hDenPt_2021_midy, 1, 1, "B");
  hAccPt_2021_Fory->Divide(hNumPt_2021_Fory, hDenPt_2021_Fory, 1, 1, "B");
  hAccPt_2021_midy_Int->Divide(hNumPt_2021_midy_Int, hDenPt_2021_midy_Int, 1, 1, "B");
  hAccPt_2021_Fory_Int->Divide(hNumPt_2021_Fory_Int, hDenPt_2021_Fory_Int, 1, 1, "B");

  cout << "--- p_{T} binning, |y| < 2.4" << endl;
  PrintAcc(aDenPt_2021_allybin, hAccPt_2021_ally);
  cout << "--- p_{T} binning, |y| < 1.6" << endl;
  PrintAcc(aDenPt_2021_midybin, hAccPt_2021_midy);
  cout << "--- p_{T} binning, 1.6 < |y| < 2.4" << endl;
  PrintAcc(aDenPt_2021_Forybin, hAccPt_2021_Fory);

  const TString flavor = (state == 1 ? "PromptJpsi" : "BtoJpsi");
  TFile *wf = new TFile(Form("./roots/acceptance_%s_GenOnly_wgt%d_PbPb_SysUp%d_%s.root",
                             flavor.Data(), wtopt, isPtWgtUp, rmk.Data()),
                        "RECREATE");
  hAccPt_2021_ally->Write();
  hAccPt_2021_midy->Write();
  hAccPt_2021_Fory->Write();
  hDenPt_2021_midy_Int->Write();
  hDenPt_2021_Fory_Int->Write();
  hNumPt_2021_midy_Int->Write();
  hNumPt_2021_Fory_Int->Write();
  hAccPt_2021_midy_Int->Write();
  hAccPt_2021_Fory_Int->Write();
  wf->Write();
  wf->Close();

  TCanvas *cAll = new TCanvas("cAll_pbpb", "", 700, 600);
  hAccPt_2021_ally->GetXaxis()->SetRangeUser(6.5, 50.0);
  hAccPt_2021_ally->GetYaxis()->SetRangeUser(0.0, 1.3);
  hAccPt_2021_ally->Draw("E");
  cAll->SaveAs(Form("./figs/AccPt_AllY_PbPb_state%d_wgt%d_%s.png", state, wtopt, rmk.Data()));

  TCanvas *cMid = new TCanvas("cMid_pbpb", "", 700, 600);
  hAccPt_2021_midy->GetXaxis()->SetRangeUser(6.5, 50.0);
  hAccPt_2021_midy->GetYaxis()->SetRangeUser(0.0, 1.3);
  hAccPt_2021_midy->Draw("E");
  cMid->SaveAs(Form("./figs/AccPt_MidY_PbPb_state%d_wgt%d_%s.png", state, wtopt, rmk.Data()));

  TCanvas *cFor = new TCanvas("cFor_pbpb", "", 700, 600);
  hAccPt_2021_Fory->GetXaxis()->SetRangeUser(3.0, 50.0);
  hAccPt_2021_Fory->GetYaxis()->SetRangeUser(0.0, 1.3);
  hAccPt_2021_Fory->Draw("E");
  cFor->SaveAs(Form("./figs/AccPt_ForY_PbPb_state%d_wgt%d_%s.png", state, wtopt, rmk.Data()));

  timer.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", timer.RealTime(), timer.CpuTime());
  cout << "output File: " << Form("./roots/acceptance_%s_GenOnly_wgt%d_PbPb_SysUp%d_%s.root",
                             flavor.Data(), wtopt, isPtWgtUp, rmk.Data()) << endl;
}
