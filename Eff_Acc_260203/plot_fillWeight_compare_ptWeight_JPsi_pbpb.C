#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../cutsAndBin.h"

using namespace std;

namespace {
  TH1D *MakeRatio(TH1D *num, TH1D *den, const char *name)
  {
    TH1D *h = (TH1D *)num->Clone(name);
    h->Reset();
    for (int ib = 1; ib <= num->GetNbinsX(); ++ib)
    {
      const double n = num->GetBinContent(ib);
      const double d = den->GetBinContent(ib);
      if (d > 0)
      {
        h->SetBinContent(ib, n / d);
      }
    }
    return h;
  }
}

void plot_fillWeight_compare_ptWeight_JPsi_pbpb(
    int state = 1,
    float cLow = 0, float cHigh = 180,
    int maxEvents = 5000000)
{
  gStyle->SetOptStat(0);

  TString sampleTag = (state == 1) ? "PR" : "NP";
  TString inputMC = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
  TString ptWeightFile = "../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
  if (state == 2)
  {
    inputMC = "/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root";
    ptWeightFile = "../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
  }

  TFile *fPtW = TFile::Open(ptWeightFile, "READ");
  if (!fPtW || fPtW->IsZombie())
  {
    cout << "[ERROR] Cannot open pt-weight file: " << ptWeightFile.Data() << endl;
    return;
  }
  TF1 *fptw = (TF1 *)fPtW->Get("dataMC_Ratio1");
  if (!fptw)
  {
    cout << "[ERROR] Cannot find TF1 dataMC_Ratio1 in: " << ptWeightFile.Data() << endl;
    return;
  }

  TChain *mytree = new TChain("hionia/myTree");
  mytree->Add(inputMC.Data());
  const int totalEntries = mytree->GetEntries();
  const int nevt = (maxEvents > 0 && maxEvents < totalEntries) ? maxEvents : totalEntries;
  cout << "[INFO] Entries in tree: " << totalEntries << ", processing: " << nevt << endl;

  const int maxBranchSize = 1000;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Float_t Gen_weight;

  Short_t Gen_QQ_size;
  Short_t Gen_mu_size;
  TClonesArray *Gen_QQ_4mom = 0;
  TClonesArray *Gen_mu_4mom = 0;
  Short_t Gen_QQ_mupl_idx[maxBranchSize];
  Short_t Gen_QQ_mumi_idx[maxBranchSize];
  Short_t Gen_mu_charge[maxBranchSize];

  Short_t Reco_QQ_size;
  TClonesArray *Reco_QQ_4mom = 0;
  TClonesArray *Reco_mu_4mom = 0;
  ULong64_t Reco_QQ_trig[maxBranchSize];
  Float_t Reco_QQ_VtxProb[maxBranchSize];
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_whichGen[maxBranchSize];
  Short_t Reco_QQ_sign[maxBranchSize];
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Short_t Reco_mu_whichGen[maxBranchSize];
  Float_t Reco_mu_dxy[maxBranchSize];
  Float_t Reco_mu_dz[maxBranchSize];
  Int_t Reco_mu_nTrkWMea[maxBranchSize];
  Int_t Reco_mu_nPixWMea[maxBranchSize];

  mytree->SetBranchAddress("Centrality", &Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers);
  mytree->SetBranchAddress("Gen_weight", &Gen_weight);

  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
  mytree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);
  mytree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge);

  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  mytree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  mytree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);

  TH1D *hGenNoPt = new TH1D("hGenNoPt", "GEN fill-weight;log_{10}(fill weight);Candidates", 120, -4, 5);
  TH1D *hGenWithPt = new TH1D("hGenWithPt", "GEN fill-weight;log_{10}(fill weight);Candidates", 120, -4, 5);
  TH1D *hRecoNoPt = new TH1D("hRecoNoPt", "RECO fill-weight;log_{10}(fill weight);Candidates", 120, -4, 5);
  TH1D *hRecoWithPt = new TH1D("hRecoWithPt", "RECO fill-weight;log_{10}(fill weight);Candidates", 120, -4, 5);

  const int kTrigSel = 12;
  const double massLow = 2.9;
  const double massHigh = 3.3;

  for (int iev = 0; iev < nevt; ++iev)
  {
    if (iev % 1000000 == 0)
      cout << "[INFO] EVENT " << iev << " / " << nevt << endl;

    mytree->GetEntry(iev);
    if (!(Centrality > cLow && Centrality < cHigh))
      continue;

    const double baseWeight = findNcoll(Centrality) * Gen_weight;
    bool HLTPass = ((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

    for (int igen = 0; igen < Gen_QQ_size; ++igen)
    {
      TLorentzVector *JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
      TLorentzVector *mupl_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      TLorentzVector *mumi_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      const double rap = fabs(JP_Gen->Rapidity());
      const double pt = JP_Gen->Pt();
      if (!(JP_Gen->M() > massLow && JP_Gen->M() < massHigh))
        continue;
      if (!(pt > 3.0 && pt < 50.0 && rap < 2.4))
        continue;
      if (!(IsAcceptanceQQ(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptanceQQ(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
        continue;
      if (!(fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
        continue;
      if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
        continue;

      const bool passHptGen0 = ((pt > 3.0 && pt < 6.5 && rap > 1.6 && rap < 2.4) || (pt > 6.5 && pt < 50.0 && rap < 2.4));
      if (!passHptGen0)
        continue;

      const double ptWeight = fptw->Eval(pt);
      if (baseWeight > 0)
      {
        hGenNoPt->Fill(log10(baseWeight));
        hGenWithPt->Fill(log10(baseWeight * ptWeight));
      }
    }

    if (!HLTPass)
      continue;

    for (int irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      TLorentzVector *JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
      TLorentzVector *mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      TLorentzVector *mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      bool HLTFilterPass = ((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
      if (!HLTFilterPass)
        continue;

      if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
        continue;
      if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
        continue;
      if (Reco_QQ_whichGen[irqq] == -1)
        continue;
      if (Reco_QQ_sign[irqq] != 0)
        continue;
      if (Reco_QQ_VtxProb[irqq] < 0.01)
        continue;

      const double rap = fabs(JP_Reco->Rapidity());
      const double pt = JP_Reco->Pt();
      if (rap > 2.4)
        continue;
      if (pt < 3.0 || pt > 50.0)
        continue;
      if (JP_Reco->M() < massLow || JP_Reco->M() > massHigh)
        continue;

      if (!(IsAcceptanceQQ(mupl_Reco->Pt(), mupl_Reco->Eta()) && IsAcceptanceQQ(mumi_Reco->Pt(), mumi_Reco->Eta())))
        continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));
      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));
      if (!(passMuonTypePl && passMuonTypeMi))
        continue;

      const bool muplSoft =
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.);
      const bool mumiSoft =
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.);
      if (!(muplSoft && mumiSoft))
        continue;

      const bool passHptReco0 = ((rap > 1.6 && rap < 2.4 && pt > 3.0 && pt < 6.5) || (rap < 2.4 && pt > 6.5 && pt < 50.0));
      if (!passHptReco0)
        continue;

      const double ptWeight = fptw->Eval(pt);
      if (baseWeight > 0)
      {
        hRecoNoPt->Fill(log10(baseWeight));
        hRecoWithPt->Fill(log10(baseWeight * ptWeight));
      }
    }
  }

  TH1D *hGenNoPtN = (TH1D *)hGenNoPt->Clone("hGenNoPtN");
  TH1D *hGenWithPtN = (TH1D *)hGenWithPt->Clone("hGenWithPtN");
  TH1D *hRecoNoPtN = (TH1D *)hRecoNoPt->Clone("hRecoNoPtN");
  TH1D *hRecoWithPtN = (TH1D *)hRecoWithPt->Clone("hRecoWithPtN");
  if (hGenNoPtN->Integral() > 0)
    hGenNoPtN->Scale(1. / hGenNoPtN->Integral());
  if (hGenWithPtN->Integral() > 0)
    hGenWithPtN->Scale(1. / hGenWithPtN->Integral());
  if (hRecoNoPtN->Integral() > 0)
    hRecoNoPtN->Scale(1. / hRecoNoPtN->Integral());
  if (hRecoWithPtN->Integral() > 0)
    hRecoWithPtN->Scale(1. / hRecoWithPtN->Integral());

  hGenNoPtN->SetLineColor(kBlue + 1);
  hGenNoPtN->SetLineWidth(2);
  hGenWithPtN->SetLineColor(kRed + 1);
  hGenWithPtN->SetLineWidth(2);
  hRecoNoPtN->SetLineColor(kBlue + 1);
  hRecoNoPtN->SetLineWidth(2);
  hRecoWithPtN->SetLineColor(kRed + 1);
  hRecoWithPtN->SetLineWidth(2);

  TH1D *hRatioGen = MakeRatio(hGenWithPtN, hGenNoPtN, "hRatioGen");
  TH1D *hRatioReco = MakeRatio(hRecoWithPtN, hRecoNoPtN, "hRatioReco");
  hRatioGen->SetTitle("GEN normalized ratio (with/without)");
  hRatioReco->SetTitle("RECO normalized ratio (with/without)");
  hRatioGen->GetYaxis()->SetRangeUser(0.6, 1.4);
  hRatioReco->GetYaxis()->SetRangeUser(0.6, 1.4);
  hRatioGen->GetXaxis()->SetTitle("log_{10}(fill weight)");
  hRatioReco->GetXaxis()->SetTitle("log_{10}(fill weight)");
  hRatioGen->GetYaxis()->SetTitle("Ratio");
  hRatioReco->GetYaxis()->SetTitle("Ratio");

  gSystem->mkdir("./figs", kTRUE);
  gSystem->mkdir("./roots", kTRUE);

  TCanvas *c = new TCanvas("c_fillWeight_ptCompare", "c_fillWeight_ptCompare", 1400, 900);
  c->Divide(2, 2);

  TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(hGenNoPtN, "without pt weight", "l");
  leg->AddEntry(hGenWithPtN, "with pt weight", "l");

  c->cd(1);
  hGenNoPtN->SetTitle("GEN fill-weight comparison");
  hGenNoPtN->Draw("hist");
  hGenWithPtN->Draw("hist same");
  leg->Draw();

  c->cd(2);
  hRatioGen->Draw("hist");

  c->cd(3);
  hRecoNoPtN->SetTitle("RECO fill-weight comparison");
  hRecoNoPtN->Draw("hist");
  hRecoWithPtN->Draw("hist same");
  leg->Draw();

  c->cd(4);
  hRatioReco->Draw("hist");

  TLatex lt;
  lt.SetNDC();
  lt.SetTextSize(0.03);
  c->cd(1);
  lt.DrawLatex(0.14, 0.86, Form("%s J/#psi, cent %.0f-%.0f%%", (state == 1 ? "Prompt" : "Non-prompt"), cLow / 2.0, cHigh / 2.0));

  TString outPng = Form("./figs/fillWeight_compare_ptWeight_JPsi_pbpb_%s_cent%.0fto%.0f.png", sampleTag.Data(), cLow, cHigh);
  c->SaveAs(outPng);

  TString outRoot = Form("./roots/fillWeight_compare_ptWeight_JPsi_pbpb_%s_cent%.0fto%.0f.root", sampleTag.Data(), cLow, cHigh);
  TFile *fout = new TFile(outRoot, "RECREATE");
  hGenNoPt->Write();
  hGenWithPt->Write();
  hRecoNoPt->Write();
  hRecoWithPt->Write();
  hGenNoPtN->Write();
  hGenWithPtN->Write();
  hRecoNoPtN->Write();
  hRecoWithPtN->Write();
  hRatioGen->Write();
  hRatioReco->Write();
  c->Write();
  fout->Write();
  fout->Close();

  cout << "[INFO] Saved figure : " << outPng.Data() << endl;
  cout << "[INFO] Saved ROOT   : " << outRoot.Data() << endl;
  cout << "[INFO] GEN mean log10(weight) noPt/withPt = "
       << hGenNoPt->GetMean() << " / " << hGenWithPt->GetMean() << endl;
  cout << "[INFO] RECO mean log10(weight) noPt/withPt = "
       << hRecoNoPt->GetMean() << " / " << hRecoWithPt->GetMean() << endl;
}
