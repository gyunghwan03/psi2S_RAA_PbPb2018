#include <iostream>
#include <vector>
#include <cmath>

#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../cutsAndBin.h"

using namespace std;

namespace {

struct ChannelConfig {
  TString tag;
  bool isPbPb;
  bool useRapSplitPtW;
  int kTrigSel;
  double massLow;
  double massHigh;
  double ptPlotMin;
  double ptPlotMax;

  TString inputPrompt;
  TString inputNonPrompt;

  TString ptwPrompt;
  TString ptwNonPrompt;
  TString ptwPromptMid;
  TString ptwPromptFor;
  TString ptwNonPromptMid;
  TString ptwNonPromptFor;
};

TH1D *MakeRatio(TH1D *num, TH1D *den, const char *name)
{
  TH1D *h = (TH1D *)num->Clone(name);
  h->Reset();
  for (int ib = 1; ib <= num->GetNbinsX(); ++ib)
  {
    const double n = num->GetBinContent(ib);
    const double d = den->GetBinContent(ib);
    if (d != 0)
    {
      h->SetBinContent(ib, n / d);
    }
  }
  return h;
}

TH1D *ProfileToHist(TProfile *p, const char *name)
{
  TH1D *h = (TH1D *)p->ProjectionX(name);
  h->Reset();
  for (int ib = 1; ib <= p->GetNbinsX(); ++ib)
  {
    h->SetBinContent(ib, p->GetBinContent(ib));
    h->SetBinError(ib, p->GetBinError(ib));
  }
  return h;
}

void PrintProfileSummary(const char *label, TProfile *p)
{
  cout << "[SUMMARY] " << label << endl;
  for (int ib = 1; ib <= p->GetNbinsX(); ++ib)
  {
    const double xlow = p->GetXaxis()->GetBinLowEdge(ib);
    const double xhigh = p->GetXaxis()->GetBinUpEdge(ib);
    const double mean = p->GetBinContent(ib);
    const double err = p->GetBinError(ib);
    const double n = p->GetBinEntries(ib);
    cout << "  bin " << ib << " [" << xlow << ", " << xhigh << ")"
         << " : <ptWeight>=" << mean
         << " +/- " << err
         << " (entries=" << n << ")" << endl;
  }
}

void PrintTripleProfileSummary(const char *label, TProfile *pPtW, TProfile *pBase, TProfile *pWeight)
{
  cout << "[SUMMARY] " << label << endl;
  const int nBins = pPtW->GetNbinsX();
  for (int ib = 1; ib <= nBins; ++ib)
  {
    const double xlow = pPtW->GetXaxis()->GetBinLowEdge(ib);
    const double xhigh = pPtW->GetXaxis()->GetBinUpEdge(ib);
    const double n = pPtW->GetBinEntries(ib);
    cout << "  bin " << ib << " [" << xlow << ", " << xhigh << ")"
         << " : pT Weight=" << pPtW->GetBinContent(ib)
         << ", Gen Weight=" << pBase->GetBinContent(ib)
         << ", Weight=" << pWeight->GetBinContent(ib)
         << " (entries=" << n << ")" << endl;
  }
}

void StyleProfileForDraw(TProfile *p, int color, int marker)
{
  p->SetLineColor(color);
  p->SetMarkerColor(color);
  p->SetMarkerStyle(marker);
  p->SetMarkerSize(1.1);
  p->SetLineWidth(2);
}

void SaveTripleProfilePlot(TProfile *pPtWeight,
                           TProfile *pGenWeight,
                           TProfile *pWeight,
                           const TString &title,
                           const TString &outPath,
                           const TString &xTitle)
{
  TCanvas *c = new TCanvas(Form("c_%s", pPtWeight->GetName()), "", 850, 700);

  StyleProfileForDraw(pPtWeight, kBlue + 1, 20);
  StyleProfileForDraw(pGenWeight, kRed + 1, 21);
  StyleProfileForDraw(pWeight, kBlack, 33);

  pPtWeight->SetTitle(title);
  pPtWeight->GetXaxis()->SetTitle(xTitle);
  pPtWeight->GetYaxis()->SetTitle("Mean value");
  pPtWeight->Draw("E1");
  pGenWeight->Draw("E1 SAME");
  pWeight->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(pPtWeight, "pT Weight", "lep");
  leg->AddEntry(pGenWeight, "Gen Weight", "lep");
  leg->AddEntry(pWeight, "Weight", "lep");
  leg->Draw();

  c->SaveAs(outPath);
  delete leg;
  delete c;
}

bool IsJpsiPbpbRegion(double rap, double pt)
{
  return ((pt > 3.0 && pt < 6.5 && rap > 1.6 && rap < 2.4) ||
          (pt > 6.5 && pt < 50.0 && rap < 2.4));
}

bool IsPsi2SPtRegion(double rap, double pt)
{
  return ((pt > 3.5 && pt < 40.0 && rap > 1.6 && rap < 2.4) ||
          (pt > 6.5 && pt < 40.0 && rap < 1.6));
}

bool PassChannelPtRegion(const TString &tag, double rap, double pt)
{
  if (tag == "JPsi_pbpb") return IsJpsiPbpbRegion(rap, pt);
  if (tag == "Jpsi_pp") return (pt > 0.0 && pt < 50.0 && rap < 2.4);
  if (tag == "psi2S_pp") return IsPsi2SPtRegion(rap, pt);
  if (tag == "psi2S_pbpb") return IsPsi2SPtRegion(rap, pt);
  return false;
}

double EvalPtWeight(const ChannelConfig &cfg, TF1 *fSingle, TF1 *fMid, TF1 *fFor, double rap, double pt)
{
  if (!cfg.useRapSplitPtW)
  {
    if (!fSingle) return 1.0;
    return fSingle->Eval(pt);
  }

  if (rap >= 1.6 && rap < 2.4)
  {
    if (!fFor) return 1.0;
    return fFor->Eval(pt);
  }

  if (rap < 2.4)
  {
    if (!fMid) return 1.0;
    return fMid->Eval(pt);
  }

  return 1.0;
}

bool BuildConfig(const TString &channel, ChannelConfig &cfg)
{
  if (channel == "JPsi_pbpb")
  {
    cfg.tag = "JPsi_pbpb";
    cfg.isPbPb = true;
    cfg.useRapSplitPtW = false;
    cfg.kTrigSel = 12;
    cfg.massLow = 2.9;
    cfg.massHigh = 3.3;
    cfg.ptPlotMin = 3.0;
    cfg.ptPlotMax = 50.0;

    cfg.inputPrompt = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
    cfg.inputNonPrompt = "/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root";

    cfg.ptwPrompt = "../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
    cfg.ptwNonPrompt = "../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
    return true;
  }

  if (channel == "Jpsi_pp")
  {
    cfg.tag = "Jpsi_pp";
    cfg.isPbPb = false;
    cfg.useRapSplitPtW = false;
    cfg.kTrigSel = 3;
    cfg.massLow = 0.0;
    cfg.massHigh = 10.0;
    cfg.ptPlotMin = 0.0;
    cfg.ptPlotMax = 50.0;

    cfg.inputPrompt = "/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root";
    cfg.inputNonPrompt = "/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root";

    cfg.ptwPrompt = "../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
    cfg.ptwNonPrompt = "../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root";
    return true;
  }

  if (channel == "psi2S_pp")
  {
    cfg.tag = "psi2S_pp";
    cfg.isPbPb = false;
    cfg.useRapSplitPtW = true;
    cfg.kTrigSel = 3;
    cfg.massLow = 3.5;
    cfg.massHigh = 4.0;
    cfg.ptPlotMin = 3.0;
    cfg.ptPlotMax = 40.0;

    cfg.inputPrompt = "/data/Oniatree/psi2S/OniatreeMC_Psi2SMM_TuneCUETP8M1_5p02TeV_pythia8_RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2.root";
    cfg.inputNonPrompt = "/data/Oniatree/miniAOD/OniatreeMC_BPsi2SMM_TuneCUETP8M1_5p02TeV_pythia8_ptHatMin2_ONIATREE.root";

    cfg.ptwPromptMid = "../compareDataToMC/ratioDataMC_pp_psi2S_DATA_ctauCut_y0_2p4_260323_2exp.root";
    cfg.ptwPromptFor = "../compareDataToMC/ratioDataMC_pp_psi2S_DATA_ctauCut_y1p6_2p4_260323_2exp.root";

    cfg.ptwNonPromptMid = "../compareDataToMC/ratioDataMC_pp_Btopsi2S_DATA_ctauCut_y0_2p4_260323_2exp.root";
    cfg.ptwNonPromptFor = "../compareDataToMC/ratioDataMC_pp_Btopsi2S_DATA_ctauCut_y1p6_2p4_260323_2exp.root";
    return true;
  }

  if (channel == "psi2S_pbpb")
  {
    cfg.tag = "psi2S_pbpb";
    cfg.isPbPb = true;
    cfg.useRapSplitPtW = true;
    cfg.kTrigSel = 12;
    cfg.massLow = 3.5;
    cfg.massHigh = 4.0;
    cfg.ptPlotMin = 3.0;
    cfg.ptPlotMax = 40.0;

    cfg.inputPrompt = "/data/Oniatree/miniAOD/OniaTree_Run2_PbPb2018MC_Prompt_forPsi2S_miniAOD112X_muonSelGlb_13May23v2.root";
    cfg.inputNonPrompt = "/data/Oniatree/miniAOD/OniaTree_Run2_PbPb2018MC_BToPsi_forPsi2S_miniAOD112X_muonSelGlb_13May23v2.root";

    cfg.ptwPromptMid = "../compareDataToMC/ratioDataMC_AA_psi2S_DATA_ctauCut_y0_2p4_260323_2exp.root";
    cfg.ptwPromptFor = "../compareDataToMC/ratioDataMC_AA_psi2S_DATA_ctauCut_y1p6_2p4_260323_2exp.root";

    cfg.ptwNonPromptMid = "../compareDataToMC/ratioDataMC_AA_Btopsi2S_DATA_ctauCut_y0_2p4_260323_2exp.root";
    cfg.ptwNonPromptFor = "../compareDataToMC/ratioDataMC_AA_Btopsi2S_DATA_ctauCut_y1p6_2p4_260323_2exp.root";
    return true;
  }

  return false;
}

} // namespace

void plot_weight_vs_pt_compare(
    TString channel = "JPsi_pbpb",
    int state = 1,
    float cLow = 0, float cHigh = 180,
    int maxEvents = 10)
{
  gStyle->SetOptStat(0);

  ChannelConfig cfg;
  if (!BuildConfig(channel, cfg))
  {
    cout << "[ERROR] Unsupported channel: " << channel.Data() << endl;
    cout << "        Use one of: JPsi_pbpb, Jpsi_pp, psi2S_pp, psi2S_pbpb" << endl;
    return;
  }

  TString sampleTag = (state == 1) ? "PR" : "NP";
  TString inputMC = (state == 1) ? cfg.inputPrompt : cfg.inputNonPrompt;

  TFile *fSingle = nullptr;
  TFile *fMid = nullptr;
  TFile *fFor = nullptr;
  TF1 *fptwSingle = nullptr;
  TF1 *fptwMid = nullptr;
  TF1 *fptwFor = nullptr;

  if (!cfg.useRapSplitPtW)
  {
    TString fpath = (state == 1) ? cfg.ptwPrompt : cfg.ptwNonPrompt;
    fSingle = TFile::Open(fpath, "READ");
    if (!fSingle || fSingle->IsZombie())
    {
      cout << "[ERROR] Cannot open pt-weight file: " << fpath.Data() << endl;
      return;
    }
    fptwSingle = (TF1 *)fSingle->Get("dataMC_Ratio1");
    if (!fptwSingle)
    {
      cout << "[ERROR] Cannot find TF1 dataMC_Ratio1 in: " << fpath.Data() << endl;
      return;
    }
  }
  else
  {
    TString fpathMid = (state == 1) ? cfg.ptwPromptMid : cfg.ptwNonPromptMid;
    TString fpathFor = (state == 1) ? cfg.ptwPromptFor : cfg.ptwNonPromptFor;

    fMid = TFile::Open(fpathMid, "READ");
    fFor = TFile::Open(fpathFor, "READ");

    if (!fMid || fMid->IsZombie() || !fFor || fFor->IsZombie())
    {
      cout << "[ERROR] Cannot open rapidity-split pt-weight files:" << endl;
      cout << "        mid = " << fpathMid.Data() << endl;
      cout << "        for = " << fpathFor.Data() << endl;
      return;
    }

    fptwMid = (TF1 *)fMid->Get("dataMC_Ratio1");
    fptwFor = (TF1 *)fFor->Get("dataMC_Ratio1");

    if (!fptwMid || !fptwFor)
    {
      cout << "[ERROR] Missing TF1 dataMC_Ratio1 in rapidity-split pt-weight files." << endl;
      return;
    }
  }

  TChain *mytree = new TChain("hionia/myTree");
  mytree->Add(inputMC.Data());
  const int totalEntries = mytree->GetEntries();
  const int nevt = (maxEvents > 0 && maxEvents < totalEntries) ? maxEvents : totalEntries;
  cout << "[INFO] Channel=" << cfg.tag.Data() << ", state=" << sampleTag.Data() << endl;
  cout << "[INFO] Entries in tree: " << totalEntries << ", processing: " << nevt << endl;

  const int maxBranchSize = 1000;

  Int_t Centrality = 0;
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

  const int nPtBins = (cfg.tag == "Jpsi_pp") ? 50 : ((cfg.tag == "JPsi_pbpb") ? 47 : 37);

  TProfile *pGenNoPt = new TProfile("pGenNoPt", "GEN; p_{T} (GeV/c); <fill weight>", nPtBins, cfg.ptPlotMin, cfg.ptPlotMax);
  TProfile *pGenWithPt = new TProfile("pGenWithPt", "GEN; p_{T} (GeV/c); <fill weight>", nPtBins, cfg.ptPlotMin, cfg.ptPlotMax);
  TProfile *pRecoNoPt = new TProfile("pRecoNoPt", "RECO; p_{T} (GeV/c); <fill weight>", nPtBins, cfg.ptPlotMin, cfg.ptPlotMax);
  TProfile *pRecoWithPt = new TProfile("pRecoWithPt", "RECO; p_{T} (GeV/c); <fill weight>", nPtBins, cfg.ptPlotMin, cfg.ptPlotMax);

  double ptBin_for[6] = {0, 3.5, 6.5, 9, 12, 40};
  double ptBin_mid[8] = {0, 6.5, 9, 12, 15, 20, 25, 40};
  double centBin_for[5] = {0, 20, 60, 100, 180};
  double centBin_mid[7] = {0, 20, 40, 60, 80, 100, 180};

  TProfile *pPtWGenPtFor = new TProfile("pPtWGenPtFor", "GEN ptWeight; p_{T} (GeV/c); <ptWeight>", 5, ptBin_for);
  TProfile *pPtWGenPtMid = new TProfile("pPtWGenPtMid", "GEN ptWeight; p_{T} (GeV/c); <ptWeight>", 7, ptBin_mid);
  TProfile *pPtWRecoPtFor = new TProfile("pPtWRecoPtFor", "RECO ptWeight; p_{T} (GeV/c); <ptWeight>", 5, ptBin_for);
  TProfile *pPtWRecoPtMid = new TProfile("pPtWRecoPtMid", "RECO ptWeight; p_{T} (GeV/c); <ptWeight>", 7, ptBin_mid);
  TProfile *pPtWGenCentFor = new TProfile("pPtWGenCentFor", "GEN ptWeight; Centrality; <ptWeight>", 4, centBin_for);
  TProfile *pPtWGenCentMid = new TProfile("pPtWGenCentMid", "GEN ptWeight; Centrality; <ptWeight>", 6, centBin_mid);
  TProfile *pPtWRecoCentFor = new TProfile("pPtWRecoCentFor", "RECO ptWeight; Centrality; <ptWeight>", 4, centBin_for);
  TProfile *pPtWRecoCentMid = new TProfile("pPtWRecoCentMid", "RECO ptWeight; Centrality; <ptWeight>", 6, centBin_mid);

  TProfile *pBaseGenPtFor = new TProfile("pBaseGenPtFor", "GEN baseWeight; p_{T} (GeV/c); <Gen Weight>", 5, ptBin_for);
  TProfile *pBaseGenPtMid = new TProfile("pBaseGenPtMid", "GEN baseWeight; p_{T} (GeV/c); <Gen Weight>", 7, ptBin_mid);
  TProfile *pBaseRecoPtFor = new TProfile("pBaseRecoPtFor", "RECO baseWeight; p_{T} (GeV/c); <Gen Weight>", 5, ptBin_for);
  TProfile *pBaseRecoPtMid = new TProfile("pBaseRecoPtMid", "RECO baseWeight; p_{T} (GeV/c); <Gen Weight>", 7, ptBin_mid);
  TProfile *pBaseGenCentFor = new TProfile("pBaseGenCentFor", "GEN baseWeight; Centrality; <Gen Weight>", 4, centBin_for);
  TProfile *pBaseGenCentMid = new TProfile("pBaseGenCentMid", "GEN baseWeight; Centrality; <Gen Weight>", 6, centBin_mid);
  TProfile *pBaseRecoCentFor = new TProfile("pBaseRecoCentFor", "RECO baseWeight; Centrality; <Gen Weight>", 4, centBin_for);
  TProfile *pBaseRecoCentMid = new TProfile("pBaseRecoCentMid", "RECO baseWeight; Centrality; <Gen Weight>", 6, centBin_mid);

  TProfile *pWGenPtFor = new TProfile("pWGenPtFor", "GEN Weight; p_{T} (GeV/c); <Weight>", 5, ptBin_for);
  TProfile *pWGenPtMid = new TProfile("pWGenPtMid", "GEN Weight; p_{T} (GeV/c); <Weight>", 7, ptBin_mid);
  TProfile *pWRecoPtFor = new TProfile("pWRecoPtFor", "RECO Weight; p_{T} (GeV/c); <Weight>", 5, ptBin_for);
  TProfile *pWRecoPtMid = new TProfile("pWRecoPtMid", "RECO Weight; p_{T} (GeV/c); <Weight>", 7, ptBin_mid);
  TProfile *pWGenCentFor = new TProfile("pWGenCentFor", "GEN Weight; Centrality; <Weight>", 4, centBin_for);
  TProfile *pWGenCentMid = new TProfile("pWGenCentMid", "GEN Weight; Centrality; <Weight>", 6, centBin_mid);
  TProfile *pWRecoCentFor = new TProfile("pWRecoCentFor", "RECO Weight; Centrality; <Weight>", 4, centBin_for);
  TProfile *pWRecoCentMid = new TProfile("pWRecoCentMid", "RECO Weight; Centrality; <Weight>", 6, centBin_mid);

  const ULong64_t trigMask = ((ULong64_t)1 << cfg.kTrigSel);

  for (int iev = 0; iev < nevt; ++iev)
  {
    if (iev % 1000000 == 0)
      cout << "[INFO] EVENT " << iev << " / " << nevt << endl;

    mytree->GetEntry(iev);

    if (cfg.isPbPb && !(Centrality > cLow && Centrality < cHigh))
      continue;

    double baseWeight = 1.0;
    if (cfg.tag == "JPsi_pbpb" || cfg.tag == "psi2S_pbpb")
      baseWeight = findNcoll(Centrality) * Gen_weight;
    else if (cfg.tag == "psi2S_pp")
      baseWeight = Gen_weight;
    else if (cfg.tag == "Jpsi_pp")
      baseWeight = 1.0;

    const bool HLTPass = ((HLTriggers & trigMask) == trigMask);

    for (int igen = 0; igen < Gen_QQ_size; ++igen)
    {
      TLorentzVector *JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
      TLorentzVector *mupl_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      TLorentzVector *mumi_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      const double rap = fabs(JP_Gen->Rapidity());
      const double pt = JP_Gen->Pt();

      if (!(JP_Gen->M() > cfg.massLow && JP_Gen->M() < cfg.massHigh)) continue;
      if (!(pt < 50.0 && rap < 2.4)) continue;
      if (!(IsAcceptanceQQ(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptanceQQ(mumi_Gen->Pt(), fabs(mumi_Gen->Eta())))) continue;
      if (!(fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4)) continue;
      if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0) continue;

      if (!PassChannelPtRegion(cfg.tag, rap, pt)) continue;

      const double ptWeight = EvalPtWeight(cfg, fptwSingle, fptwMid, fptwFor, rap, pt);
      const double totalWeight = baseWeight * ptWeight;
      pGenNoPt->Fill(pt, baseWeight);
      pGenWithPt->Fill(pt, totalWeight);
      if (rap > 1.6 && rap < 2.4 && pt > 3.5 && pt < 40.0)
      {
        pPtWGenPtFor->Fill(pt, ptWeight);
        pBaseGenPtFor->Fill(pt, baseWeight);
        pWGenPtFor->Fill(pt, totalWeight);
        if (cfg.isPbPb) pPtWGenCentFor->Fill(Centrality, ptWeight);
        if (cfg.isPbPb) pBaseGenCentFor->Fill(Centrality, baseWeight);
        if (cfg.isPbPb) pWGenCentFor->Fill(Centrality, totalWeight);
      }
      if (rap < 1.6 && pt > 6.5 && pt < 40.0)
      {
        pPtWGenPtMid->Fill(pt, ptWeight);
        pBaseGenPtMid->Fill(pt, baseWeight);
        pWGenPtMid->Fill(pt, totalWeight);
        if (cfg.isPbPb) pPtWGenCentMid->Fill(Centrality, ptWeight);
        if (cfg.isPbPb) pBaseGenCentMid->Fill(Centrality, baseWeight);
        if (cfg.isPbPb) pWGenCentMid->Fill(Centrality, totalWeight);
      }
    }

    if (!HLTPass) continue;

    for (int irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      TLorentzVector *JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
      TLorentzVector *mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      TLorentzVector *mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      const bool HLTFilterPass = ((Reco_QQ_trig[irqq] & trigMask) == trigMask);
      if (!HLTFilterPass) continue;

      if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      if (Reco_QQ_whichGen[irqq] == -1) continue;
      if (Reco_QQ_sign[irqq] != 0) continue;
      if (Reco_QQ_VtxProb[irqq] < 0.01) continue;

      const double rap = fabs(JP_Reco->Rapidity());
      const double pt = JP_Reco->Pt();

      if (!(JP_Reco->M() > cfg.massLow && JP_Reco->M() < cfg.massHigh)) continue;
      if (!(pt < 50.0 && rap < 2.4)) continue;
      if (!(IsAcceptanceQQ(mupl_Reco->Pt(), mupl_Reco->Eta()) && IsAcceptanceQQ(mumi_Reco->Pt(), mumi_Reco->Eta()))) continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));
      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));
      if (!(passMuonTypePl && passMuonTypeMi)) continue;

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
      if (!(muplSoft && mumiSoft)) continue;

      if (!PassChannelPtRegion(cfg.tag, rap, pt)) continue;

      const double ptWeight = EvalPtWeight(cfg, fptwSingle, fptwMid, fptwFor, rap, pt);
      const double totalWeight = baseWeight * ptWeight;
      pRecoNoPt->Fill(pt, baseWeight);
      pRecoWithPt->Fill(pt, totalWeight);
      if (rap > 1.6 && rap < 2.4 && pt > 3.5 && pt < 40.0)
      {
        pPtWRecoPtFor->Fill(pt, ptWeight);
        pBaseRecoPtFor->Fill(pt, baseWeight);
        pWRecoPtFor->Fill(pt, totalWeight);
        if (cfg.isPbPb) pPtWRecoCentFor->Fill(Centrality, ptWeight);
        if (cfg.isPbPb) pBaseRecoCentFor->Fill(Centrality, baseWeight);
        if (cfg.isPbPb) pWRecoCentFor->Fill(Centrality, totalWeight);
      }
      if (rap < 1.6 && pt > 6.5 && pt < 40.0)
      {
        pPtWRecoPtMid->Fill(pt, ptWeight);
        pBaseRecoPtMid->Fill(pt, baseWeight);
        pWRecoPtMid->Fill(pt, totalWeight);
        if (cfg.isPbPb) pPtWRecoCentMid->Fill(Centrality, ptWeight);
        if (cfg.isPbPb) pBaseRecoCentMid->Fill(Centrality, baseWeight);
        if (cfg.isPbPb) pWRecoCentMid->Fill(Centrality, totalWeight);
      }
    }
  }

  PrintTripleProfileSummary("GEN values vs pt (forward bins)", pPtWGenPtFor, pBaseGenPtFor, pWGenPtFor);
  PrintTripleProfileSummary("GEN values vs pt (mid bins)", pPtWGenPtMid, pBaseGenPtMid, pWGenPtMid);
  PrintTripleProfileSummary("RECO values vs pt (forward bins)", pPtWRecoPtFor, pBaseRecoPtFor, pWRecoPtFor);
  PrintTripleProfileSummary("RECO values vs pt (mid bins)", pPtWRecoPtMid, pBaseRecoPtMid, pWRecoPtMid);
  if (cfg.isPbPb)
  {
    PrintTripleProfileSummary("GEN values vs centrality (forward bins)", pPtWGenCentFor, pBaseGenCentFor, pWGenCentFor);
    PrintTripleProfileSummary("GEN values vs centrality (mid bins)", pPtWGenCentMid, pBaseGenCentMid, pWGenCentMid);
    PrintTripleProfileSummary("RECO values vs centrality (forward bins)", pPtWRecoCentFor, pBaseRecoCentFor, pWRecoCentFor);
    PrintTripleProfileSummary("RECO values vs centrality (mid bins)", pPtWRecoCentMid, pBaseRecoCentMid, pWRecoCentMid);
  }

  TH1D *hGenNoPt = ProfileToHist(pGenNoPt, "hGenNoPt");
  TH1D *hGenWithPt = ProfileToHist(pGenWithPt, "hGenWithPt");
  TH1D *hRecoNoPt = ProfileToHist(pRecoNoPt, "hRecoNoPt");
  TH1D *hRecoWithPt = ProfileToHist(pRecoWithPt, "hRecoWithPt");

  TH1D *hRatioGen = MakeRatio(hGenWithPt, hGenNoPt, "hRatioGen");
  TH1D *hRatioReco = MakeRatio(hRecoWithPt, hRecoNoPt, "hRatioReco");

  hGenNoPt->SetLineColor(kBlue + 1);
  hGenNoPt->SetLineWidth(2);
  hGenWithPt->SetLineColor(kRed + 1);
  hGenWithPt->SetLineWidth(2);
  hRecoNoPt->SetLineColor(kBlue + 1);
  hRecoNoPt->SetLineWidth(2);
  hRecoWithPt->SetLineColor(kRed + 1);
  hRecoWithPt->SetLineWidth(2);

  hRatioGen->SetTitle("GEN ratio (with/without)");
  hRatioReco->SetTitle("RECO ratio (with/without)");
  hRatioGen->GetYaxis()->SetRangeUser(0.6, 1.4);
  hRatioReco->GetYaxis()->SetRangeUser(0.6, 1.4);
  hRatioGen->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatioReco->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatioGen->GetYaxis()->SetTitle("Ratio");
  hRatioReco->GetYaxis()->SetTitle("Ratio");

  gSystem->mkdir("./figs", kTRUE);
  gSystem->mkdir("./roots", kTRUE);

  TCanvas *c = new TCanvas("c_weight_vs_pt_compare", "c_weight_vs_pt_compare", 1400, 900);
  c->Divide(2, 2);

  TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(hGenNoPt, "without pt weight", "l");
  leg->AddEntry(hGenWithPt, "with pt weight", "l");

  c->cd(1);
  hGenNoPt->SetTitle("GEN mean fill weight vs p_{T}");
  hGenNoPt->GetYaxis()->SetTitle("<fill weight>");
  hGenNoPt->Draw("hist");
  hGenWithPt->Draw("hist same");
  leg->Draw();

  c->cd(2);
  hRatioGen->Draw("hist");

  c->cd(3);
  hRecoNoPt->SetTitle("RECO mean fill weight vs p_{T}");
  hRecoNoPt->GetYaxis()->SetTitle("<fill weight>");
  hRecoNoPt->Draw("hist");
  hRecoWithPt->Draw("hist same");
  leg->Draw();

  c->cd(4);
  hRatioReco->Draw("hist");

  TLatex lt;
  lt.SetNDC();
  lt.SetTextSize(0.03);
  c->cd(1);
  if (cfg.isPbPb)
    lt.DrawLatex(0.14, 0.86, Form("%s %s, cent %.0f-%.0f%%", (state == 1 ? "Prompt" : "Non-prompt"), cfg.tag.Data(), cLow / 2.0, cHigh / 2.0));
  else
    lt.DrawLatex(0.14, 0.86, Form("%s %s", (state == 1 ? "Prompt" : "Non-prompt"), cfg.tag.Data()));

  TString outBase;
  if (cfg.isPbPb)
    outBase = Form("fillWeight_vsPt_compare_%s_%s_cent%.0fto%.0f", cfg.tag.Data(), sampleTag.Data(), cLow, cHigh);
  else
    outBase = Form("fillWeight_vsPt_compare_%s_%s", cfg.tag.Data(), sampleTag.Data());

  TString outPng = Form("./figs/%s.pdf", outBase.Data());
  c->SaveAs(outPng);

  TString outRoot = Form("./roots/%s.root", outBase.Data());
  TFile *fout = new TFile(outRoot, "RECREATE");
  pGenNoPt->Write();
  pGenWithPt->Write();
  pRecoNoPt->Write();
  pRecoWithPt->Write();
  pPtWGenPtFor->Write();
  pPtWGenPtMid->Write();
  pPtWRecoPtFor->Write();
  pPtWRecoPtMid->Write();
  pPtWGenCentFor->Write();
  pPtWGenCentMid->Write();
  pPtWRecoCentFor->Write();
  pPtWRecoCentMid->Write();
  pBaseGenPtFor->Write();
  pBaseGenPtMid->Write();
  pBaseRecoPtFor->Write();
  pBaseRecoPtMid->Write();
  pBaseGenCentFor->Write();
  pBaseGenCentMid->Write();
  pBaseRecoCentFor->Write();
  pBaseRecoCentMid->Write();
  pWGenPtFor->Write();
  pWGenPtMid->Write();
  pWRecoPtFor->Write();
  pWRecoPtMid->Write();
  pWGenCentFor->Write();
  pWGenCentMid->Write();
  pWRecoCentFor->Write();
  pWRecoCentMid->Write();
  hGenNoPt->Write();
  hGenWithPt->Write();
  hRecoNoPt->Write();
  hRecoWithPt->Write();
  hRatioGen->Write();
  hRatioReco->Write();
  c->Write();
  fout->Write();
  fout->Close();

  SaveTripleProfilePlot(pPtWGenPtFor,
                        pBaseGenPtFor,
                        pWGenPtFor,
                        "GEN values vs p_{T} (forward bins)",
                        Form("./figs/%s_GEN_ptWeight_vs_pt_forward.pdf", outBase.Data()),
                        "p_{T} (GeV/c)");
  SaveTripleProfilePlot(pPtWGenPtMid,
                        pBaseGenPtMid,
                        pWGenPtMid,
                        "GEN values vs p_{T} (mid bins)",
                        Form("./figs/%s_GEN_ptWeight_vs_pt_mid.pdf", outBase.Data()),
                        "p_{T} (GeV/c)");
  SaveTripleProfilePlot(pPtWRecoPtFor,
                        pBaseRecoPtFor,
                        pWRecoPtFor,
                        "RECO values vs p_{T} (forward bins)",
                        Form("./figs/%s_RECO_ptWeight_vs_pt_forward.pdf", outBase.Data()),
                        "p_{T} (GeV/c)");
  SaveTripleProfilePlot(pPtWRecoPtMid,
                        pBaseRecoPtMid,
                        pWRecoPtMid,
                        "RECO values vs p_{T} (mid bins)",
                        Form("./figs/%s_RECO_ptWeight_vs_pt_mid.pdf", outBase.Data()),
                        "p_{T} (GeV/c)");

  if (cfg.isPbPb)
  {
    SaveTripleProfilePlot(pPtWGenCentFor,
                          pBaseGenCentFor,
                          pWGenCentFor,
                          "GEN values vs centrality (forward bins)",
                          Form("./figs/%s_GEN_ptWeight_vs_cent_forward.pdf", outBase.Data()),
                          "Centrality");
    SaveTripleProfilePlot(pPtWGenCentMid,
                          pBaseGenCentMid,
                          pWGenCentMid,
                          "GEN values vs centrality (mid bins)",
                          Form("./figs/%s_GEN_ptWeight_vs_cent_mid.pdf", outBase.Data()),
                          "Centrality");
    SaveTripleProfilePlot(pPtWRecoCentFor,
                          pBaseRecoCentFor,
                          pWRecoCentFor,
                          "RECO values vs centrality (forward bins)",
                          Form("./figs/%s_RECO_ptWeight_vs_cent_forward.pdf", outBase.Data()),
                          "Centrality");
    SaveTripleProfilePlot(pPtWRecoCentMid,
                          pBaseRecoCentMid,
                          pWRecoCentMid,
                          "RECO values vs centrality (mid bins)",
                          Form("./figs/%s_RECO_ptWeight_vs_cent_mid.pdf", outBase.Data()),
                          "Centrality");
  }

  cout << "[INFO] Saved figure : " << outPng.Data() << endl;
  cout << "[INFO] Saved ROOT   : " << outRoot.Data() << endl;
}
