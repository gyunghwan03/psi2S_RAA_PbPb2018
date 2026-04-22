#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tnp_weight_lowptPbPb_num_den_new.h"
#include <TAttMarker.h>

using namespace std;
// v4: compute ctau cuts from the current sample (quantile) instead of reading pre-made histograms.
void get_Eff_JPsi_pbpb_ctauCut_v4_backup(
    int state = 2,
    bool isTnP = true, int isPtWeight = 0, bool isMC = true,
    float ptLow = 3.0, float ptHigh = 50.0,
    float yLow = 0.0, float yHigh = 2.4,
    float cLow = 0, float cHigh = 180)
{

  gStyle->SetOptStat(0);
  int kTrigSel = 12; // jpsi=12,Upsilon=13

  float muPtCut = 0; // 3.5, 1.8

  TString ptSys;
  if (isPtWeight == 0)
    ptSys = "nomi";
  else if (isPtWeight == 1)
    ptSys = "up";
  else if (isPtWeight == -1)
    ptSys = "down";
  // jpsi
  float massLow = 0.0;
  float massHigh = 10.0;
  // float massLow = 2.6;
  // float massHigh = 3.5;

  double min = 0;
  double max = ptHigh;
  const int numBins = 9; // 50;//(max-min)/binwidth;  //31//

  TString inputPrompt = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
  TString inputNonPrompt = "/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root";

  // Trees: use prompt sample to define ctau cuts when running nonprompt.
  TChain *mytree = new TChain("hionia/myTree");
  TChain *mytreePrompt = nullptr;
  if (state == 2)
  {
    mytree->Add(inputNonPrompt.Data());
    mytreePrompt = new TChain("hionia/myTree");
    mytreePrompt->Add(inputPrompt.Data());
  }
  else
  {
    mytree->Add(inputPrompt.Data());
  }

  cout << "dmoon chk entries (target) : " << mytree->GetEntries() << endl;
  if (mytreePrompt)
    cout << "dmoon chk entries (prompt for cut) : " << mytreePrompt->GetEntries() << endl;

  // pT reweighting function
  TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
  TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root", "read");

  if (state == 2)
  {
    fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
    fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root", "read");
  }
  TF1 *fptw1 = (TF1 *)fPtW1->Get("dataMC_Ratio1");
  TF1 *fptw2 = (TF1 *)fPtW2->Get("dataMC_Ratio1");

  // double ptBin_all[16] = {0,6.5,7.5,8.5,9.5,11,13,15,17.5,20,22.5,25,27.5,30,40,50};
  double ptBin_all[18] = {0, 3., 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 22.5, 25, 27.5, 30, 50};
  // double ptBin_all[17] = {0,3.,4.5,5.5,6.5,7.5,8.5,9.5,11,13,15,17.5,20,22.5,25,27.5,50};
  double ptBin_for[6] = {0, 3.5, 6.5, 9, 12, 40};
  double ptBin_mid[8] = {0, 6.5, 9, 12, 15, 20, 25, 40};
  double centBin_all[16] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 180};
  double centBin_for[5] = {0, 20, 60, 100, 180};
  double centBin_mid[7] = {0, 20, 40, 60, 80, 100, 180};
  float yBin[7] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

  TH1D *hpt_reco_0 = new TH1D("hpt_reco_0", "hpt_reco_0", 17, ptBin_all);
  TH1D *hpt_reco_1 = new TH1D("hpt_reco_1", "hpt_reco_1", 5, ptBin_for);
  TH1D *hpt_reco_2 = new TH1D("hpt_reco_2", "hpt_reco_2", 7, ptBin_mid);

  TH1D *hpt_gen_0 = new TH1D("hpt_gen_0", "hpt_gen_0", 17, ptBin_all);
  TH1D *hpt_gen_1 = new TH1D("hpt_gen_1", "hpt_gen_1", 5, ptBin_for);
  TH1D *hpt_gen_2 = new TH1D("hpt_gen_2", "hpt_gen_2", 7, ptBin_mid);

  TH1D *hcent_reco_0 = new TH1D("hcent_reco_0", "hcent_reco_0", 15, centBin_all);
  TH1D *hcent_reco_1 = new TH1D("hcent_reco_1", "hcent_reco_1", 4, centBin_for);
  TH1D *hcent_reco_2 = new TH1D("hcent_reco_2", "hcent_reco_2", 6, centBin_mid);

  TH1D *hcent_gen_0 = new TH1D("hcent_gen_0", "hcent_gen_0", 15, centBin_all);
  TH1D *hcent_gen_1 = new TH1D("hcent_gen_1", "hcent_gen_1", 4, centBin_for);
  TH1D *hcent_gen_2 = new TH1D("hcent_gen_2", "hcent_gen_2", 6, centBin_mid);

  TH1D *hInt_gen_1 = new TH1D("hInt_gen_1", "", 1, 0, 50);
  TH1D *hInt_gen_2 = new TH1D("hInt_gen_2", "", 1, 0, 50);

  TH1D *hInt_reco_1 = new TH1D("hInt_reco_1", "", 1, 0, 50);
  TH1D *hInt_reco_2 = new TH1D("hInt_reco_2", "", 1, 0, 50);

  // Precompute and cache all needed l_cut values; filled later after ctau scan
  double l_cut_for_pt[4];
  double l_cut_mid_pt[6];
  double l_cut_for_cent[4];
  double l_cut_mid_cent[6];
  double l_cut_for_int = -1e6;
  double l_cut_mid_int = -1e6;

  struct Candidate
  {
    double pt;
    double rap;
    double cent;
    double ctau;
    double weight;
  };

  std::vector<std::pair<double, double>> ctau_for_pt[4];
  std::vector<std::pair<double, double>> ctau_mid_pt[6];
  std::vector<std::pair<double, double>> ctau_for_cent[4];
  std::vector<std::pair<double, double>> ctau_mid_cent[6];
  std::vector<std::pair<double, double>> ctau_for_int;
  std::vector<std::pair<double, double>> ctau_mid_int;


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// for systematics  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // parameters for all
  double p10 = 0.0, p11 = 0.0, p12 = 0.0, p13 = 0.0;
  // parameters for forward
  double p20 = 0.0, p21 = 0.0, p22 = 0.0, p23 = 0.0;
  // errors for all
  double e10 = 0.0, e11 = 0.0, e12 = 0.0, e13 = 0.0;
  // errors for forward
  double e20 = 0.0, e21 = 0.0, e22 = 0.0, e23 = 0.0;

  p10 = fptw1->GetParameter(0);
  p11 = fptw1->GetParameter(1);
  p12 = fptw1->GetParameter(2);
  p13 = fptw1->GetParameter(3);

  e10 = fptw1->GetParError(0);
  e11 = fptw1->GetParError(1);
  e12 = fptw1->GetParError(2);
  e13 = fptw1->GetParError(3);

  p20 = fptw2->GetParameter(0);
  p21 = fptw2->GetParameter(1);
  p22 = fptw2->GetParameter(2);
  p23 = fptw2->GetParameter(3);

  e20 = fptw2->GetParError(0);
  e21 = fptw2->GetParError(1);
  e22 = fptw2->GetParError(2);
  e23 = fptw2->GetParError(3);

  TF1 *fptw1sys = new TF1("fptw1sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 6.5, 50.0);
  TF1 *fptw2sys = new TF1("fptw2sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 3.0, 50.0);

  cout << "p10 : " << p10 << ", e10 : " << e10 << ", p20 : " << p20 << ", e20 : " << e20 << endl;

  // sys for nominal
  if (isPtWeight == 0)
  {
    p10 = p10;
    p11 = p11;
    p12 = p12;
    p13 = p13;

    p20 = p20;
    p21 = p21;
    p22 = p22;
    p23 = p23;
  }

  // sys for up
  if (isPtWeight == 1)
  {
    p10 = p10 + e10;
    p11 = p11 + e11;
    p12 = p12 + e12;
    p13 = p13 + e13;

    p20 = p20 + e10;
    p21 = p21 + e11;
    p22 = p22 + e12;
    p23 = p23 + e13;
  }

  // sys for down
  if (isPtWeight == -1)
  {
    p10 = p10 - e10;
    p11 = p11 - e11;
    p12 = p12 - e12;
    p13 = p13 - e13;

    p20 = p20 - e10;
    p21 = p21 - e11;
    p22 = p22 - e12;
    p23 = p23 - e13;
  }

  // cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

  fptw1sys->SetParameter(0, p10);
  fptw1sys->SetParameter(1, p11);
  fptw1sys->SetParameter(2, p12);
  fptw1sys->SetParameter(3, p13);

  fptw2sys->SetParameter(0, p20);
  fptw2sys->SetParameter(1, p21);
  fptw2sys->SetParameter(2, p22);
  fptw2sys->SetParameter(3, p23);

  TF1 *f1 = new TF1("f1", "[0]*TMath::Erf((x-[1])/[2])", 3, 50);
  f1->SetParameters(214, 22, 14);
  f1->SetParLimits(0, 0, 500);
  f1->SetParLimits(1, -1, 500);
  f1->SetParLimits(2, 0, 500);

  TH1D *hy_gen = new TH1D("hy_gen", "hy_gen", 6, yBin);
  TH1D *hy_reco = new TH1D("hy_reco", "hy_reco", 6, yBin);
  hy_gen->Sumw2();
  hy_reco->Sumw2();

  hpt_reco_0->Sumw2();
  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();

  hpt_gen_0->Sumw2();
  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();

  const int maxBranchSize = 1000;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Short_t Gen_QQ_size;
  Short_t Gen_mu_size;
  TClonesArray *Gen_QQ_4mom;
  TClonesArray *Gen_mu_4mom;
  ULong64_t Gen_QQ_trig[maxBranchSize];  //[Gen_QQ_size]
  Float_t Gen_QQ_VtxProb[maxBranchSize]; //[Gen_QQ_size]
  Short_t Gen_QQ_mupl_idx[maxBranchSize];
  Short_t Gen_QQ_mumi_idx[maxBranchSize];
  Short_t Gen_mu_charge[maxBranchSize];
  TBranch *b_Centrality; //! 
  TBranch *b_HLTriggers; //! 
  TBranch *b_Gen_QQ_size; //! 
  TBranch *b_Gen_mu_size; //! 
  TBranch *b_Gen_QQ_4mom; //! 
  TBranch *b_Gen_mu_4mom; //! 
  TBranch *b_Gen_QQ_trig; //! 
  TBranch *b_Gen_QQ_VtxProb; //! 
  TBranch *b_Gen_QQ_mupl_idx; //! 
  TBranch *b_Gen_QQ_mumi_idx; //! 
  TBranch *b_Gen_mu_charge; //!

  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;
  TClonesArray *Reco_QQ_4mom;
  TClonesArray *Reco_mu_4mom;
  ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
  ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_mu_size]
  Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
  Float_t Reco_QQ_ctau3D[maxBranchSize];  //[Reco_QQ_size]
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_whichGen[maxBranchSize];   //[Reco_QQ_size]
  Float_t Reco_mu_dxy[maxBranchSize];         //[Reco_mu_size]
  Float_t Reco_mu_dz[maxBranchSize];          //[Reco_mu_size]
  Int_t Reco_mu_nTrkWMea[maxBranchSize];      //[Reco_mu_size]
  Int_t Reco_mu_nPixWMea[maxBranchSize];      //[Reco_mu_size]
  Short_t Reco_QQ_sign[maxBranchSize];        //[Reco_QQ_size]
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Short_t Reco_mu_whichGen[maxBranchSize];
  Float_t Gen_weight;
  TBranch *b_Reco_QQ_size; //! 
  TBranch *b_Reco_mu_size; //! 
  TBranch *b_Reco_QQ_4mom; //! 
  TBranch *b_Reco_mu_4mom; //! 
  TBranch *b_Reco_QQ_trig; //! 
  TBranch *b_Reco_mu_trig; //! 
  TBranch *b_Reco_QQ_VtxProb; //! 
  TBranch *b_Reco_QQ_ctau3D; //! 
  TBranch *b_Reco_QQ_mupl_idx; //! 
  TBranch *b_Reco_QQ_mumi_idx; //! 
  TBranch *b_Reco_QQ_whichGen; //! 
  TBranch *b_Reco_mu_dxy; //! 
  TBranch *b_Reco_mu_dz; //! 
  TBranch *b_Reco_mu_nTrkWMea; //! 
  TBranch *b_Reco_mu_nPixWMea; //! 
  TBranch *b_Reco_QQ_sign; //! 
  TBranch *b_Reco_mu_SelectionType; //! 
  TBranch *b_Reco_mu_whichGen; //! 
  TBranch *b_Gen_weight; //!

  Gen_QQ_4mom = nullptr;
  Gen_mu_4mom = nullptr;
  Reco_QQ_4mom = nullptr;
  Reco_mu_4mom = nullptr;

  auto setBranches = [&](TChain *chain)
  {
    chain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    chain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
    chain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    chain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
    chain->SetBranchAddress("Gen_QQ_trig", Gen_QQ_trig, &b_Gen_QQ_trig);
    chain->SetBranchAddress("Gen_QQ_VtxProb", Gen_QQ_VtxProb, &b_Gen_QQ_VtxProb);
    chain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
    chain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
    chain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);

    chain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
    chain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
    chain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
    chain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
    chain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
    chain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
    chain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
    chain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
    chain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
    chain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
    chain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
    chain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

    chain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
    chain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
    chain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
    chain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
    chain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
    chain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
    chain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
    chain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
    chain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
  };

  setBranches(mytree);
  if (mytreePrompt)
    setBranches(mytreePrompt);

  TLorentzVector *JP_Gen = new TLorentzVector;
  TLorentzVector *mupl_Gen = new TLorentzVector;
  TLorentzVector *mumi_Gen = new TLorentzVector;

  TLorentzVector *JP_Reco = new TLorentzVector;
  TLorentzVector *mupl_Reco = new TLorentzVector;
  TLorentzVector *mumi_Reco = new TLorentzVector;

  double weight = 1;
  double tnp_weight = 1;
  double tnp_trig_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double pt_weight = 1;

  double pt, y, eta, pt1, pt2, eta1, eta2, phi1, phi2;
  double dcosMuTheta;
  double ptimb;
  //  TH2D* hpt_tnp_trig = new TH2D("hpt_tnp_trig","hpt_tnp_trig",numBins,ptBin,40,0,2);

  int kL2filter = 16; // jpsi=16,Upsilon=38
  int kL3filter = 17; // jpsi=17,Upsilon=39

  int count = 0;
  int counttnp = 0;
  int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;

  auto weightedQuantile = [&](const std::vector<std::pair<double, double>> &samples, double quantile, const std::string &label)
  {
    if (samples.empty())
    {
      cout << "[WARN] No ctau entries for " << label << endl;
      return -1e6;
    }
    std::vector<std::pair<double, double>> ordered = samples;
    std::sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    double totalW = 0.0;
    for (const auto &v : ordered)
      totalW += v.second;

    const double target = quantile * totalW;
    double running = 0.0;
    double value = ordered.back().first;
    for (const auto &v : ordered)
    {
      running += v.second;
      if (running >= target)
      {
        value = v.first;
        break;
      }
    }
    cout << "[l_cut] " << label << " = " << value << " (N=" << ordered.size() << ", sumW=" << totalW << ")" << endl;
    return value;
  };

  auto processCandidates = [&](TChain *chain, bool fillGen, const std::function<void(const Candidate &)> &consumer)
  {
    int nevt = chain->GetEntries();
    cout << "Total Events (loop target) : " << nevt << endl;
    int maxLoop = nevt; // cap for quick tests; set to nevt for full stats
    for (int iev = 0; iev < maxLoop && iev < nevt; ++iev)
    {
      if (iev % 100000 == 0)
        cout << ">>>>> EVENT " << iev << " / " << nevt << " (" << (int)(100. * iev / (double)nevt) << "%)" << endl;

      chain->GetEntry(iev);
      weight = findNcoll(Centrality) * Gen_weight;

      if (fillGen)
      {
        for (int igen = 0; fillGen && igen < Gen_QQ_size; igen++)
        {
          JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
          mupl_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
          mumi_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

          Double_t Rapidity_g = fabs(JP_Gen->Rapidity());
          if (!(JP_Gen->M() > massLow && JP_Gen->M() < massHigh))
            continue;

          if (!(JP_Gen->Pt() > 3 && JP_Gen->Pt() < 50 && fabs(JP_Gen->Rapidity()) < 2.4 && IsAcceptanceQQ(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptanceQQ(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
            continue;

          if (!(fabs(JP_Gen->Rapidity()) < 2.4 && fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
            continue;
          if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
            continue;

          pt_weight = 1;
          if (fabs(JP_Gen->Rapidity()) <= 1.6 && JP_Gen->Pt() >= 6.5)
            pt_weight = fptw1sys->Eval(JP_Gen->Pt());
          if (fabs(JP_Gen->Rapidity()) > 1.6 && fabs(JP_Gen->Rapidity()) < 2.4)
            pt_weight = fptw2sys->Eval(JP_Gen->Pt());
          if (JP_Gen->Pt() > 6.5)
            hy_gen->Fill(Rapidity_g, weight * pt_weight);
          if ((Centrality > cLow && Centrality < cHigh))
          {
            if (JP_Gen->Pt() > 3. && JP_Gen->Pt() < 6.5 && Rapidity_g > 1.6 && Rapidity_g < 2.4)
            {
              hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight);
            }
            if (JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 50 && Rapidity_g < 2.4)
            {
              hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight);
            }
            if (Rapidity_g > 1.6 && Rapidity_g < 2.4)
            {
              hpt_gen_1->Fill(JP_Gen->Pt(), weight * pt_weight);
            }
            if (Rapidity_g < 1.6)
            {
              hpt_gen_2->Fill(JP_Gen->Pt(), weight * pt_weight);
            }
          }
          if (Rapidity_g > 1.6 && Rapidity_g < 2.4)
          {
            hcent_gen_1->Fill(Centrality, weight * pt_weight);
            if (JP_Gen->Pt() > 3.5 && JP_Gen->Pt() < 40)
            {
              hInt_gen_1->Fill(1, weight * pt_weight);
            }
          }
          if (Rapidity_g <= 1.6)
          {
            hcent_gen_2->Fill(Centrality, weight * pt_weight);
            hInt_gen_2->Fill(1, weight * pt_weight);
          }
          if (Rapidity_g < 2.4 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 50)
          {
            hcent_gen_0->Fill(Centrality, weight * pt_weight);
          }
        }
      }

      bool HLTPass = false;
      if ((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
        HLTPass = true;

      for (Int_t irqq = 0; irqq < Reco_QQ_size; irqq++)
      {
        JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
        mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
        mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

        bool HLTFilterPass = false;
        if ((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
          HLTFilterPass = true;
        if (HLTFilterPass == false)
          continue;

        if (isMC)
        {
          if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
            continue;
          if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
            continue;
          if (Reco_QQ_whichGen[irqq] == -1)
            continue;
        }
        if (Reco_QQ_sign[irqq] != 0)
          continue;

        if (abs(JP_Reco->Rapidity()) > yHigh || abs(JP_Reco->Rapidity()) < yLow)
          continue;
        if (JP_Reco->Pt() < ptLow || JP_Reco->Pt() > ptHigh)
          continue;
        if (JP_Reco->M() < massLow || JP_Reco->M() > massHigh)
          continue;
        pt1 = mupl_Reco->Pt();
        pt2 = mumi_Reco->Pt();
        eta1 = mupl_Reco->Eta();
        eta2 = mumi_Reco->Eta();

        if (!(IsAcceptanceQQ(pt1, eta1) && IsAcceptanceQQ(pt2, eta2)))
          continue;

        bool passMuonTypePl = true;
        passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
        passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

        bool passMuonTypeMi = true;
        passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
        passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));
        if (passMuonTypePl == false || passMuonTypeMi == false)
          continue;

        phi1 = mupl_Reco->Phi();
        phi2 = mumi_Reco->Phi();

        TLorentzVector mudiff_Reco(*mupl_Reco - *mumi_Reco);
        dcosMuTheta = (mupl_Reco->P() * mupl_Reco->P() + mumi_Reco->P() * mumi_Reco->P() - mudiff_Reco.P() * mudiff_Reco.P()) / (2 * mupl_Reco->P() * mumi_Reco->P());
        ptimb = abs(pt1 - pt2) / (pt1 + pt2);

        bool muplSoft = ((Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) && (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) && (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) && (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) && passMuonTypePl);

        bool mumiSoft = ((Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) && (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) && (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) && (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) && passMuonTypeMi);

        if (!(muplSoft && mumiSoft))
          continue;

          if(isPtWeight){
        //pt_weight = f1->Eval(JP_Reco->Pt());
        if(abs(JP_Reco->Rapidity())>=1.6 && JP_Reco->Rapidity()<2.4) pt_weight = fptw1sys->Eval(JP_Reco->Pt());
        else if(abs(JP_Reco->Rapidity())<1.6 && JP_Reco->Pt()>=6.5) pt_weight = fptw2sys->Eval(JP_Reco->Pt());
        weight = weight * pt_weight;
      }

        if (isTnP)
        {
          tnp_weight = 1;
          tnp_trig_weight_mupl = -1;
          tnp_trig_weight_mumi = -1;
          tnp_weight = tnp_weight * tnp_weight_muid_pbpb(pt1, eta1, 0) * tnp_weight_muid_pbpb(pt2, eta2, 0);
          tnp_weight = tnp_weight * tnp_weight_trk_pbpb(eta1, 0) * tnp_weight_trk_pbpb(eta2, 0);

          if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
            continue;

          bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
          bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
          bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
          bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
          if (mupl_L2Filter == false || mumi_L2Filter == false)
          {
            cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl;
            cout << endl;
          }

          bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
          bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
          bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
          bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
          bool SelDone = false;

          if (mupl_isL2 && mumi_isL3)
          {
            tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
            tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
            SelDone = true;
          }
          else if (mupl_isL3 && mumi_isL2)
          {
            tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
            tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
            SelDone = true;
          }
          else if (mupl_isL3 && mumi_isL3)
          {
            int t[2] = {-1, 1};
            int l = rand() % (2);
            if (t[l] == -1)
            {
              tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
              tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
            }
            else if (t[l] == 1)
            {
              tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
              tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
            }
            SelDone = true;
          }
          if (SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1))
          {
            cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl;
            continue;
          }
          tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
        }

        pt_weight = 1;
        if (fabs(JP_Reco->Rapidity()) < 1.6 && JP_Reco->Pt() > 6.5)
          pt_weight = fptw1sys->Eval(JP_Reco->Pt());
        if (fabs(JP_Reco->Rapidity()) > 1.6 && fabs(JP_Reco->Rapidity()) < 2.4)
          pt_weight = fptw2sys->Eval(JP_Reco->Pt());

        Double_t Rapidity = fabs(JP_Reco->Rapidity());
        if (HLTPass == true && HLTFilterPass == true)
        {
          Candidate cand{JP_Reco->Pt(), Rapidity, static_cast<double>(Centrality), static_cast<double>(Reco_QQ_ctau3D[irqq]), weight * tnp_weight * pt_weight};
          consumer(cand);
        }
      }
    }
  };

  auto collectCtau = [&](const Candidate &cand)
  {
    if ((cand.cent > cLow && cand.cent < cHigh) && (cand.rap > 1.6 && cand.rap < 2.4))
    {
      for (int i = 0; i < 4; i++)
      {
        double ptLowBin = ptBin_for[i + 1];
        double ptHighBin = ptBin_for[i + 2];
        if (cand.pt > ptLowBin && cand.pt < ptHighBin)
          ctau_for_pt[i].push_back({cand.ctau, cand.weight});
      }
    }

    if ((cand.cent > cLow && cand.cent < cHigh) && (cand.rap < 1.6))
    {
      for (int i = 0; i < 6; i++)
      {
        double ptLowBin = ptBin_mid[i + 1];
        double ptHighBin = ptBin_mid[i + 2];
        if (cand.pt > ptLowBin && cand.pt < ptHighBin)
          ctau_mid_pt[i].push_back({cand.ctau, cand.weight});
      }
    }

    if (cand.rap > 1.6 && cand.rap < 2.4 && cand.pt > 3.5 && cand.pt < 40)
    {
      for (int i = 0; i < 4; i++)
      {
        double cLowBin = centBin_for[i];
        double cHighBin = centBin_for[i + 1];
        if (cand.cent > cLowBin && cand.cent < cHighBin)
          ctau_for_cent[i].push_back({cand.ctau, cand.weight});
      }
      if (cand.cent > 0 && cand.cent < 180)
        ctau_for_int.push_back({cand.ctau, cand.weight});
    }

    if (cand.rap < 1.6 && cand.pt > 6.5 && cand.pt < 40)
    {
      for (int i = 0; i < 6; i++)
      {
        double cLowBin = centBin_mid[i];
        double cHighBin = centBin_mid[i + 1];
        if (cand.cent > cLowBin && cand.cent < cHighBin)
          ctau_mid_cent[i].push_back({cand.ctau, cand.weight});
      }
      if (cand.cent > 0 && cand.cent < 180)
        ctau_mid_int.push_back({cand.ctau, cand.weight});
    }
  };

  TChain *chainCut = mytree;
  TChain *chainEval = mytree;
  if (state == 2 && mytreePrompt)
    chainCut = mytreePrompt;

  processCandidates(chainCut, false, collectCtau);

  // Prompt: keep 90% with ctau < l_cut.
  // Nonprompt: choose l_cut where prompt-like residual above cut would be 2%, i.e., 98th percentile.
  const double targetQuantile = (state == 1) ? 0.9 : 0.98;
  for (int i = 0; i < 4; i++)
    l_cut_for_pt[i] = weightedQuantile(ctau_for_pt[i], targetQuantile, Form("pt_for_bin_%d", i));
  for (int i = 0; i < 6; i++)
    l_cut_mid_pt[i] = weightedQuantile(ctau_mid_pt[i], targetQuantile, Form("pt_mid_bin_%d", i));
  for (int i = 0; i < 4; i++)
    l_cut_for_cent[i] = weightedQuantile(ctau_for_cent[i], targetQuantile, Form("cent_for_bin_%d", i));
  for (int i = 0; i < 6; i++)
    l_cut_mid_cent[i] = weightedQuantile(ctau_mid_cent[i], targetQuantile, Form("cent_mid_bin_%d", i));
  l_cut_for_int = weightedQuantile(ctau_for_int, targetQuantile, "cent_for_integrated");
  l_cut_mid_int = weightedQuantile(ctau_mid_int, targetQuantile, "cent_mid_integrated");

  auto fillHistograms = [&](const Candidate &cand)
  {
    if (cand.pt > 6.5)
      hy_reco->Fill(cand.rap, cand.weight);

    if ((cand.cent > cLow && cand.cent < cHigh))
    {
      if (cand.rap < 2.4 && cand.rap > 1.6 && cand.pt > 3. && cand.pt < 6.5)
        hpt_reco_0->Fill(cand.pt, cand.weight);
      if (cand.rap < 2.4 && cand.pt > 6.5 && cand.pt < 50)
        hpt_reco_0->Fill(cand.pt, cand.weight);

      if (cand.rap > 1.6 && cand.rap < 2.4)
      {
        for (int i = 0; i < 4; i++)
        {
          double ptLowBin = ptBin_for[i + 1];
          double ptHighBin = ptBin_for[i + 2];
          if (cand.pt > ptLowBin && cand.pt < ptHighBin)
          {
            double l_cut = l_cut_for_pt[i];
            if (state == 1)
            {
              if (cand.ctau < l_cut)
                hpt_reco_1->Fill(cand.pt, cand.weight);
            }
            else if (state == 2)
            {
              if (cand.ctau > l_cut)
                hpt_reco_1->Fill(cand.pt, cand.weight);
            }
          }
        }
      }

      if (cand.rap < 1.6)
      {
        for (int i = 0; i < 6; i++)
        {
          double ptLowBin = ptBin_mid[i + 1];
          double ptHighBin = ptBin_mid[i + 2];
          if (cand.pt > ptLowBin && cand.pt < ptHighBin)
          {
            double l_cut = l_cut_mid_pt[i];
            if (state == 1)
            {
              if (cand.ctau < l_cut)
                hpt_reco_2->Fill(cand.pt, cand.weight);
            }
            else if (state == 2)
            {
              if (cand.ctau > l_cut)
                hpt_reco_2->Fill(cand.pt, cand.weight);
            }
          }
        }
      }
    }

    if (cand.rap > 1.6 && cand.rap < 2.4)
    {
      if (cand.pt > 3.5 && cand.pt < 40)
      {
        for (int i = 0; i < 4; i++)
        {
          double cLowBin = centBin_for[i];
          double cHighBin = centBin_for[i + 1];
          if (cand.cent > cLowBin && cand.cent < cHighBin)
          {
            double l_cut = l_cut_for_cent[i];
            if (state == 1)
            {
              if (cand.ctau < l_cut)
                hcent_reco_1->Fill(cand.cent, cand.weight);
            }
            else if (state == 2)
            {
              if (cand.ctau > l_cut)
                hcent_reco_1->Fill(cand.cent, cand.weight);
            }
          }
        }
        if (cand.cent > 0 && cand.cent < 180)
        {
          double l_cut = l_cut_for_int;
          if (state == 1)
          {
            if (cand.ctau < l_cut)
              hInt_reco_1->Fill(1, cand.weight);
          }
          else if (state == 2)
          {
            if (cand.ctau > l_cut)
              hInt_reco_1->Fill(1, cand.weight);
          }
        }
      }
    }
    if (cand.rap < 1.6 && cand.pt > 6.5 && cand.pt < 40)
    {
      for (int i = 0; i < 6; i++)
      {
        double cLowBin = centBin_mid[i];
        double cHighBin = centBin_mid[i + 1];
        if (cand.cent > cLowBin && cand.cent < cHighBin)
        {
          double l_cut = l_cut_mid_cent[i];
          if (state == 1)
          {
            if (cand.ctau < l_cut)
              hcent_reco_2->Fill(cand.cent, cand.weight);
          }
          else if (state == 2)
          {
            if (cand.ctau > l_cut)
              hcent_reco_2->Fill(cand.cent, cand.weight);
          }
        }
      }
      if (cand.cent > 0 && cand.cent < 180)
      {
        double l_cut = l_cut_mid_int;
        if (state == 1)
        {
          if (cand.ctau < l_cut)
            hInt_reco_2->Fill(1, cand.weight);
        }
        else if (state == 2)
        {
          if (cand.ctau > l_cut)
            hInt_reco_2->Fill(1, cand.weight);
        }
      }
    }
    if (cand.rap < 2.4 && cand.pt > 6.5 && cand.pt < 50)
    {
      hcent_reco_0->Fill(cand.cent, cand.weight);
    }
  };

  processCandidates(chainEval, true, fillHistograms);

  //cout << "count " << count << endl;
  //cout << "counttnp " << counttnp << endl;

  // Divide
  TH1D *hpt_eff_0;
  TH1D *hpt_eff_1;
  TH1D *hpt_eff_2;

  TH1D *hcent_eff_0;
  TH1D *hcent_eff_1;
  TH1D *hcent_eff_2;

  TH1D *hInt_eff_1;
  TH1D *hInt_eff_2;

  hpt_eff_0 = (TH1D *)hpt_reco_0->Clone("hpt_eff_0");
  hpt_eff_1 = (TH1D *)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D *)hpt_reco_2->Clone("hpt_eff_2");

  hpt_eff_0->Divide(hpt_eff_0, hpt_gen_0, 1, 1, "B");
  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1, "B");
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1, "B");

  hcent_eff_0 = (TH1D *)hcent_reco_0->Clone("hcent_eff_0");
  hcent_eff_1 = (TH1D *)hcent_reco_1->Clone("hcent_eff_1");
  hcent_eff_2 = (TH1D *)hcent_reco_2->Clone("hcent_eff_2");

  hcent_eff_0->Divide(hcent_eff_0, hcent_gen_0, 1, 1, "B");
  hcent_eff_1->Divide(hcent_eff_1, hcent_gen_1, 1, 1, "B");
  hcent_eff_2->Divide(hcent_eff_2, hcent_gen_2, 1, 1, "B");

  hInt_eff_1 = (TH1D *)hInt_reco_1->Clone("hInt_eff_1");
  hInt_eff_2 = (TH1D *)hInt_reco_2->Clone("hInt_eff_2");

  hInt_eff_1->Divide(hInt_eff_1, hInt_gen_1, 1, 1, "B");
  hInt_eff_2->Divide(hInt_eff_2, hInt_gen_2, 1, 1, "B");

  TH1D *hy_eff;
  hy_eff = (TH1D *)hy_reco->Clone("hy_eff");
  hy_eff->Divide(hy_eff, hy_gen, 1, 1, "B");

  hpt_eff_0->SetTitle("Eff: Rapidity 0.0-2.4");
  hpt_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hpt_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");

  hcent_eff_0->SetTitle("Eff: Rapidity 0.0-2.4");
  hcent_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hcent_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");

  hInt_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hInt_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");
  /*
  f1->SetLineColor(kBlack);
  f1->SetLineWidth(2);
  hpt_eff_1->Fit(f1);
  f1->SetLineColor(kRed+1);
  f1->SetLineWidth(2);
  hpt_eff_2->Fit(f1);
  f1->SetLineColor(kOrange+1);
  f1->SetLineWidth(2);
  hpt_eff_3->Fit(f1);
  f1->SetLineColor(kGreen+1);
  f1->SetLineWidth(2);
  hpt_eff_4->Fit(f1);
  f1->SetLineColor(kBlue+1);
  f1->SetLineWidth(2);
  hpt_eff_5->Fit(f1);
  f1->SetLineColor(kViolet+1);
  f1->SetLineWidth(2);
  hpt_eff_6->Fit(f1);
  */

  // draw same
  // gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend = new TLegend(0.6, 0.84);
  TLatex *lt2 = new TLatex();
  lt2->SetNDC();
  lt2->SetTextSize(0.03);
  // auto legend2 = new TLegend(0.6,0.84);

  TLatex *lt3 = new TLatex();
  lt3->SetNDC();
  lt3->SetTextSize(0.03);
  auto legend2 = new TLegend(0.6, 0.84);

  gStyle->SetOptFit(0);

  TCanvas *cpt_eff_0 = new TCanvas("cpt_eff_0", "cpt_eff_0", 0, 0, 900, 800);
  cpt_eff_0->cd();
  hpt_eff_0->SetMarkerStyle(24);
  hpt_eff_0->SetMarkerColor(1);
  hpt_eff_0->SetLineColor(1);
  hpt_eff_0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hpt_eff_0->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_0->GetYaxis()->SetRangeUser(0., 0.9);
  hpt_eff_0->Draw("E");
  lt3->SetTextSize(0.03);
  lt3->DrawLatex(0.13, 0.84, "Prompt #Jpsi (PbPb)");
  legend2->AddEntry("hpt_eff_0", Form("|y|: 0.0-2.4, %0.0f-%0.0f %%", cLow / 2, cHigh / 2), "lep");
  legend2->SetBorderSize(0);
  legend2->Draw("E");
  cpt_eff_0->SaveAs("./pt50/figs/pbpb/Eff_pt_pbpb_0to2p4.png");

  TCanvas *cpt_eff = new TCanvas("cpt_eff", "cpt_eff", 0, 0, 900, 800);
  cpt_eff->cd();
  hpt_eff_1->SetMarkerStyle(24);
  hpt_eff_1->SetMarkerColor(1);
  hpt_eff_1->SetLineColor(1);
  hpt_eff_1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hpt_eff_1->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_1->GetYaxis()->SetRangeUser(0., 1.3);
  hpt_eff_1->Draw("E");

  // eeeeeee
  // hpt_eff_2->SetMarkerStyle(24);
  // hpt_eff_2->SetMarkerColor(1);
  // hpt_eff_2->SetLineColor(1);
  // hpt_eff_2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  // hpt_eff_2->GetYaxis()->SetTitle("Efficiency");
  // hpt_eff_2->GetYaxis()->SetRangeUser(0.,1.3);
  // hpt_eff_2->Draw("E");
  // eeeeeee
  lt1->SetTextSize(0.03);
  if (state == 1)
    lt1->DrawLatex(0.13, 0.84, "Prompt J/#Psi (PbPb)");
  else if (state == 2)
    lt1->DrawLatex(0.13, 0.84, "Non-Prompt J/#Psi (PbPb)");
  // drawsame
  hpt_eff_2->SetMarkerStyle(25);
  hpt_eff_2->SetMarkerColor(kRed + 1);
  hpt_eff_2->SetLineColor(kRed + 1);
  hpt_eff_2->Draw("same");
  legend->AddEntry("hpt_eff_1", Form("|y|: 1.6-2.4, %0.0f-%0.0f %%", cLow / 2, cHigh / 2), "lep");
  legend->AddEntry("hpt_eff_2", Form("|y|: 0-1.6, %0.0f-%0.0f %%", cLow / 2, cHigh / 2), "lep");
  legend->SetBorderSize(0);
  legend->Draw("E");
  cpt_eff->SaveAs("Eff_pt_noweight.png");

  TCanvas *ccent_eff = new TCanvas("ccent_eff", "ccent_eff", 0, 0, 900, 800);
  ccent_eff->cd();
  hcent_eff_1->SetMarkerStyle(24);
  hcent_eff_1->SetMarkerColor(1);
  hcent_eff_1->SetLineColor(1);
  hcent_eff_1->GetXaxis()->SetTitle("Centrality [%]");
  hcent_eff_1->GetYaxis()->SetTitle("Efficiency");
  hcent_eff_1->GetYaxis()->SetRangeUser(0., 1.3);
  hcent_eff_1->Draw("E");
  hcent_eff_2->SetMarkerStyle(25);
  hcent_eff_2->SetMarkerColor(kRed + 1);
  hcent_eff_2->SetLineColor(kRed + 1);
  hcent_eff_2->Draw("same");
  hcent_eff_0->SetMarkerStyle(26);
  hcent_eff_0->SetMarkerColor(kBlue + 1);
  hcent_eff_0->SetLineColor(kBlue + 1);
  hcent_eff_0->Draw("same");

  lt2->SetTextSize(0.03);
  if (state == 1)
    lt2->DrawLatex(0.13, 0.84, "Prompt J/#Psi (PbPb)");
  else if (state == 2)
    lt2->DrawLatex(0.13, 0.84, "Non-Prompt J/#Psi (PbPb)");

  legend->AddEntry("hcent_eff_1", Form("|y|: 1.6-2.4, p_{T} : 3-50 GeV/c"), "lep");
  legend->AddEntry("hcent_eff_2", Form("|y|: 0-1.6, p_{T} : 6.5-50 GeV/c"), "lep");
  legend->AddEntry("hcent_eff_0", Form("|y|: 0-2.4, p_{T} : 6.5-50 GeV/c"), "lep");
  legend->SetBorderSize(0);
  legend->Draw("E");
  ccent_eff->SaveAs("./figs/Eff_cent_noweight.png");

  TCanvas *cy_eff = new TCanvas("cy_eff", "cy_eff", 0, 0, 900, 800);
  cy_eff->cd();
  hy_eff->SetMarkerStyle(24);
  hy_eff->SetMarkerColor(1);
  hy_eff->SetLineColor(1);
  hy_eff->GetXaxis()->SetTitle("|y|");
  hy_eff->GetYaxis()->SetTitle("Efficiency");
  hy_eff->GetYaxis()->SetRangeUser(0., 1.2);
  hy_eff->Draw("E");
  lt1->SetTextSize(0.03);
  if (state == 1)
    lt1->DrawLatex(0.13, 0.84, "Prompt J/#Psi (PbPb)");
  else if (state == 2)
    lt1->DrawLatex(0.13, 0.84, "Non-Prompt J/#Psi (PbPb)");
  cy_eff->SaveAs("./figs/Eff_absy_noweight.png");

  // Save efficiency files for later use.

  hpt_eff_0->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_2p4", isTnP, 1, cLow, cHigh));
  hpt_eff_1->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4", isTnP, 1, cLow, cHigh));
  hpt_eff_2->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6  ", isTnP, 1, cLow, cHigh));
  hcent_eff_0->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_6p5_to_50_absy0_2p4", isTnP, 1));
  hcent_eff_1->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_3_to_40_absy1p6_2p4", isTnP, 1));
  hcent_eff_2->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_6p5_to_40_absy0_1p6  ", isTnP, 1));
  hy_eff->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d", isTnP, 1));
  hInt_eff_1->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4", isTnP, 1, cLow, cHigh));
  hInt_eff_2->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6  ", isTnP, 1, cLow, cHigh));

  // TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_drawsame1.root",cLow,cHigh,isPtWeight,isTnP);
  //TString outFileName = Form("./roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260105.root", cLow, cHigh, ptSys.Data(), isTnP);
  TString outFileName = Form("./roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260106_test.root", cLow, cHigh, ptSys.Data(), isTnP);
  if (state == 2)
    outFileName = Form("./roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_nprompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260106_test.root", cLow, cHigh, ptSys.Data(), isTnP);
  TFile *outFile = new TFile(outFileName, "RECREATE");
  hpt_eff_0->Write();
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hcent_eff_0->Write();
  hcent_eff_1->Write();
  hcent_eff_2->Write();
  hInt_eff_1->Write();
  hInt_eff_2->Write();
  hy_eff->Write();
  cpt_eff->Write();
  // if(isTnP) hpt_tnp_trig->Write();

  hpt_reco_0->Write();
  hpt_reco_1->Write();
  hpt_reco_2->Write();

  hpt_gen_0->Write();
  hpt_gen_1->Write();
  hpt_gen_2->Write();

  hcent_reco_0->Write();
  hcent_reco_1->Write();
  hcent_reco_2->Write();

  hcent_gen_0->Write();
  hcent_gen_1->Write();
  hcent_gen_2->Write();

  hInt_reco_1->Write();
  hInt_reco_2->Write();

  hInt_gen_1->Write();
  hInt_gen_2->Write();

  hy_gen->Write();
  hy_reco->Write();

  outFile->Close();

  cout << "Save file : " << outFileName << endl;
}


