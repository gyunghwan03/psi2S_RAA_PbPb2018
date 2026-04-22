#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
// #include "Style_jaebeom.h"
// #include "tnp_weight_lowptPbPb.h"
#include "../tnp_weight_lowptPbPb_num_den_new.h"
#include "../cutsAndBin.h"
#include <TAttMarker.h>
#include <TSystem.h>

using namespace std;
// v2 : change efficiency to include ctau cut
// v3 : trying to reduce time
// v4 : ctau Cut from decay_length_OniaTree_v2_PbPb_1S.C
// v5 : use single pT reweight function from all-rapidity file
// v6 : use rapidity-split pT reweight functions (mid/forward) like pp J/psi macro
// v7 : align centrality/reco selection with Gwak efficiency_1d baseline

static double EvalPtWeight(bool applyPtWeight, double rap, double pt, TF1* fptwMid, TF1* fptwFor)
{
  if(!applyPtWeight) return 1.0;
  double absRap = fabs(rap);
  if(absRap < 1.6 && fptwMid) return fptwMid->Eval(pt);
  if(absRap >= 1.6 && absRap < 2.4 && fptwFor) return fptwFor->Eval(pt);
  return 1.0;
}

static TF1* BuildPtWeightSys(TF1* src, const char* name, int isPtWeight)
{
  if(!src) return nullptr;
  TF1 *f = static_cast<TF1 *>(src->Clone(name));
  f->SetName(name);
  if(isPtWeight == 0) return f;

  const int nPar = f->GetNpar();
  for(int i = 0; i < nPar; ++i)
  {
    const double p = f->GetParameter(i);
    const double e = f->GetParError(i);
    if(isPtWeight == 1) f->SetParameter(i, p + e);
    else if(isPtWeight == -1) f->SetParameter(i, p - e);
  }
  return f;
}

static TString ResolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}
	
void get_Eff_JPsi_pbpb_ctauCut_v7(
    int state = 2,
    bool applyPtWeight = true,
    bool isTnP = true, int isPtWeight = 0, // 0 : nominal , 1 : up, -1 : down
    bool isMC = true,
    // bool isTnP = true, bool isPtWeight = false,
    float ptLow = 3.0, float ptHigh = 50.0,
    float yLow = 0.0, float yHigh = 2.4,
    float cLow = 0, float cHigh = 180,
    int maxEvents = -1,
    int reportEvery = 1000000
    // bool isTnP = false, bool isPtWeight = false, int state=1
    // bool isTnP = true, bool isPtWeight = false, int state=1
    // bool isTnP = false, bool isPtWeight = true, int state=1
)
{
  const TString baseDir = ResolveBaseDir();

  const int nCores = 28;
  ROOT::EnableImplicitMT(nCores);
  cout << "[INFO] ROOT ImplicitMT enabled with " << nCores
       << " threads (I/O/decompression-level parallelism)." << endl;
  cout << "[INFO] This macro still runs an explicit GetEntry event loop; "
          "ImplicitMT does not fully parallelize the event loop itself." << endl;

  gStyle->SetOptStat(0);
  int kTrigSel = 12; // jpsi=12,Upsilon=13

  float muPtCut = 0; // 3.5, 1.8

  TString fname;
  if(state==1) fname = "PR";
  else if(state==2) fname = "NP";

  TString ptSys;
  if (isPtWeight == 0)
    ptSys = "nomi";
  else if (isPtWeight == 1)
    ptSys = "up";
  else if (isPtWeight == -1)
    ptSys = "down";
  // jpsi
  //float massLow = 0.0;
  //float massHigh = 10.0;
  float massLow = 2.9;
  float massHigh = 3.3;

  double min = 0;
  double max = ptHigh;
  const int numBins = 9; // 50;//(max-min)/binwidth;  //31//

  // input files
  // PbPb
  TString inputMC1 = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
  if (state == 2)
    inputMC1 = "/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root"; // PbPb_non prompt
  TChain *mytree = new TChain("hionia/myTree");
  mytree->Add(inputMC1.Data());
  mytree->SetCacheSize(128 * 1024 * 1024);
  mytree->SetCacheLearnEntries(20);
  // TFile *inf = new TFile("../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root","READ");
  // TTree *mytree = (TTree*)inf->Get("myTree");
  // if(state==1) mytree->Add(inputMC2.Data());
  cout << "dmoon chk entries : " << mytree->GetEntries() << endl;
  // TString outFileName = "mc_eff_vs_pt_cent_prompt_pbpb_Jpsi.root";
  // if(state==2) outFileName = "mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi.root";

  // pT reweighting function
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y1p6_2p4_241016.root", "read");
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y0_1p6_241008.root", "read");
  // ratioDataMC_AA_Jpsi_DATA_y0.0-2.4_210915.root
  // TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_250116.root","read"); //v2
  // before
  // TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_241126.root", "read"); // v4

  // affter 250121
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_All_y_211218.root", "read");
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_241126.root", "read"); // v4
  // Jpsi v2 weighting function
  //TFile *fPtW1 = new TFile("./ratioDataMC_AA_Jpsi_DATA_All_y_211218.root", "read");
  //TFile *fPtW2 = new TFile("./ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root", "read");
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Psi2S_DATA_y0_1p6_230521.root","read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_PbPb_Psi2S_DATA_y1p6_2p4_230520.root","read");
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y0_1p6_230321.root", "read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y1p6_2p4_230621.root", "read");
  // PbPb J/psi weight function (rapidity split, same location priority as pp macro)
  TFile *fPtWMid = nullptr;
  TFile *fPtWFor = nullptr;
  if(state==1){
    fPtWMid = TFile::Open("/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260410.root","READ");
    if(!fPtWMid || fPtWMid->IsZombie()){
      fPtWMid = TFile::Open("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310.root","READ");
    }
    fPtWFor = TFile::Open("/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_260410.root","READ");
    if(!fPtWFor || fPtWFor->IsZombie()){
      fPtWFor = TFile::Open("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root","READ");
    }
  }
  else if(state==2){
    fPtWMid = TFile::Open("/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260410.root","READ");
    if(!fPtWMid || fPtWMid->IsZombie()){
      fPtWMid = TFile::Open("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310.root","READ");
    }
    fPtWFor = TFile::Open("/data/users/hanseul/jpsi_2603/compareDataToMC/roots/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_260410.root","READ");
    if(!fPtWFor || fPtWFor->IsZombie()){
      fPtWFor = TFile::Open("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root","READ");
    }
  }

  // if (state == 2) {
  //  TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoPsi2S_DATA_y0_1p6_240515.root", "read");
  //  TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoPsi2S_DATA_y1p6_2p4_240515.root", "read");
  //}
  TF1 *fptwMid = (fPtWMid && !fPtWMid->IsZombie()) ? (TF1*)fPtWMid->Get("dataMC_Ratio1") : nullptr;
  TF1 *fptwFor = (fPtWFor && !fPtWFor->IsZombie()) ? (TF1*)fPtWFor->Get("dataMC_Ratio1") : nullptr;

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
  TH1D *hpt_reco_0_noW = new TH1D("hpt_reco_0_noW", "hpt_reco_0_noW", 17, ptBin_all);
  TH1D *hpt_reco_1_noW = new TH1D("hpt_reco_1_noW", "hpt_reco_1_noW", 5, ptBin_for);
  TH1D *hpt_reco_2_noW = new TH1D("hpt_reco_2_noW", "hpt_reco_2_noW", 7, ptBin_mid);

  TH1D *hpt_gen_0_noW = new TH1D("hpt_gen_0_noW", "hpt_gen_0_noW", 17, ptBin_all);
  TH1D *hpt_gen_1_noW = new TH1D("hpt_gen_1_noW", "hpt_gen_1_noW", 5, ptBin_for);
  TH1D *hpt_gen_2_noW = new TH1D("hpt_gen_2_noW", "hpt_gen_2_noW", 7, ptBin_mid);

  TH1D *hInt_gen_1_noW = new TH1D("hInt_gen_1_noW", "", 1, 0, 50);
  TH1D *hInt_gen_2_noW = new TH1D("hInt_gen_2_noW", "", 1, 0, 50);
  TH1D *hInt_reco_1_noW = new TH1D("hInt_reco_1_noW", "", 1, 0, 50);
  TH1D *hInt_reco_2_noW = new TH1D("hInt_reco_2_noW", "", 1, 0, 50);

  // Precompute ctau cut values for forward pt bins to avoid opening files inside event loop
  // --- Precompute and cache all needed l_cut values ---
  double l_cut_for_pt[5];
  double l_cut_mid_pt[7];

  TString ctauCutFile = Form("%s/roots_1S_PbPb/ctau3D_cut_ptBin_%sMC_y0.0-2.4.root", baseDir.Data(), fname.Data());
  TFile *f_ctau = TFile::Open(ctauCutFile, "READ");
  TH1D *h_Lcut_for = nullptr;
  TH1D *h_Lcut_mid = nullptr;
  if(f_ctau && !f_ctau->IsZombie()){
    h_Lcut_for = (TH1D*)f_ctau->Get("hpt_ctau_1");
    h_Lcut_mid = (TH1D*)f_ctau->Get("hpt_ctau_2");
  } else {
    cout << "[WARN] ctau cut file not found: " << ctauCutFile << endl;
  }

  for(int i=0; i<5; ++i){
    if(h_Lcut_for) l_cut_for_pt[i] = h_Lcut_for->GetBinContent(i+1);
    else l_cut_for_pt[i] = -1e6;
  }
  cout << "L_cut forward pt bins: ";
  for(int i=0; i<5; ++i) cout << l_cut_for_pt[i] << " ";
  cout << endl;
  cout << "Mapping: l_cut_for_pt[0]=" << l_cut_for_pt[0] << " (pt 0-3.5)" << endl;
  cout << "         l_cut_for_pt[1]=" << l_cut_for_pt[1] << " (pt 3.5-6.5)" << endl;
  cout << "         l_cut_for_pt[2]=" << l_cut_for_pt[2] << " (pt 6.5-9)" << endl;
  cout << "         l_cut_for_pt[3]=" << l_cut_for_pt[3] << " (pt 9-12)" << endl;
  cout << "         l_cut_for_pt[4]=" << l_cut_for_pt[4] << " (pt 12-40)" << endl;

  for(int i=0; i<7; ++i){
    if(h_Lcut_mid) l_cut_mid_pt[i] = h_Lcut_mid->GetBinContent(i+1);
    else l_cut_mid_pt[i] = -1e6;
  }
  if(f_ctau){
    f_ctau->Close();
    delete f_ctau;
  }
  cout << "L_cut mid pt bins: ";
  for(int i=0; i<7; ++i) cout << l_cut_mid_pt[i] << " ";
  cout << endl;

  // v7: centrality efficiency uses the same pt-binned ctau criteria as pT efficiency (Gwak baseline)

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// for systematics  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  TF1 *fptwsysMid = BuildPtWeightSys(fptwMid, "fptwsysMid", isPtWeight);
  TF1 *fptwsysFor = BuildPtWeightSys(fptwFor, "fptwsysFor", isPtWeight);

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
  TBranch *b_Centrality;                 //!
  TBranch *b_HLTriggers;                 //!
  TBranch *b_Gen_QQ_size;                //!
  TBranch *b_Gen_mu_size;                //!
  TBranch *b_Gen_QQ_4mom;                //!
  TBranch *b_Gen_mu_4mom;                //!
  TBranch *b_Gen_QQ_trig;                //!
  TBranch *b_Gen_QQ_VtxProb;             //!

  Gen_QQ_4mom = 0;
  Gen_mu_4mom = 0;
  mytree->SetBranchStatus("*", 0);
  mytree->SetBranchStatus("Gen_QQ_size", 1);
  mytree->SetBranchStatus("Gen_mu_size", 1);
  mytree->SetBranchStatus("Gen_QQ_4mom", 1);
  mytree->SetBranchStatus("Gen_mu_4mom", 1);
  mytree->SetBranchStatus("Gen_QQ_mupl_idx", 1);
  mytree->SetBranchStatus("Gen_QQ_mumi_idx", 1);
  mytree->SetBranchStatus("Gen_mu_charge", 1);
  mytree->SetBranchStatus("Centrality", 1);
  mytree->SetBranchStatus("HLTriggers", 1);
  mytree->SetBranchStatus("Reco_QQ_size", 1);
  mytree->SetBranchStatus("Reco_mu_size", 1);
  mytree->SetBranchStatus("Reco_QQ_4mom", 1);
  mytree->SetBranchStatus("Reco_mu_4mom", 1);
  mytree->SetBranchStatus("Reco_QQ_trig", 1);
  mytree->SetBranchStatus("Reco_mu_trig", 1);
  mytree->SetBranchStatus("Reco_QQ_VtxProb", 1);
  mytree->SetBranchStatus("Reco_QQ_ctau3D", 1);
  mytree->SetBranchStatus("Reco_QQ_mupl_idx", 1);
  mytree->SetBranchStatus("Reco_QQ_mumi_idx", 1);
  mytree->SetBranchStatus("Reco_QQ_whichGen", 1);
  mytree->SetBranchStatus("Reco_mu_dxy", 1);
  mytree->SetBranchStatus("Reco_mu_dz", 1);
  mytree->SetBranchStatus("Reco_mu_nTrkWMea", 1);
  mytree->SetBranchStatus("Reco_mu_nPixWMea", 1);
  mytree->SetBranchStatus("Reco_QQ_sign", 1);
  mytree->SetBranchStatus("Reco_mu_SelectionType", 1);
  mytree->SetBranchStatus("Reco_mu_whichGen", 1);
  mytree->SetBranchStatus("Gen_weight", 1);
  mytree->AddBranchToCache("*", kTRUE);
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id
  Short_t Gen_QQ_mupl_idx[maxBranchSize];
  Short_t Gen_QQ_mumi_idx[maxBranchSize];
  TBranch *b_Gen_QQ_mupl_idx;
  TBranch *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);

  Short_t Gen_mu_charge[maxBranchSize];
  TBranch *b_Gen_mu_charge; //!
  mytree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);

  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;
  TClonesArray *Reco_QQ_4mom;
  TClonesArray *Reco_mu_4mom;
  ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
  ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_QQ_size]
  Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
  Float_t Reco_QQ_ctau3D[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_QQ_size;                //!
  TBranch *b_Reco_mu_size;                //!
  TBranch *b_Reco_QQ_4mom;                //!
  TBranch *b_Reco_mu_4mom;                //!
  TBranch *b_Reco_QQ_trig;                //!
  TBranch *b_Reco_mu_trig;                //!
  TBranch *b_Reco_QQ_VtxProb;             //!
  TBranch *b_Reco_QQ_ctau3D;              //!

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);

  //  muon id
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_whichGen[maxBranchSize];   //[Reco_QQ_size]
  TBranch *b_Reco_QQ_mupl_idx;
  TBranch *b_Reco_QQ_mumi_idx;
  TBranch *b_Reco_QQ_whichGen;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
  mytree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);

  Float_t Reco_mu_dxy[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dxy;             //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  Float_t Reco_mu_dz[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dz;             //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nTrkWMea;           //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nPixWMea;           //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Short_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_QQ_sign;             //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t Reco_mu_SelectionType[maxBranchSize];
  TBranch *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Short_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);

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
  const ULong64_t kTrigSelMask = (1ULL << kTrigSel);
  const ULong64_t kL2FilterMask = (1ULL << kL2filter);
  const ULong64_t kL3FilterMask = (1ULL << kL3filter);
  const int kMuonSelGlobal = (1 << 1);
  const int kMuonSelTMOneStation = (1 << 3);

  int count = 0;
  int counttnp = 0;
  Long64_t totalEntries = mytree->GetEntries();
  Long64_t nevt = totalEntries;
  if(maxEvents > 0 && nevt > maxEvents) nevt = maxEvents;
  // const int nevt = mytree->GetEntries();
  //nevt = 1000;
  //nevt = 5000000;
  cout << "Total Events : " << nevt;
  if(nevt < totalEntries) cout << " (limited from " << totalEntries << ")";
  cout << endl;
  //for (int iev = 0; iev < nevt; ++iev)
   for(Long64_t iev=0; iev<nevt; ++iev)
  {
    if (reportEvery > 0 && iev % reportEvery == 0)
      cout << ">>>>> EVENT " << iev << " / " << nevt << " (" << (int)(100. * iev / nevt) << "%)" << endl;

    mytree->GetEntry(iev);
    weight = findNcoll(Centrality) * Gen_weight;
    /// Gen_QQ_size = 1;
    // if(Gen_QQ_size > 0) cout<<"weight : "<<weight<<", Gen_QQ_size : "<<Gen_QQ_size<<", Reco_QQ_size : "<<Reco_QQ_size<<endl;

    for (int igen = 0; igen < Gen_QQ_size; igen++)
    {
      JP_Gen = (TLorentzVector *)Gen_QQ_4mom->At(igen);
      mupl_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      mumi_Gen = (TLorentzVector *)Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      Double_t Rapidity_g = fabs(JP_Gen->Rapidity());
      // cout<<"Rapidity_g : "<<Rapidity_g<<", pt : "<<JP_Gen->Pt()<<endl;
      if (!(JP_Gen->M() > massLow && JP_Gen->M() < massHigh))
        continue;

      if (!(JP_Gen->Pt() > 3 && JP_Gen->Pt() < 50 && fabs(JP_Gen->Rapidity()) < 2.4 && IsAcceptanceQQ(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptanceQQ(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
        continue;

      if (!(fabs(JP_Gen->Rapidity()) < 2.4 && fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
        continue;
      if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
        continue;

      pt_weight = EvalPtWeight(applyPtWeight, JP_Gen->Rapidity(), JP_Gen->Pt(), fptwsysMid, fptwsysFor);
      if (JP_Gen->Pt() > 6.5)
        hy_gen->Fill(Rapidity_g, weight * pt_weight);
      if ((Centrality >= cLow && Centrality < cHigh))
      {
        if (JP_Gen->Pt() > 3. && JP_Gen->Pt() < 6.5 && Rapidity_g > 1.6 && Rapidity_g < 2.4) { hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight); hpt_gen_0_noW->Fill(JP_Gen->Pt(), 1.0); }
        if (JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 50 && Rapidity_g < 2.4) { hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight); hpt_gen_0_noW->Fill(JP_Gen->Pt(), 1.0);}
        if (Rapidity_g > 1.6 && Rapidity_g < 2.4 && JP_Gen->Pt() > 3.5 && JP_Gen->Pt() < 40) { hpt_gen_1->Fill(JP_Gen->Pt(), weight * pt_weight); hpt_gen_1_noW->Fill(JP_Gen->Pt(), 1.0);}
        if (Rapidity_g < 1.6 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 40) { hpt_gen_2->Fill(JP_Gen->Pt(), weight * pt_weight); hpt_gen_2_noW->Fill(JP_Gen->Pt(), 1.0); }
      }
      if (Rapidity_g > 1.6 && Rapidity_g < 2.4 && JP_Gen->Pt() > 3.5 && JP_Gen->Pt() < 40)
      {
        hcent_gen_1->Fill(Centrality, weight * pt_weight);
        hInt_gen_1->Fill(1, weight * pt_weight);
        hInt_gen_1_noW->Fill(1, 1.0);
      }
      if (Rapidity_g <= 1.6 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 40)
      {
        hcent_gen_2->Fill(Centrality, weight * pt_weight);
        hInt_gen_2->Fill(1, weight * pt_weight);
        hInt_gen_2_noW->Fill(1, 1.0);
      }
      if (Rapidity_g < 2.4 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 50)
      {
        hcent_gen_0->Fill(Centrality, weight * pt_weight);
      }
    }

    bool HLTPass = false;
    if ((HLTriggers & kTrigSelMask) == kTrigSelMask)
      HLTPass = true;
    // if (HLTPass == false)
    //   continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; irqq++) 
    {
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      bool HLTFilterPass=false;
      if ((Reco_QQ_trig[irqq] & kTrigSelMask) == kTrigSelMask) HLTFilterPass=true;
      if(HLTFilterPass==false) continue;

      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      if(Reco_QQ_sign[irqq]!=0) continue;  

      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;
      
      if(abs(JP_Reco->Rapidity())>yHigh || abs(JP_Reco->Rapidity())<yLow) continue;
      if(JP_Reco->Pt()<ptLow || JP_Reco->Pt()>ptHigh) continue;
      if(JP_Reco->M()<massLow || JP_Reco->M()>massHigh) continue;
      pt1 = mupl_Reco->Pt();
      pt2 = mumi_Reco->Pt();
      eta1 = mupl_Reco->Eta();
      eta2 = mumi_Reco->Eta();

      if(!(IsAcceptanceQQ(pt1,eta1) && IsAcceptanceQQ(pt2,eta2))) continue;
      //if(!settree_.SoftMuIdCut(irqq)) continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & kMuonSelGlobal);
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & kMuonSelTMOneStation);

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & kMuonSelGlobal);
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & kMuonSelTMOneStation);
      if(passMuonTypePl==false || passMuonTypeMi==false) continue;

      phi1 = mupl_Reco->Phi();
      phi2 = mumi_Reco->Phi();

      TLorentzVector mudiff_Reco(*mupl_Reco - *mumi_Reco);
      dcosMuTheta = (mupl_Reco->P()*mupl_Reco->P() + mumi_Reco->P()*mumi_Reco->P() - mudiff_Reco.P()*mudiff_Reco.P())/(2*mupl_Reco->P()*mumi_Reco->P());
      ptimb = abs(pt1-pt2)/(pt1+pt2);

      bool muplSoft = (  //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl        //                       &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //                        &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      if ( !(muplSoft && mumiSoft) ) 
        continue;   

	  if(isTnP){
		  tnp_weight = 1;
		  tnp_trig_weight_mupl = -1;
		  tnp_trig_weight_mumi = -1;
		  tnp_weight = tnp_weight * tnp_weight_muid_pbpb(pt1, eta1, 0) * tnp_weight_muid_pbpb(pt2, eta2, 0); //mu id
		  tnp_weight = tnp_weight * tnp_weight_trk_pbpb(eta1, 0) * tnp_weight_trk_pbpb(eta2, 0); //inner tracker

		  if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & kL2FilterMask) == kL2FilterMask && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & kL2FilterMask) == kL2FilterMask ) ) continue;

		  bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & kL2FilterMask) == kL2FilterMask) ? true : false;
		  bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & kL3FilterMask) == kL3FilterMask) ? true : false;
		  bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & kL2FilterMask) == kL2FilterMask) ? true : false;
		  bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & kL3FilterMask) == kL3FilterMask) ? true : false;
		  if(mupl_L2Filter == false || mumi_L2Filter == false){ cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl; cout << endl;}

		  bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
		  bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
		  bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
		  bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
		  bool SelDone = false;

		  if( mupl_isL2 && mumi_isL3){
			  tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0); 
			  tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0); 
			  SelDone = true;
		  }   
		  else if( mupl_isL3 && mumi_isL2){
			  tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0); 
			  tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0); 
			  SelDone = true;
		  }   
		  else if( mupl_isL3 && mumi_isL3){
			  int t[2] = {-1,1}; // mupl, mumi
			  int l = rand() % (2);
			  //pick up what will be L2
			  if(t[l]==-1){
				  tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
				  tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
			  }
			  else if(t[l]==1){
				  tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
				  tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
			  }
			  SelDone = true;
		  }
		  if(SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
		  tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
	  }

      pt_weight = EvalPtWeight(applyPtWeight, JP_Reco->Rapidity(), JP_Reco->Pt(), fptwsysMid, fptwsysFor);

      // cout<<"pt_Weight in reco : "<<pt_weight<<endl;
      Double_t Rapidity = fabs(JP_Reco->Rapidity());
      if (HLTPass == true && HLTFilterPass == true)
      {
        if (JP_Reco->Pt() > 6.5)
          hy_reco->Fill(Rapidity, weight * tnp_weight * pt_weight);

        bool passCtauPtBinned = false;
        if ((Centrality >= cLow && Centrality < cHigh))
        {
          if (Rapidity < 2.4 && Rapidity > 1.6 && JP_Reco->Pt() > 3. && JP_Reco->Pt() < 6.5)
          {
            hpt_reco_0->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
            hpt_reco_0_noW->Fill(JP_Reco->Pt(), 1.0);
          }
          if (Rapidity < 2.4 && JP_Reco->Pt() > 6.5 && JP_Reco->Pt() < 50)
          {
            hpt_reco_0->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
            hpt_reco_0_noW->Fill(JP_Reco->Pt(), 1.0);
          }
          if (Rapidity > 1.6 && Rapidity < 2.4)
          {
            double pt = JP_Reco->Pt();
            // Start from i=1 to skip the first bin [0-3.5 GeV]
            for(int i=1; i < 5 ; i++) {
              double ptLow = ptBin_for[i];
              double ptHigh = ptBin_for[i+1];
              if(pt >= ptLow && pt < ptHigh) {
                double l_cut = l_cut_for_pt[i];
                if(state==1) { 
                  if(Reco_QQ_ctau3D[irqq] < l_cut) {
                    passCtauPtBinned = true;
                    hpt_reco_1->Fill(pt, weight * tnp_weight * pt_weight);
                    hpt_reco_1_noW->Fill(pt, 1.0);
                    //if(count < 10) cout << "Forward: pT=" << pt << " (bin " << ptLow << "-" << ptHigh << "), l_cut[" << i << "]=" << l_cut << ", ctau3D=" << Reco_QQ_ctau3D[irqq] << " -> FILLED" << endl;
                    cout.precision(4);
                    //cout << "pT = " << pt << "\ty = " << Rapidity << "\tweight = " << weight << "\tpT weight = " << pt_weight << "\tTnP weight = " << tnp_weight << endl;
                  }
                }
                else if(state==2) { 
                  if(Reco_QQ_ctau3D[irqq] > l_cut) {
                    passCtauPtBinned = true;
                    hpt_reco_1->Fill(pt, weight * tnp_weight * pt_weight);
                    hpt_reco_1_noW->Fill(pt, 1.0);
                    //if(count < 10) cout << "Forward: pT=" << pt << " (bin " << ptLow << "-" << ptHigh << "), l_cut[" << i << "]=" << l_cut << ", ctau3D=" << Reco_QQ_ctau3D[irqq] << " -> FILLED" << endl;
                  }
                }
                break;
              }
            }
          }
          if (Rapidity < 1.6)
          {
            double pt = JP_Reco->Pt();
            // Start from i=1 to skip the first bin [0-6.5 GeV]
            for(int i=1; i < 7 ; i++) {
              double ptLow = ptBin_mid[i];
              double ptHigh = ptBin_mid[i+1];
              if(pt >= ptLow && pt < ptHigh) {
                double l_cut = l_cut_mid_pt[i];
                if(state==1) { 
                  if(Reco_QQ_ctau3D[irqq] < l_cut) { passCtauPtBinned = true; hpt_reco_2->Fill(pt, weight * tnp_weight * pt_weight); hpt_reco_2_noW->Fill(pt, 1.0); }
                  //cout << "pT : " << pt << ",\tl_cut : " << l_cut << ",\tctau3D : " << Reco_QQ_ctau3D[irqq] << endl;
                }
                else if(state==2) { 
                  if(Reco_QQ_ctau3D[irqq] > l_cut) { passCtauPtBinned = true; hpt_reco_2->Fill(pt, weight * tnp_weight * pt_weight); hpt_reco_2_noW->Fill(pt, 1.0); } 
                }
                break;
              }
            }
          }

          if (passCtauPtBinned && Rapidity > 1.6 && Rapidity < 2.4)
          {
            if (JP_Reco->Pt() > 3.5 && JP_Reco->Pt() < 40)
            {
              hcent_reco_1->Fill(Centrality, weight * tnp_weight * pt_weight);
              if(Centrality >= 0 && Centrality < 180)
              {
                hInt_reco_1->Fill(1, weight * tnp_weight * pt_weight);
                hInt_reco_1_noW->Fill(1, 1.0);
              }
            }
          }

          if (passCtauPtBinned && Rapidity < 1.6 && JP_Reco->Pt() > 6.5 && JP_Reco->Pt() < 40)
          {
            hcent_reco_2->Fill(Centrality, weight * tnp_weight * pt_weight);
            if(Centrality >= 0 && Centrality < 180)
            {
              hInt_reco_2->Fill(1, weight * tnp_weight * pt_weight);
              hInt_reco_2_noW->Fill(1, 1.0);
            }
          }
        }
      }
        if (Rapidity < 2.4 && JP_Reco->Pt() > 6.5 && JP_Reco->Pt() < 50)
        {
          hcent_reco_0->Fill(Centrality, weight * tnp_weight * pt_weight);
        }
      }
    }

  //cout << "count " << count << endl;
  //cout << "counttnp " << counttnp << endl;

  // Divide
  TH1D *hpt_eff_0;
  TH1D *hpt_eff_1;
  TH1D *hpt_eff_2;
  TH1D *hpt_eff_0_noW;
  TH1D *hpt_eff_1_noW;
  TH1D *hpt_eff_2_noW;

  TH1D *hcent_eff_0;
  TH1D *hcent_eff_1;
  TH1D *hcent_eff_2;

  TH1D *hInt_eff_1;
  TH1D *hInt_eff_2;
  TH1D *hInt_eff_1_noW;
  TH1D *hInt_eff_2_noW;

  hpt_eff_0 = (TH1D *)hpt_reco_0->Clone("hpt_eff_0");
  hpt_eff_1 = (TH1D *)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D *)hpt_reco_2->Clone("hpt_eff_2");
  hpt_eff_0_noW = (TH1D *)hpt_reco_0_noW->Clone("hpt_eff_0_noW");
  hpt_eff_1_noW = (TH1D *)hpt_reco_1_noW->Clone("hpt_eff_1_noW");
  hpt_eff_2_noW = (TH1D *)hpt_reco_2_noW->Clone("hpt_eff_2_noW");

  hpt_eff_0->Divide(hpt_eff_0, hpt_gen_0, 1, 1);
  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1);
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1);
  hpt_eff_0_noW->Divide(hpt_eff_0_noW, hpt_gen_0_noW, 1, 1);
  hpt_eff_1_noW->Divide(hpt_eff_1_noW, hpt_gen_1_noW, 1, 1);
  hpt_eff_2_noW->Divide(hpt_eff_2_noW, hpt_gen_2_noW, 1, 1);

  hcent_eff_0 = (TH1D *)hcent_reco_0->Clone("hcent_eff_0");
  hcent_eff_1 = (TH1D *)hcent_reco_1->Clone("hcent_eff_1");
  hcent_eff_2 = (TH1D *)hcent_reco_2->Clone("hcent_eff_2");

  hcent_eff_0->Divide(hcent_eff_0, hcent_gen_0, 1, 1);
  hcent_eff_1->Divide(hcent_eff_1, hcent_gen_1, 1, 1);
  hcent_eff_2->Divide(hcent_eff_2, hcent_gen_2, 1, 1);

  hInt_eff_1 = (TH1D *)hInt_reco_1->Clone("hInt_eff_1");
  hInt_eff_2 = (TH1D *)hInt_reco_2->Clone("hInt_eff_2");
  hInt_eff_1_noW = (TH1D *)hInt_reco_1_noW->Clone("hInt_eff_1_noW");
  hInt_eff_2_noW = (TH1D *)hInt_reco_2_noW->Clone("hInt_eff_2_noW");

  hInt_eff_1->Divide(hInt_eff_1, hInt_gen_1, 1, 1);
  hInt_eff_2->Divide(hInt_eff_2, hInt_gen_2, 1, 1);
  hInt_eff_1_noW->Divide(hInt_eff_1_noW, hInt_gen_1_noW, 1, 1);
  hInt_eff_2_noW->Divide(hInt_eff_2_noW, hInt_gen_2_noW, 1, 1);

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
  hpt_eff_0_noW->SetTitle("Eff (no weight): Rapidity 0.0-2.4");
  hpt_eff_1_noW->SetTitle("Eff (no weight): Rapidity 1.6-2.4");
  hpt_eff_2_noW->SetTitle("Eff (no weight): Rapidity 0.0-1.6");
  hInt_eff_1_noW->SetTitle("Integrated Eff (no weight): Rapidity 1.6-2.4");
  hInt_eff_2_noW->SetTitle("Integrated Eff (no weight): Rapidity 0.0-1.6");
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
  hpt_eff_2->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6", isTnP, 1, cLow, cHigh));
  hcent_eff_0->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_6p5_to_50_absy0_2p4", isTnP, 1));
  hcent_eff_1->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_3_to_40_absy1p6_2p4", isTnP, 1));
  hcent_eff_2->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_6p5_to_40_absy0_1p6", isTnP, 1));
  hy_eff->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d", isTnP, 1));
  hInt_eff_1->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4", isTnP, 1, cLow, cHigh));
  hInt_eff_2->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6", isTnP, 1, cLow, cHigh));
  hpt_eff_0_noW->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_2p4_noW", isTnP, 1, cLow, cHigh));
  hpt_eff_1_noW->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4_noW", isTnP, 1, cLow, cHigh));
  hpt_eff_2_noW->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6_noW", isTnP, 1, cLow, cHigh));
  hInt_eff_1_noW->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4_noW", isTnP, 1, cLow, cHigh));
  hInt_eff_2_noW->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6_noW", isTnP, 1, cLow, cHigh));

  // TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_drawsame1.root",cLow,cHigh,isPtWeight,isTnP);
  //TString outFileName = Form("./roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260105.root", cLow, cHigh, ptSys.Data(), isTnP);
  TString outFileName;
  gSystem->mkdir(Form("%s/roots", baseDir.Data()), true);
  if (state == 1 && applyPtWeight) outFileName = Form("%s/roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260422_v7.root", baseDir.Data(), cLow, cHigh, ptSys.Data(), isTnP);
  else if (state == 1 && !applyPtWeight) outFileName = Form("%s/roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_JPsi_PtW%d_tnp%d_ctauCut_260422_v7_noPtW.root", baseDir.Data(), cLow, cHigh, 1, isTnP);
  else if (state == 2 && applyPtWeight) outFileName = Form("%s/roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_nprompt_pbpb_JPsi_PtW%s_tnp%d_ctauCut_260422_v7.root", baseDir.Data(), cLow, cHigh, ptSys.Data(), isTnP);
  else if (state == 2 && !applyPtWeight) outFileName = Form("%s/roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_nprompt_pbpb_JPsi_PtW%d_tnp%d_ctauCut_260422_v7_noPtW.root", baseDir.Data(), cLow, cHigh, 1, isTnP);
  TFile *outFile = new TFile(outFileName, "RECREATE");
  hpt_eff_0->Write();
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hcent_eff_0->Write();
  hcent_eff_1->Write();
  hcent_eff_2->Write();
  hInt_eff_1->Write();
  hInt_eff_2->Write();
  hpt_eff_0_noW->Write();
  hpt_eff_1_noW->Write();
  hpt_eff_2_noW->Write();
  hInt_eff_1_noW->Write();
  hInt_eff_2_noW->Write();
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

  hpt_reco_0_noW->Write();
  hpt_reco_1_noW->Write();
  hpt_reco_2_noW->Write();
  hpt_gen_0_noW->Write();
  hpt_gen_1_noW->Write();
  hpt_gen_2_noW->Write();
  hInt_reco_1_noW->Write();
  hInt_reco_2_noW->Write();
  hInt_gen_1_noW->Write();
  hInt_gen_2_noW->Write();

  hy_gen->Write();
  hy_reco->Write();

  outFile->Close();

  cout << "Save file : " << outFileName.Data() << endl;
}
