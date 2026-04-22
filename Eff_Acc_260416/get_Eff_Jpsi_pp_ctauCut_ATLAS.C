#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
//#include "Style_jaebeom.h"
//#include "tnp_weight_lowptPbPb.h"
#include "../tnp_weight_pp.h"
#include <TAttMarker.h>

using namespace std;
// v4 : change systematic up and down

void get_Eff_Jpsi_pp_ctauCut_ATLAS(
  int state=1,
  bool isTnP = true, int isPtWeight = 0, // 0 : nominal , 1 : up, -1 : down
  float ptLow = 6.5, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4
  //float cLow = 0, float cHigh = 20, 
  //float cLow = 0, float cHigh = 200, 
  //bool isTnP = false, bool isPtWeight = false, int state=1
  //bool isTnP = true, bool isPtWeight = false, int state=1
  //bool isTnP = false, bool isPtWeight = true, int state=1
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 3;	//jpsi=12,Upsilon=13 pp Jpsi=3

  float muPtCut = 0; //3.5, 1.8

  //jpsi
  float massLow = 0.0;
  float massHigh = 10.0;
  //float massLow = 2.6;
  //float massHigh = 3.5;

  double min = 0;
  double max = ptHigh;
  const int numBins = 16; //50;//(max-min)/binwidth;  //31//

  TString ptSys;
  if(isPtWeight==0) ptSys = "nomi";
  else if(isPtWeight==1) ptSys = "up";
  else if(isPtWeight==-1) ptSys = "down";

  //input files
  //PbPb
  //TString inputMC1 = "/work2/Oniatree/JPsi/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  //TString inputMC1 = "/work2/Oniatree/JPsi/OniatreeMC_BJPsiMM_TuneCUETP8M1_5p02TeV_pythia8_230214.root";	//pp_non prompt
  TString inputMC1;
  if(state==1) inputMC1 = "/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root"; // pp prompt
  else if (state==2) inputMC1 = "/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root";
  TChain* mytree;
  mytree = new TChain("hionia/myTree"); 
  //TChain* mytree = new TChain("hionia/myTree"); 
  mytree->Add(inputMC1.Data());
  //TFile *inf = new TFile("../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root","READ");
  //TTree *mytree = (TTree*)inf->Get("myTree");
  //if(state==1) mytree->Add(inputMC2.Data());
  cout<<"dmoon chk entries : "<<mytree->GetEntries()<<endl;
  //TString outFileName = "mc_eff_vs_pt_cent_prompt_pbpb_Jpsi.root"; 
  //if(state==2) outFileName = "mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi.root";

  // pT reweighting function
  // ratioDataMC_AA_Jpsi_DATA_y0.0-2.4_210915.root
  // TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y0_1p6_230321.root","read");
  // TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y1p6_2p4_230420.root","read");
  // TFile *fPtW1 = new TFile("../../compareDataToMC/ratioDataMC_pp_BtoPsi2S_DATA_y0_1p6_230629.root", "read");
  // TFile *fPtW2 = new TFile("../../compareDataToMC/ratioDataMC_pp_BtoPsi2S_DATA_y1p6_2p4_230629.root", "read");
  //jpsi
  TFile *fPtW1;
  TFile *fPtW2;
  if(state==1){
    fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
    fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
  }
  else if(state==2){
    fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
    fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
  }
  //TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
  //TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
  //TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y1p6_2p4_241016.root", "read");
  //TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y0_1p6_241008.root", "read");
  //TFile *fPtW1 = new TFile("./ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
  //TFile *fPtW2 = new TFile("./ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root","read");
  //TFile *fPtW1 = new TFile("./ratioDataMC_pp_BtoJPsi_DATA_y0_2p4_240402_new.root", "read");
  //TFile *fPtW1 = new TFile("../roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
  //psi2s
  //TFile *fPtW1 = new TFile("../../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y0_1p6_230321.root", "read");
  //TFile *fPtW2 = new TFile("../../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y1p6_2p4_230420.root", "read");
  // if(state==2) TFile *fPtW = new TFile("ratioDataMC_AA_btojpsi_DATA_1s.root","read");
  TF1 *fptw1 = (TF1 *)fPtW1->Get("dataMC_Ratio1");
  TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");


  double ptBin_all[18] = {0, 3., 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 22.5, 25, 27.5, 30, 50};
  double ptBin_for[6] = {0,3.5,6.5,9,12,40};
  double ptBin_mid[8] = {0,6.5,9,12,15,20,25,40};
  double centBin_for[4] = {0,40,80,180};
  double centBin_mid[7] = {0,10,20,30,40,50,90};
  float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH1D* hpt_reco_0 = new TH1D("hpt_reco_0","hpt_reco_0",17,ptBin_all);
  TH1D* hpt_reco_1 = new TH1D("hpt_reco_1","hpt_reco_1",5,ptBin_for);
  TH1D* hpt_reco_2 = new TH1D("hpt_reco_2","hpt_reco_2",7,ptBin_mid);

  TH1D* hpt_gen_0 = new TH1D("hpt_gen_0","hpt_gen_0",17,ptBin_all);
  TH1D* hpt_gen_1 = new TH1D("hpt_gen_1","hpt_gen_1",5,ptBin_for);
  TH1D* hpt_gen_2 = new TH1D("hpt_gen_2","hpt_gen_2",7,ptBin_mid);

  TH1D* hInt_gen_1 = new TH1D("hInt_gen_1", "",1,0,40);
  TH1D* hInt_gen_2 = new TH1D("hInt_gen_2", "",1,0,40);

  TH1D* hInt_reco_1 = new TH1D("hInt_reco_1", "",1,0,40);
  TH1D* hInt_reco_2 = new TH1D("hInt_reco_2", "",1,0,40);

  //TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",3,50);
  //f1->SetParameters(214,22,14);
  //f1->SetParLimits(0,0,500);
  //f1->SetParLimits(1,-1,500);
  //f1->SetParLimits(2,0,500);


  TH1D* hy_gen = new TH1D("hy_gen","hy_gen",6,yBin);
  TH1D* hy_reco = new TH1D("hy_reco","hy_reco",6,yBin);
  hy_gen ->Sumw2();
  hy_reco ->Sumw2();

  hpt_reco_0->Sumw2();
  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();

  hpt_gen_0->Sumw2();
  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();

  hInt_gen_1->Sumw2();
  hInt_gen_2->Sumw2();

  hInt_reco_1->Sumw2();
  hInt_reco_2->Sumw2();


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

  TF1 *fptw1sys = new TF1("fptw1sys","( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",6.5, 50.0);
  TF1 *fptw2sys = new TF1("fptw2sys","( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",3.0, 50.0);

  cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

  // sys for nominal
  if(isPtWeight == 0) {
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
  if(isPtWeight == 1) {
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
  if(isPtWeight == -1) {
    p10 = p10 - e10;
    p11 = p11 - e11;
    p12 = p12 - e12;
    p13 = p13 - e13;

    p20 = p20 - e10;
    p21 = p21 - e11;
    p22 = p22 - e12;
    p23 = p23 - e13;
  }

  //cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

  fptw1sys->SetParameter(0, p10);
  fptw1sys->SetParameter(1, p11);
  fptw1sys->SetParameter(2, p12);
  fptw1sys->SetParameter(3, p13);

  fptw2sys->SetParameter(0, p20);
  fptw2sys->SetParameter(1, p21);
  fptw2sys->SetParameter(2, p22);
  fptw2sys->SetParameter(3, p23);

  const int maxBranchSize = 1000;
  ULong64_t       HLTriggers;
  Short_t           Gen_QQ_size;
  Short_t           Gen_mu_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Gen_QQ_trig[maxBranchSize];   //[Gen_QQ_size]
  Float_t         Gen_QQ_VtxProb[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_QQ_trig;   //!
  TBranch        *b_Gen_QQ_VtxProb;   //!

  Gen_QQ_4mom = 0; Gen_mu_4mom = 0;
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id 
  Short_t           Gen_QQ_mupl_idx[maxBranchSize];
  Short_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);

  Short_t           Gen_mu_charge[maxBranchSize];
  TBranch        *b_Gen_mu_charge;   //!
  mytree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);


  Short_t           Reco_QQ_size;
  Short_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Reco_QQ_4mom = 0; Reco_mu_4mom = 0;
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Short_t           Reco_QQ_mupl_idx[maxBranchSize];
  Short_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);


  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Short_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Short_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);

  
  TLorentzVector* JP_Gen= new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;

  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  double weight = 1;
  double tnp_weight = 1;
  double tnp_trig_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double tnp_trig_weight_muplL2_num = -1;
  double tnp_trig_weight_muplL3_num = -1;
  double tnp_trig_weight_mumiL2_num = -1;
  double tnp_trig_weight_mumiL3_num = -1;
  double tnp_trig_weight_muplL2_den = -1;
  double tnp_trig_weight_muplL3_den = -1;
  double tnp_trig_weight_mumiL2_den = -1;
  double tnp_trig_weight_mumiL3_den = -1;
  double tnp_trig_weight_num = 1;
  double tnp_trig_weight_den = 1;
  double pt_weight = 1;
  
  double tnp_trig_dimu=-1;
//  TH2D* hpt_tnp_trig = new TH2D("hpt_tnp_trig","hpt_tnp_trig",numBins,ptBin,40,0,2);

  int kL2filter = 4;	//jpsi=16,Upsilon=38
  int kL3filter = 5;	//jpsi=17,Upsilon=39

  int count =0;
  int counttnp =0;
  int nevt = mytree->GetEntries();
  int nevt_2 = nevt / 100;
  //const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  double counts[10]={0};
  //nevt = 100;
  for (int iev = 0; iev < nevt; ++iev)
  //for(int iev=0; iev<300000 ; ++iev)
  {
    if (iev % 100000 == 0)
      //cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() << " (" << (int)(100. * iev / mytree->GetEntries()) << "%)" << endl;
      cout << ">>>>> EVENT " << iev << " / " << nevt << " (" << (int)(100. * iev / nevt) << "%)" << endl;

    mytree->GetEntry(iev);
    //weight = findNcoll(Centrality) * Gen_weight;
    weight = 1.;
    //weight = Gen_weight;
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

      // if(! (JP_Gen->Pt() > 3 && JP_Gen->Pt()<50 && fabs(JP_Gen->Rapidity())<2.4 && IsAcceptanceQQ(mupl_Gen->Pt(),fabs(mupl_Gen->Eta()))&&IsAcceptanceQQ(mumi_Gen->Pt(),fabs(mumi_Gen->Eta()))) ) continue;
      if (!(JP_Gen->Pt() < 50 && fabs(JP_Gen->Rapidity()) < 2.4 && IsAcceptanceQQ(mupl_Gen->Pt(), fabs(mupl_Gen->Eta())) && IsAcceptanceQQ(mumi_Gen->Pt(), fabs(mumi_Gen->Eta()))))
        continue;

      if (!(fabs(JP_Gen->Rapidity()) < 2.4 && fabs(mupl_Gen->Eta()) < 2.4 && fabs(mumi_Gen->Eta()) < 2.4))
        continue;
      if (Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0)
        continue;

      pt_weight = 1;
      if (isPtWeight && fabs(JP_Gen->Rapidity()) > 1.6 && fabs(JP_Gen->Rapidity() < 2.4))
        pt_weight = fptw1sys->Eval(JP_Gen->Pt());
      if (isPtWeight && fabs(JP_Gen->Rapidity()) < 1.6)
        pt_weight = fptw2sys->Eval(JP_Gen->Pt());

	  //if( (double)(pt_weight) < 0 && JP_Gen->Pt()>3.0) cout << ">>>>> Event " << iev << " / " << nevt << " (" << (double)(pt_weight) << ")" << "\tpt : "<< JP_Gen->Pt() << "\ty : " << JP_Gen->Rapidity() << endl;

      if (JP_Gen->Pt() > 6.5)
        hy_gen->Fill(Rapidity_g, weight * pt_weight);
      if (Rapidity_g < 2.4)
      {
        hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight);
        //hpt_gen_0->Fill(JP_Gen->Pt(), weight * pt_weight);
      }
      if (Rapidity_g > 1.6 && Rapidity_g < 2.4)
      {
        hpt_gen_1->Fill(JP_Gen->Pt(), weight * pt_weight);
        hInt_gen_1->Fill(1, weight * pt_weight);
      }
      if (Rapidity_g < 2 && JP_Gen->Pt() > 9 && JP_Gen->Pt() < 40)
      {
        hpt_gen_2->Fill(JP_Gen->Pt(), weight * pt_weight);
        hInt_gen_2->Fill(1, weight * pt_weight);
      }

      // if(! (Centrality > cLow && Centrality < cHigh)) continue;
      // cout<<"pt_Weight in gen : "<<pt_weight<<endl;

      // if(pt_weight < 0) cout<<"Rapidity_g : "<<Rapidity_g<<", pt : "<<JP_Gen->Pt()<<", weight : "<<weight<<", pt_weight : "<<pt_weight<<endl;
    }

    bool HLTPass = false;
    if ((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
      HLTPass = true;

    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      // cout<<"Reco QQ pt : "<<JP_Reco->Pt()<<endl;
      Double_t Rapidity_r = fabs(JP_Reco->Rapidity());

      bool HLTFilterPass = false;
	  counts[0]++;
      if ((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))
        HLTFilterPass = true;

      if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
        continue;
	  counts[1]++;
      if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
        continue;
	  counts[2]++;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));

      bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) &&
          passMuonTypePl //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
      );

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) &&
          passMuonTypeMi //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
      );

      if (!(muplSoft && mumiSoft))
        continue;
      if (Reco_QQ_VtxProb[irqq] < 0.01)
        continue;
      if (Reco_QQ_sign[irqq] != 0)
        continue;
	  counts[2]++;

      Double_t Rapidity = fabs(JP_Reco->Rapidity());

      // if(! (JP_Reco->Pt()>3&&JP_Reco->Pt()<50&&fabs(JP_Reco->Rapidity())<2.4&&IsAcceptanceQQ(mupl_Reco->Pt(),fabs(mupl_Reco->Eta()))&&IsAcceptanceQQ(mumi_Reco->Pt(),fabs(mumi_Reco->Eta()))) ) continue;
      if (!(JP_Reco->Pt() < 50 && fabs(JP_Reco->Rapidity()) < 2.4 && IsAcceptanceQQ(mupl_Reco->Pt(), fabs(mupl_Reco->Eta())) && IsAcceptanceQQ(mumi_Reco->Pt(), fabs(mumi_Reco->Eta()))))
        continue;
	  counts[2]++;

      if (!(fabs(mupl_Reco->Eta()) < 2.4 && fabs(mumi_Reco->Eta()) < 2.4 && fabs(JP_Reco->Rapidity()) < 2.4 && JP_Reco->M() > massLow && JP_Reco->M() < massHigh))
        continue;
	  counts[3]++;

      if (HLTPass == true && HLTFilterPass == true)
        count++;
      if (isTnP)
      {
        tnp_weight = 1;
        tnp_trig_weight = 1;
        tnp_trig_weight_mupl = -1;
        tnp_trig_weight_mumi = -1;
        tnp_trig_weight_muplL2_num = -1;
        tnp_trig_weight_muplL3_num = -1;
        tnp_trig_weight_mumiL2_num = -1;
        tnp_trig_weight_mumiL3_num = -1;
        tnp_trig_weight_muplL2_den = -1;
        tnp_trig_weight_muplL3_den = -1;
        tnp_trig_weight_mumiL2_den = -1;
        tnp_trig_weight_mumiL3_den = -1;
        tnp_trig_weight_num = 1;
        tnp_trig_weight_den = 1;

        tnp_weight = tnp_weight * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mupl_Reco->Pt(), mupl_Reco->Eta())) * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mumi_Reco->Pt(), mumi_Reco->Eta())); // mu id

        tnp_weight = tnp_weight * tnp_trig_weight;

        // if(SelDone == false){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
        // if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
        if (HLTPass == true && HLTFilterPass == true)
        {
          counttnp++;
          tnp_trig_dimu = tnp_trig_weight;
          // hpt_tnp_trig->Fill(JP_Reco->Pt(),tnp_trig_dimu);
        }
      }
      //weight = 1;
      pt_weight = 1;
      if (isPtWeight && fabs(JP_Reco->Rapidity()) > 1.6 && fabs(JP_Reco->Rapidity()) < 2.4)
        pt_weight = fptw1sys->Eval(JP_Reco->Pt());
      if (isPtWeight && fabs(JP_Reco->Rapidity()) < 1.6)
        pt_weight = fptw2sys->Eval(JP_Reco->Pt());

      // cout<<"pt_Weight in reco : "<<pt_weight<<endl;
      if (HLTPass == true && HLTFilterPass == true)
      {
        hy_reco->Fill(Rapidity_r, weight * tnp_weight * pt_weight);
        if (Rapidity < 2.4)
        {
          hpt_reco_0->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
	  counts[4]++;
          //hpt_reco_0->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
        }
        if (Rapidity > 1.6 && Rapidity < 2.4)
        {
          hpt_reco_1->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
          hInt_reco_1->Fill(1, weight * tnp_weight * pt_weight);
        }
        if (Rapidity < 2 && JP_Reco->Pt() > 9 && JP_Reco->Pt() < 40)
        {
          hpt_reco_2->Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
          hInt_reco_2->Fill(1, weight * tnp_weight * pt_weight);
        }

        // if(! (Centrality > cLow && Centrality < cHigh)) continue;
      }
    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  for (int i=0; i<5; i++) {
	  cout << "counts[" << i << "]" << " : " << counts[i] << endl;
  }

  // Divide
  TH1D *hpt_eff_0;
  TH1D *hpt_eff_1;
  TH1D* hpt_eff_2;


  hpt_eff_0 = (TH1D*)hpt_reco_0->Clone("hpt_eff_0");
  hpt_eff_1 = (TH1D*)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D*)hpt_reco_2->Clone("hpt_eff_2");


  hpt_eff_0->Divide(hpt_eff_0, hpt_gen_0, 1, 1, "B");
  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1, "B");
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1, "B");

  TH1D* hInt_eff_1;
  TH1D* hInt_eff_2;


  hInt_eff_1 = (TH1D*)hInt_reco_1->Clone("hInt_eff_1");
  hInt_eff_2 = (TH1D*)hInt_reco_2->Clone("hInt_eff_2");


  hInt_eff_1->Divide(hInt_eff_1, hInt_gen_1, 1, 1, "B");
  hInt_eff_2->Divide(hInt_eff_2, hInt_gen_2, 1, 1, "B");

  TH1D* hy_eff;
  hy_eff = (TH1D*)hy_reco->Clone("hy_eff");
  hy_eff->Divide(hy_eff, hy_gen, 1, 1, "B");

  hpt_eff_0->SetTitle("Eff: Rapidity 0.0-2.4");
  hpt_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hpt_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");

  
  //f1->SetLineColor(kBlack);
  //f1->SetLineWidth(2);
  //hpt_eff_1->Fit(f1);
  //f1->SetLineColor(kRed+1);
  //f1->SetLineWidth(2);
  //hpt_eff_2->Fit(f1);
  //f1->SetLineColor(kOrange+1);
  //f1->SetLineWidth(2);
  


//draw same
  //gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend = new TLegend(0.6, 0.84);

  TLatex *lt3 = new TLatex();
  lt3->SetNDC();
  lt3->SetTextSize(0.03);

  gStyle->SetOptFit(0);

  TCanvas * cpt_eff_0 = new TCanvas("cpt_eff_0","cpt_eff_0",0,0,900,800);
  cpt_eff_0->cd();
  hpt_eff_0->SetMarkerStyle(24);
  hpt_eff_0->SetMarkerColor(1);
  hpt_eff_0->SetLineColor(1);
  hpt_eff_0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hpt_eff_0->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_0->GetYaxis()->SetRangeUser(0.,0.9);
  hpt_eff_0->Draw("E");
  lt3->SetTextSize(0.03);
  lt3->DrawLatex(0.13,0.84,"Non-Prompt #Jpsi (pp)");
  //cpt_eff_0->SaveAs("./figs/pp/Eff_pt_pp_0to2p4_NP.png");

  TCanvas * cpt_eff = new TCanvas("cpt_eff","cpt_eff",0,0,900,800);
  cpt_eff->cd();
  hpt_eff_1->SetMarkerStyle(24);
  hpt_eff_1->SetMarkerColor(1);
  hpt_eff_1->SetLineColor(1);
  hpt_eff_1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hpt_eff_1->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_1->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_1->Draw("E");
  //eeeeeee
  //hpt_eff_2->SetMarkerStyle(24);
  //hpt_eff_2->SetMarkerColor(1);
  //hpt_eff_2->SetLineColor(1);
  //hpt_eff_2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  //hpt_eff_2->GetYaxis()->SetTitle("Efficiency");
  //hpt_eff_2->GetYaxis()->SetRangeUser(0.,1.3);
  //hpt_eff_2->Draw("E");
  //eeeeeee
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #Jpsi (pp)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #Jpsi (pp)");
  //drawsame
  /*
  hpt_eff_2->SetMarkerStyle(25);
  hpt_eff_2->SetMarkerColor(kRed+1);
  hpt_eff_2->SetLineColor(kRed+1);
  hpt_eff_2->Draw("same");
  legend->AddEntry("hpt_eff_1","|y|: 1.6-2.4","lep");
  legend->AddEntry("hpt_eff_2","|y|: 0-1.6","lep");
  legend->SetBorderSize( 0);
  legend->Draw("E");
  cpt_eff->SaveAs("./figs/pp/Eff_pt_pp_PR.png");

  TCanvas * cy_eff = new TCanvas("cy_eff","cy_eff",0,0,900,800);
  cy_eff->cd();
  hy_eff->SetMarkerStyle(24);
  hy_eff->SetMarkerColor(1);
  hy_eff->SetLineColor(1);
  hy_eff->GetXaxis()->SetTitle("|y|");
  hy_eff->GetYaxis()->SetTitle("Efficiency");
  hy_eff->GetYaxis()->SetRangeUser(0.,1.2);
  hy_eff->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #Jpsi (pp)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #Jpsi (pp)");
  cy_eff->SaveAs("./figs/pp/Eff_absy_noweight_PR.png");*/

  //Save efficiency files for later use.

  hpt_eff_0 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_absy0_2p4",isTnP, 1));
  hpt_eff_1 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_absy1p6_2p4",isTnP, 1));
  hpt_eff_2 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_absy0_2",isTnP, 1));
  hInt_eff_1 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_absy1p6_2p4",isTnP, 1));
  hInt_eff_2 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_absy0_2",isTnP, 1));
  hy_eff ->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d",isTnP, 1));

  // TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_drawsame1.root",cLow,cHigh,isPtWeight,isTnP);
  //TString outFileName = Form("./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW%d_tnp%d_20230416.root", isPtWeight, isTnP);
  //if (state == 2)
  TString  outFileName;
  if(state==1) outFileName = Form("./roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW%s_tnp%d_ctauCut_ATLAS.root", ptSys.Data(), isTnP);
  else if(state==2) outFileName = Form("./roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtW%s_tnp%d_ctauCut_ATLAS.root", ptSys.Data(), isTnP);
  TFile *outFile = new TFile(outFileName, "RECREATE");
  hpt_eff_0->Write();
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hInt_eff_1->Write();
  hInt_eff_2->Write();
  hy_eff->Write();
  cpt_eff->Write();
  cpt_eff_0->Write();
  //  if(isTnP) hpt_tnp_trig->Write();

  hpt_reco_0->Write();
  hpt_reco_1->Write();
  hpt_reco_2->Write();

  hpt_gen_0->Write();
  hpt_gen_1->Write();
  hpt_gen_2->Write();

  hInt_reco_1->Write();
  hInt_reco_2->Write();

  hInt_gen_1->Write();
  hInt_gen_2->Write();
  hy_gen->Write();
  hy_reco->Write();

  outFile->Close();

}
