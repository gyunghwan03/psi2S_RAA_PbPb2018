#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
//#include "Style_jaebeom.h"
//#include "tnp_weight_lowptPbPb.h"
#include "../tnp_weight_lowptPbPb_num_den_new.h"
#include <TAttMarker.h>

using namespace std;
// v5 : remove weighting factor systematics

void get_Eff_psi_pbpb_ctauCut(
  int state=1,
  float ptLow = 3.0, float ptHigh = 40.0,
  float yLow = 0.0, float yHigh = 2.4,
  //float cLow = 0, float cHigh = 20, 
  float cLow = 0, float cHigh = 180, 
  //bool isTnP = false, bool isPtWeight = false, int state=1
  //bool isTnP = true, bool isPtWeight = false, int state=1
  //bool isTnP = false, bool isPtWeight = true, int state=1
  bool isTnP = true, bool isPtWeight = true
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 12;	//jpsi=12,Upsilon=13

  float muPtCut = 0; //3.5, 1.8

  //jpsi
  //float massLow = 0.0;
  //float massHigh = 10.0;
  float massLow = 3.4;
  float massHigh = 4.1;

  double min = 0;
  double max = ptHigh;
  const int numBins = 9; //40;//(max-min)/binwidth;  //31//


  //input files
  //PbPb
  TString inputMC1 = "/data/Oniatree/miniAOD/OniaTree_Run2_PbPb2018MC_Prompt_forPsi2S_miniAOD112X_muonSelGlb_13May23v2.root";
  if(state==2) inputMC1 = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_psi2S_BtoPsi2S_pThat2_HydjetMB_Glb_5p02TeV_merged.root";	//PbPb_non prompt
  TChain* mytree = new TChain("hionia/myTree"); 
  mytree->Add(inputMC1.Data());
  //TFile *inf = new TFile("../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root","READ");
  //TTree *mytree = (TTree*)inf->Get("myTree");
  //if(state==1) mytree->Add(inputMC2.Data());
  cout<<"dmoon chk entries : "<<mytree->GetEntries()<<endl;
  //TString outFileName = "mc_eff_vs_pt_cent_prompt_pbpb_Jpsi.root"; 
  //if(state==2) outFileName = "mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi.root"; 

  //pT reweighting function
  //ratioDataMC_AA_Jpsi_DATA_y0.0-2.4_210915.root
  TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_psi2S_DATA_ctauCut_y0_1p6_251103.root","read");
  TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_psi2S_DATA_ctauCut_y1p6_2p4_251103.root","read");
  //TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root","read");
  //TFile *fPtW2 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
  //if(state==2) TFile *fPtW = new TFile("ratioDataMC_AA_btojpsi_DATA_1s.root","read");
  TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
  TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");


  double ptBin_for[6] = {0,3.5,6.5,9,12,40};
  double ptBin_mid[8] = {0,6.5,9,12,15,20,25,40};
  double centBin_for[5] = {0,20,60,100,180};
  double centBin_mid[7] = {0,20,40,60,80,100,180};
  float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH1D* hpt_reco_1 = new TH1D("hpt_reco_1","hpt_reco_1",5,ptBin_for);
  TH1D* hpt_reco_2 = new TH1D("hpt_reco_2","hpt_reco_2",7,ptBin_mid);

  TH1D* hpt_gen_1 = new TH1D("hpt_gen_1","hpt_gen_1",5,ptBin_for);
  TH1D* hpt_gen_2 = new TH1D("hpt_gen_2","hpt_gen_2",7,ptBin_mid);

  TH1D* hcent_reco_1 = new TH1D("hcent_reco_1", "hcent_reco_1",4,centBin_for);
  TH1D* hcent_reco_2 = new TH1D("hcent_reco_2", "hcent_reco_2",6,centBin_mid);

  TH1D* hcent_gen_1 = new TH1D("hcent_gen_1", "hcent_gen_1",4,centBin_for);
  TH1D* hcent_gen_2 = new TH1D("hcent_gen_2", "hcent_gen_2",6,centBin_mid);

  TH1D* hInt_gen_1 = new TH1D("hInt_gen_1", "",1,0,40);
  TH1D* hInt_gen_2 = new TH1D("hInt_gen_2", "",1,0,40);

  TH1D* hInt_reco_1 = new TH1D("hInt_reco_1", "",1,0,40);
  TH1D* hInt_reco_2 = new TH1D("hInt_reco_2", "",1,0,40);

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",3,40);
  f1->SetParameters(214,22,14);
  f1->SetParLimits(0,0,500);
  f1->SetParLimits(1,-1,500);
  f1->SetParLimits(2,0,500);


  TH1D* hy_gen = new TH1D("hy_gen","hy_gen",6,yBin);
  TH1D* hy_reco = new TH1D("hy_reco","hy_reco",6,yBin);
  hy_gen ->Sumw2();
  hy_reco ->Sumw2();

  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();

  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();


  const int maxBranchSize = 1000;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Short_t           Gen_QQ_size;
  Short_t           Gen_mu_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Gen_QQ_trig[maxBranchSize];   //[Gen_QQ_size]
  Float_t         Gen_QQ_VtxProb[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Centrality;   //!
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
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
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

  int kL2filter = 16;	//jpsi=16,Upsilon=38
  int kL3filter = 17;	//jpsi=17,Upsilon=39

  int count =0;
  int counttnp =0;
  int nevt = mytree->GetEntries();
  //const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  //nevt = 100;
  for(int iev=0; iev<nevt ; ++iev)
  //for(int iev=0; iev<300000 ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    weight = findNcoll(Centrality) * Gen_weight;
	///Gen_QQ_size = 1;
	//if(Gen_QQ_size > 0) cout<<"weight : "<<weight<<", Gen_QQ_size : "<<Gen_QQ_size<<", Reco_QQ_size : "<<Reco_QQ_size<<endl;

    for(int igen = 0; igen<Gen_QQ_size; igen++){
	    JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
	    mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
	    mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

	    Double_t Rapidity_g = fabs(JP_Gen->Rapidity());
		//cout<<"Rapidity_g : "<<Rapidity_g<<", pt : "<<JP_Gen->Pt()<<endl;
        if(! (JP_Gen->M()>massLow && JP_Gen->M()<massHigh) ) continue;

	    if(! (JP_Gen->Pt() > 0 && JP_Gen->Pt()<40 && fabs(JP_Gen->Rapidity())<2.4 && IsAcceptanceQQ(mupl_Gen->Pt(),fabs(mupl_Gen->Eta()))&&IsAcceptanceQQ(mumi_Gen->Pt(),fabs(mumi_Gen->Eta()))) ) continue;

	    if(! ( fabs(JP_Gen->Rapidity())<2.4 && fabs(mupl_Gen->Eta())<2.4 && fabs(mumi_Gen->Eta())<2.4) )  continue;
	    if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

	    pt_weight = 1;
	    //if(isPtWeight) pt_weight = fptw->Eval(JP_Gen->Pt());
		if(isPtWeight && fabs(JP_Gen->Rapidity()) < 1.6 ) pt_weight = fptw1->Eval(JP_Gen->Pt());
		if(isPtWeight && fabs(JP_Gen->Rapidity()) > 1.6 && fabs(JP_Gen->Rapidity()) < 2.4 ) pt_weight = fptw2->Eval(JP_Gen->Pt());
	    if(JP_Gen->Pt() > 6.5) hy_gen->Fill(Rapidity_g, weight*pt_weight);
	    if( (Centrality > cLow && Centrality < cHigh)) { 
			if(Rapidity_g > 1.6 && Rapidity_g <2.4){ hpt_gen_1->Fill(JP_Gen->Pt(),weight*pt_weight);}
			if(Rapidity_g < 1.6 ) hpt_gen_2->Fill(JP_Gen->Pt(),weight*pt_weight); }
		if(Rapidity_g > 1.6 && Rapidity_g <2.4){
			//hcent_gen_1->Fill( Centrality,weight*pt_weight);
			if(JP_Gen->Pt() > 3.5 && JP_Gen->Pt() < 40) { hInt_gen_1->Fill(1,weight*pt_weight); hcent_gen_1->Fill( Centrality,weight*pt_weight); }
		}
		if(Rapidity_g < 1.6 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 40 ){
			hcent_gen_2->Fill( Centrality,weight*pt_weight);
			hInt_gen_2->Fill(1,weight*pt_weight);
		}
		


    }


    bool HLTPass = false;
    if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTPass=true;

    for (Short_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
	    JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
	    mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
	    mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
		//cout<<"Reco QQ pt : "<<JP_Reco->Pt()<<endl;

      
      bool HLTFilterPass=false;
      if( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTFilterPass=true;

      if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      bool muplSoft = (  //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl        //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      if ( !(muplSoft && mumiSoft) ) continue;   
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;
      if(Reco_QQ_sign[irqq]!=0) continue;  

      Double_t Rapidity = fabs(JP_Reco->Rapidity());

      if(! (JP_Reco->Pt()>3.5&&JP_Reco->Pt()<40&&fabs(JP_Reco->Rapidity())<2.4&&IsAcceptanceQQ(mupl_Reco->Pt(),fabs(mupl_Reco->Eta()))&&IsAcceptanceQQ(mumi_Reco->Pt(),fabs(mumi_Reco->Eta()))) ) continue;

      if(!( fabs(mupl_Reco->Eta())<2.4 && fabs(mumi_Reco->Eta())<2.4 && fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->M()>massLow && JP_Reco->M()<massHigh)) continue;


      if(HLTPass==true && HLTFilterPass==true) count++;
      if(isTnP){
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


       tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); //mu id
       tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0); //inner tracker

       //Trigger part
       if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ){
//         cout << "irqq : " << irqq << " - iev : " << iev << endl;
//         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
         continue;
       }
       bool mupl_L2Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
       bool mupl_L3Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
       bool mumi_L2Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
       bool mumi_L3Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
       if(mupl_L2Filter == false || mumi_L2Filter == false){ cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl; cout << endl;} 

       bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
       bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
       bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
       bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
       bool SelDone = false;

       //if( mupl_isL2 && mumi_isL3){
       //        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
       //        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
       //        SelDone = true;
       //}
       //else if( mupl_isL3 && mumi_isL2){
       //        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
       //        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
       //        SelDone = true;
       //}
       //else if( mupl_isL3 && mumi_isL3){
       //        int t[2] = {-1,1}; // mupl, mumi
       //        int l = rand() % (2); 
       //        //pick up what will be L2
       //        if(t[l]==-1){
       // 	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
       // 	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
       //        }
       //        else if(t[l]==1){
       // 	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
       // 	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
       //        }
       //        else {cout << "ERROR :: No random selection done !!!!" << endl; continue;}
       //        SelDone = true;
       //}    

       if( mupl_isL2 && mumi_isL3){
	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
	       SelDone = true;
	       tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;

       }
       else if( mupl_isL3 && mumi_isL2){
	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
	       SelDone = true;
	       tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
       }
       else if( mupl_isL3 && mumi_isL3){

	       tnp_trig_weight_muplL2_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
	       tnp_trig_weight_muplL3_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
	       tnp_trig_weight_mumiL2_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
	       tnp_trig_weight_mumiL3_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

	       tnp_trig_weight_muplL2_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
	       tnp_trig_weight_muplL3_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
	       tnp_trig_weight_mumiL2_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
	       tnp_trig_weight_mumiL3_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

	       tnp_trig_weight_num = tnp_trig_weight_muplL2_num * tnp_trig_weight_mumiL3_num + tnp_trig_weight_mumiL2_num * tnp_trig_weight_muplL3_num - tnp_trig_weight_muplL3_num * tnp_trig_weight_mumiL3_num;
	       tnp_trig_weight_den = tnp_trig_weight_muplL2_den * tnp_trig_weight_mumiL3_den + tnp_trig_weight_mumiL2_den * tnp_trig_weight_muplL3_den - tnp_trig_weight_muplL3_den * tnp_trig_weight_mumiL3_den;
	       tnp_trig_weight = tnp_trig_weight_num/tnp_trig_weight_den;
       }

       tnp_weight = tnp_weight * tnp_trig_weight;



       //if(SelDone == false){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
       //if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
       if(HLTPass==true && HLTFilterPass==true){
	       counttnp++;
	       tnp_trig_dimu = tnp_trig_weight;
	       //hpt_tnp_trig->Fill(JP_Reco->Pt(),tnp_trig_dimu);
       }
      }

      pt_weight = 1;
      //if(isPtWeight) pt_weight = fptw->Eval(JP_Reco->Pt());
	  if(isPtWeight && fabs(JP_Reco->Rapidity()) < 1.6) pt_weight = fptw1->Eval(JP_Reco->Pt());
      if(isPtWeight && fabs(JP_Reco->Rapidity()) > 1.6 && fabs(JP_Reco->Rapidity()) < 2.4) pt_weight = fptw2->Eval(JP_Reco->Pt());

      //cout<<"pt_Weight in reco : "<<pt_weight<<endl;
      if(HLTPass==true && HLTFilterPass==true){
          if(JP_Reco->Pt() > 6.5 ) hy_reco->Fill(Rapidity, weight* tnp_weight* pt_weight);
          
          if( (Centrality > cLow && Centrality < cHigh)) {
		  if(Rapidity > 1.6 && Rapidity < 2.4) hpt_reco_1->Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
          if(Rapidity < 1.6) hpt_reco_2->Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight); }

		  if(Rapidity > 1.6 && Rapidity <2.4){
			  //hcent_reco_1->Fill( Centrality,weight*tnp_weight*pt_weight);
			  if (JP_Reco->Pt() > 3.5 && JP_Reco->Pt() < 40 ) { hInt_reco_1->Fill(1,weight*tnp_weight*pt_weight); hcent_reco_1->Fill( Centrality,weight*tnp_weight*pt_weight); }
		  }
		  if(Rapidity < 1.6 && JP_Reco->Pt() > 6.5 && JP_Reco->Pt() < 40 ){
			  hcent_reco_2->Fill( Centrality,weight*tnp_weight*pt_weight);
			  hInt_reco_2->Fill(1,weight*tnp_weight*pt_weight);
		  }

      }
    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  

  //Divide
  TH1D* hpt_eff_1;
  TH1D* hpt_eff_2;

  TH1D* hcent_eff_1;
  TH1D* hcent_eff_2;

  TH1D* hInt_eff_1;
  TH1D *hInt_eff_2;

  hpt_eff_1 = (TH1D*)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D*)hpt_reco_2->Clone("hpt_eff_2");


  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1, "B");
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1, "B");

  hcent_eff_1 = (TH1D*)hcent_reco_1->Clone("hcent_eff_1");
  hcent_eff_2 = (TH1D*)hcent_reco_2->Clone("hcent_eff_2");


  hcent_eff_1->Divide(hcent_eff_1, hcent_gen_1, 1, 1, "B");
  hcent_eff_2->Divide(hcent_eff_2, hcent_gen_2, 1, 1, "B");
  
  hInt_eff_1 = (TH1D*)hInt_reco_1->Clone("hInt_eff_1");
  hInt_eff_2 = (TH1D*)hInt_reco_2->Clone("hInt_eff_2");


  hInt_eff_1->Divide(hInt_eff_1, hInt_gen_1, 1, 1, "B");
  hInt_eff_2->Divide(hInt_eff_2, hInt_gen_2, 1, 1, "B");

  TH1D* hy_eff;
  hy_eff = (TH1D*)hy_reco->Clone("hy_eff");
  hy_eff->Divide(hy_eff, hy_gen, 1, 1, "B");

  hpt_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hpt_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");

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


//draw same
  //gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend = new TLegend(0.6,0.84);
  TLatex *lt2 = new TLatex();
  lt2->SetNDC();
  lt2->SetTextSize(0.03);
  //auto legend2 = new TLegend(0.6,0.84);

  gStyle->SetOptFit(0);
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
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #psi(2S) (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #psi(2S) (PbPb)");
  //drawsame
  hpt_eff_2->SetMarkerStyle(25);
  hpt_eff_2->SetMarkerColor(kRed+1);
  hpt_eff_2->SetLineColor(kRed+1);
  hpt_eff_2->Draw("same");
  legend->AddEntry("hpt_eff_1",Form("|y|: 1.6-2.4, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_2",Form("|y|: 0-1.6, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->SetBorderSize( 0);
  legend->Draw("E");
  cpt_eff->SaveAs("Eff_pt_noweight2.png");

  TCanvas * ccent_eff = new TCanvas("ccent_eff","ccent_eff",0,0,900,800);
  ccent_eff->cd();
  hcent_eff_1->SetMarkerStyle(24);
  hcent_eff_1->SetMarkerColor(1);
  hcent_eff_1->SetLineColor(1);
  hcent_eff_1->GetXaxis()->SetTitle("Centrality [%]");
  hcent_eff_1->GetYaxis()->SetTitle("Efficiency");
  hcent_eff_1->GetYaxis()->SetRangeUser(0.,1.3);
  hcent_eff_1->Draw("E");
  hcent_eff_2->SetMarkerStyle(25);
  hcent_eff_2->SetMarkerColor(kRed+1);
  hcent_eff_2->SetLineColor(kRed+1);
  hcent_eff_2->Draw("same");

  lt2->SetTextSize(0.03);
  if(state==1) lt2->DrawLatex(0.13,0.84,"Prompt #psi(2S) (PbPb)");
  else if(state==2) lt2->DrawLatex(0.13,0.84,"Non-Prompt #psi(2S) (PbPb)");

  legend->AddEntry("hcent_eff_1",Form("|y|: 1.6-2.4, p_{T} : 3-40 GeV/c"),"lep");
  legend->AddEntry("hcent_eff_2",Form("|y|: 0-1.6, p_{T} : 6.5-40 GeV/c"),"lep");
  legend->SetBorderSize( 0);
  legend->Draw("E");
  ccent_eff->SaveAs("./figs/Eff_cent_noweight2.png");

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
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #psi(2S) (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #psi(2S) (PbPb)");
  cy_eff->SaveAs("./figs/Eff_absy_noweight2.png");

  //Save efficiency files for later use.

  hpt_eff_1 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_2 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6  ",isTnP, isPtWeight, cLow, cHigh));
  hcent_eff_1 ->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_3_to_40_absy1p6_2p4",isTnP, isPtWeight));
  hcent_eff_2 ->SetName(Form("mc_eff_vs_cent_TnP%d_PtW%d_pt_6p5_to_40_absy0_1p6  ",isTnP, isPtWeight));
  hy_eff ->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d",isTnP, isPtWeight));
  hInt_eff_1 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_2p4",isTnP, isPtWeight, cLow, cHigh));
  hInt_eff_2 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_1p6  ",isTnP, isPtWeight, cLow, cHigh));

  //TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_drawsame1.root",cLow,cHigh,isPtWeight,isTnP);
  TString outFileName = Form("./roots/mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_ctauCut_251103.root",cLow,cHigh,isPtWeight,isTnP);
  if(state==2) outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_nprompt_pbpb_psi2S_PtW%d_tnp%d_ctauCut_20251103.root",cLow,cHigh,isPtWeight,isTnP);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hcent_eff_1->Write();
  hcent_eff_2->Write();
  hInt_eff_1->Write();
  hInt_eff_2->Write();
  hy_eff->Write();
  cpt_eff->Write();
  //if(isTnP) hpt_tnp_trig->Write();

  hpt_reco_1->Write();
  hpt_reco_2->Write();

  hpt_gen_1->Write();
  hpt_gen_2->Write();

  hcent_reco_1->Write();
  hcent_reco_2->Write();

  hcent_gen_1->Write();
  hcent_gen_2->Write();

  hInt_reco_1->Write();
  hInt_reco_2->Write();

  hInt_gen_1->Write();
  hInt_gen_2->Write();

  hy_gen->Write();
  hy_reco->Write();

  outFile->Close();

}
