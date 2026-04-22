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
// v5 : remove weighting factor systematics

void get_Eff_psi2S_pp_ctauCut_NP(
  float ptLow = 3.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  //float cLow = 0, float cHigh = 20, 
  //float cLow = 0, float cHigh = 200, 
  //bool isTnP = false, bool isPtWeight = false, int state=1
  //bool isTnP = true, bool isPtWeight = false, int state=1
  //bool isTnP = false, bool isPtWeight = true, int state=1
  bool isTnP = true, bool isPtWeight = true,
  int state=2
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 3;	//jpsi=12,Upsilon=13 pp Jpsi=3

  float muPtCut = 0; //3.5, 1.8

  //psi2S mass window
  float massLow = 3.4;
  float massHigh = 4.1;

  double min = 0;
  double max = ptHigh;
  const int numBins = 9; //50;//(max-min)/binwidth;  //31//


  //input files
  //pp
  TString inputMC1 = "/data/Oniatree/miniAOD/OniatreeMC_BPsi2SMM_TuneCUETP8M1_5p02TeV_pythia8_ptHatMin2_ONIATREE.root";	//pp_non prompt
  TChain* mytree = new TChain("hionia/myTree"); 
  mytree->Add(inputMC1.Data());
  cout<<"dmoon chk entries : "<<mytree->GetEntries()<<endl;

  //pT reweighting function
  TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoPsi2S_DATA_ctauCut_y0_1p6_251103.root","read");
  TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoPsi2S_DATA_ctauCut_y1p6_2p4_251103.root","read");
  TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
  TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");


  double ptBin_for[6] = {0,3.5,6.5,9,12,40};
  double ptBin_mid[9] = {0,6.5,9,12,15,20,25,30,40};
  double centBin_for[4] = {0,40,80,180};
  double centBin_mid[7] = {0,10,20,30,40,50,90};
  float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH1D* hpt_reco_1 = new TH1D("hpt_reco_1","hpt_reco_1",5,ptBin_for);
  TH1D* hpt_reco_2 = new TH1D("hpt_reco_2","hpt_reco_2",8,ptBin_mid);

  TH1D* hpt_gen_1 = new TH1D("hpt_gen_1","hpt_gen_1",5,ptBin_for);
  TH1D* hpt_gen_2 = new TH1D("hpt_gen_2","hpt_gen_2",8,ptBin_mid);

  TH1D* hInt_gen_1 = new TH1D("hInt_gen_1", "",1,0,50);
  TH1D* hInt_gen_2 = new TH1D("hInt_gen_2", "",1,0,50);

  TH1D* hInt_reco_1 = new TH1D("hInt_reco_1", "",1,0,50);
  TH1D* hInt_reco_2 = new TH1D("hInt_reco_2", "",1,0,50);

  TH1D* hy_gen = new TH1D("hy_gen","hy_gen",6,yBin);
  TH1D* hy_reco = new TH1D("hy_reco","hy_reco",6,yBin);
  hy_gen ->Sumw2();
  hy_reco ->Sumw2();

  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();

  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();

  hInt_gen_1->Sumw2();
  hInt_gen_2->Sumw2();

  hInt_reco_1->Sumw2();
  hInt_reco_2->Sumw2();


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

  // ctau3D variables for PbPb-like cuts
  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctau3DErr[maxBranchSize];
  TBranch *b_Reco_QQ_ctau3D;
  TBranch *b_Reco_QQ_ctau3DErr;
  mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  mytree->SetBranchAddress("Reco_QQ_ctau3DErr", Reco_QQ_ctau3DErr, &b_Reco_QQ_ctau3DErr);

  
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

  int kL2filter = 16;	//jpsi=16,Upsilon=38
  int kL3filter = 17;	//jpsi=17,Upsilon=39

  int counts[10]={0,0,0,0,0,0,0,0,0,0};
  int count =0;
  int counttnp =0;
  int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    weight = Gen_weight;  // pp doesn't use Ncoll

    for(int igen = 0; igen<Gen_QQ_size; igen++){
	    JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
	    mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
	    mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

	    Double_t Rapidity_g = fabs(JP_Gen->Rapidity());
        if(! (JP_Gen->M()>massLow && JP_Gen->M()<massHigh) ) continue;

	    if(! (JP_Gen->Pt() > 0 && JP_Gen->Pt()<50 && fabs(JP_Gen->Rapidity())<2.4 && IsAcceptanceQQ(mupl_Gen->Pt(),fabs(mupl_Gen->Eta()))&&IsAcceptanceQQ(mumi_Gen->Pt(),fabs(mumi_Gen->Eta()))) ) continue;

	    if(! ( fabs(JP_Gen->Rapidity())<2.4 && fabs(mupl_Gen->Eta())<2.4 && fabs(mumi_Gen->Eta())<2.4) )  continue;
	    if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

	    pt_weight = 1;
		if(isPtWeight && fabs(JP_Gen->Rapidity()) < 1.6 ) pt_weight = fptw1->Eval(JP_Gen->Pt());
		if(isPtWeight && fabs(JP_Gen->Rapidity()) > 1.6 && fabs(JP_Gen->Rapidity()) < 2.4 ) pt_weight = fptw2->Eval(JP_Gen->Pt());
	    
		if(JP_Gen->Pt() > 6.5) hy_gen->Fill(Rapidity_g, weight*pt_weight);
	    
		if(Rapidity_g > 1.6 && Rapidity_g <2.4 && JP_Gen->Pt()>3.5 && JP_Gen->Pt()<50){ 
			hpt_gen_1->Fill(JP_Gen->Pt(),weight*pt_weight); 
			counts[5]++;
		}
		if(Rapidity_g < 1.6 && JP_Gen->Pt()>6.5 && JP_Gen->Pt()<50 ) { 
			hpt_gen_2->Fill(JP_Gen->Pt(),weight*pt_weight); 
		}
		
		if(Rapidity_g > 1.6 && Rapidity_g <2.4 && JP_Gen->Pt()>3.5 && JP_Gen->Pt()<50){
			hInt_gen_1->Fill(1,weight*pt_weight); 
		}
		if(Rapidity_g < 1.6 && JP_Gen->Pt() > 6.5 && JP_Gen->Pt() < 50){
			hInt_gen_2->Fill(1,weight*pt_weight);
		}
    }


    bool HLTPass = false;
    if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTPass=true;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
	    JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
	    mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
	    mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      
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

      bool muplSoft = (  
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl
          ) ; 

      bool mumiSoft = ( 
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi
          ) ; 

      if ( !(muplSoft && mumiSoft) ) continue;   
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;
      if(Reco_QQ_sign[irqq]!=0) continue;  

      Double_t Rapidity = fabs(JP_Reco->Rapidity());

      if(! (JP_Reco->Pt()<50&&fabs(JP_Reco->Rapidity())<2.4&&IsAcceptanceQQ(mupl_Reco->Pt(),fabs(mupl_Reco->Eta()))&&IsAcceptanceQQ(mumi_Reco->Pt(),fabs(mumi_Reco->Eta()))) ) continue;

      if(!( fabs(mupl_Reco->Eta())<2.4 && fabs(mumi_Reco->Eta())<2.4 && fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->M()>massLow && JP_Reco->M()<massHigh)) continue;

      // Apply PbPb-like ctau3D cuts
      float ctauCut3D = 1000;
      float ctauCut3Derr = 3;
      float mass = JP_Reco->M();
      float pt = JP_Reco->Pt();
      
      // Forward rapidity cuts
      if(Rapidity > 1.6 && Rapidity < 2.4){
        if(pt > 3.5 && pt < 6.5){
          if(mass > 3.5 && mass < 4.0) ctauCut3D = 0.023;
        } else if(pt > 6.5 && pt < 50){
          if(mass > 3.5 && mass < 4.0) ctauCut3D = 0.014;
        }
      }
      // Mid rapidity cuts  
      else if(Rapidity < 1.6){
        if(pt > 6.5 && pt < 9){
          if(mass > 3.5 && mass < 4.0) ctauCut3D = 0.019;
        } else if(pt > 9 && pt < 12){
          if(mass > 3.5 && mass < 4.0) ctauCut3D = 0.016;
        } else if(pt > 12 && pt < 50){
          if(mass > 3.5 && mass < 4.0) ctauCut3D = 0.014;
        }
      }
      
      // Apply ctau3D cut
      if(fabs(Reco_QQ_ctau3D[irqq]) > ctauCut3D) continue;
      if(fabs(Reco_QQ_ctau3D[irqq]/Reco_QQ_ctau3DErr[irqq]) > ctauCut3Derr) continue;


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


       tnp_weight = tnp_weight * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mupl_Reco->Pt(), mupl_Reco->Eta())) * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mumi_Reco->Pt(), mumi_Reco->Eta())); //mu id

       //Trigger part
       if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ){
         continue;
       }

       tnp_weight = tnp_weight * tnp_trig_weight;

       if(HLTPass==true && HLTFilterPass==true){
	       counttnp++;
	       tnp_trig_dimu = tnp_trig_weight;
       }
      }
	  weight = 1;
      pt_weight = 1;
      if(isPtWeight && fabs(JP_Reco->Rapidity()) < 1.6) pt_weight = fptw1->Eval(JP_Reco->Pt());
      if(isPtWeight && fabs(JP_Reco->Rapidity()) > 1.6 && fabs(JP_Reco->Rapidity()) < 2.4) pt_weight = fptw2->Eval(JP_Reco->Pt());

      if(HLTPass==true && HLTFilterPass==true){
		  if(Rapidity > 1.6 && Rapidity < 2.4 && JP_Reco->Pt()>3.5 && JP_Reco->Pt()<50) { 
			  hpt_reco_1->Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight); 
			  hInt_reco_1->Fill(1, weight*tnp_weight*pt_weight); 
		  }
	      if(Rapidity < 1.6 && JP_Reco->Pt()>6.5 && JP_Reco->Pt()<50) { 
			  hpt_reco_2->Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight); 
			  hInt_reco_2->Fill(1, weight*tnp_weight*pt_weight); 
		  }
		  if(JP_Reco->Pt() > 6.5 ) hy_reco->Fill(Rapidity, weight* tnp_weight* pt_weight);
      }
    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  

  //Divide
  TH1D* hpt_eff_1;
  TH1D* hpt_eff_2;


  hpt_eff_1 = (TH1D*)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D*)hpt_reco_2->Clone("hpt_eff_2");


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

  hpt_eff_1->SetTitle("Eff: Rapidity 1.6-2.4");
  hpt_eff_2->SetTitle("Eff: Rapidity 0.0-1.6");

  
  //gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend = new TLegend(0.6,0.84);

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
  
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #psi(2S) (pp)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #psi(2S) (pp)");
  
  hpt_eff_2->SetMarkerStyle(25);
  hpt_eff_2->SetMarkerColor(kRed+1);
  hpt_eff_2->SetLineColor(kRed+1);
  hpt_eff_2->Draw("same");
  legend->AddEntry("hpt_eff_1","|y|: 1.6-2.4","lep");
  legend->AddEntry("hpt_eff_2","|y|: 0-1.6","lep");
  legend->SetBorderSize( 0);
  legend->Draw("E");
  cpt_eff->SaveAs("./figs_2S_pp/Eff_pt_ctauCut_251103.png");

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
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt #psi(2S) (pp)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt #psi(2S) (pp)");
  cy_eff->SaveAs("./figs_2S_pp/Eff_absy_ctauCut_251103.png");

  //Save efficiency files for later use.

  hpt_eff_1 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_absy1p6_2p4",isTnP, isPtWeight));
  hpt_eff_2 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_absy0_1p6",isTnP, isPtWeight));
  hInt_eff_1 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_absy1p6_2p4",isTnP, isPtWeight));
  hInt_eff_2 ->SetName(Form("mc_eff_Integrated_TnP%d_PtW%d_absy0_1p6",isTnP, isPtWeight));
  hy_eff ->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d",isTnP, isPtWeight));

  TString outFileName = Form("./roots_2S_pp/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW%d_tnp%d_ctauCut_251103.root",isPtWeight,isTnP);
  if(state==2) outFileName = Form("./roots_2S_pp/mc_eff_vs_pt_rap_nprompt_pp_psi2S_PtW%d_tnp%d_ctauCut_251103.root",isPtWeight,isTnP);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hInt_eff_1->Write();
  hInt_eff_2->Write();
  hy_eff->Write();
  cpt_eff->Write();

  hpt_reco_1->Write();
  hpt_reco_2->Write();

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
