#include <iostream>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBin.h"
//#include "Style.h"
//#include "TreeSetting.h"
#include "tnp_weight_lowptPbPb.h"

using namespace std;

void makeMuonSkimTree(bool isMC = true, bool isPR = true,
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true, bool isSwitch=false, int kTrigSel = kTrigUps, 
  int fitErrVar = 0, int tnp_ids = 0, int tnp_trks = 0, int tnp_trgs = 0 
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel_;
  if(kTrigSel == kTrigJpsi) kTrigSel_ = 0;
  if(kTrigSel == kTrigUps) kTrigSel_ = 1;
  else if(kTrigSel == kTrigL1DBOS40100) kTrigSel_ = 2; 

  float massLow = 2.0;
  float massHigh = 4.8;

  int CentSwitch = 87;

  string ftrigSel = fTrigName[kTrigSel_].Data();
//  if(isSwitch) ftrigSel += "_UpsAndL1OS40100";
//  else if(!isSwitch) ftrigSel += Form("_%s",fTrigName[kTrigSel_].Data());

  //input files
//  TString mcfilename = "Oniatree_Y3SMC_v210618_pdgfix.root";
//  if(state==1) mcfilename = "OniaTree_Ups1SMC_HydjetDrumMB_5p02TeV_v210621.root";
//  else if(state==2) mcfilename = "OniaTree_Ups2SMC_HydjetDrumMB_5p02TeV_v210628.root";
  TString inputMC = (isPR) ? "/eos/cms/store/group/phys_heavyions/soohwan/OniaTree/OniaTree_Run2_PbPb2018MC_Prompt_forPsi2S_miniAOD112X_muonSelGlb_13May23v2.root" : "/eos/cms/store/group/phys_heavyions/soohwan/OniaTree/OniaTree_Run2_PbPb2018MC_BToPsi_forPsi2S_miniAOD112X_muonSelGlb_13May23v2.root";
  //TString inputMC = "/home/samba.old/CMS_Files/UpsilonAnalysis/Ups3S_PbPb2018/OniaTree/MC/Oniatree_MC_103X_2021May.root";
  TString inputData = "/home/CMS/DataFiles/PbPb2018/OniaTree/MC/OniaTree_MC_Ups3S_PbPb2018_HydjetDrumMB_5p02TeV_merged.root";
  TChain* mytree = new TChain("hionia/myTree"); 
  if(isMC) mytree->Add(inputMC.Data());
  else if(!isMC) mytree->Add(inputData.Data());

  //SetBranchAddress
  //SetTree settree_;
  //settree_.TreeSetting(mytree,isMC,0);


  //pT reweighting function
  string rapDef = (yHigh == 2.4) ? "Forward_y_211218" : "y0_1p6_211201";
  if( !( yHigh== float(2.4) || yHigh == float(1.6) ) ){
	  std::cout << "rapidity range does not correspond to neither forward nor midrange definition!\n";
	  return;
  }

  TFile *fPtW = new TFile(Form("Eff_Acc/ratioDataMC_AA_Jpsi_DATA_%s.root", rapDef.c_str()),"read");
  TF1* f1 = (TF1*) fPtW->Get("dataMC_Ratio1");

  string fitVarStr = "Nom";
  if(fitErrVar==1){
    for(int ip=0; ip<f1->GetNpar(); ip++){ f1->SetParameter(ip, f1->GetParameter(ip) + f1->GetParError(ip));}
    fitVarStr = "Up";
  }
  else if(fitErrVar==-1){
    for(int ip=0; ip<f1->GetNpar(); ip++){ f1->SetParameter(ip, f1->GetParameter(ip) - f1->GetParError(ip));}
    fitVarStr = "Do";
  }

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Float_t         SumET_HF;
  Short_t           Reco_QQ_size;
  Short_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  Short_t         Reco_QQ_whichGen[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_SumET_HF;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Reco_QQ_whichGen;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
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

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Short_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  //Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  //TBranch        *b_Reco_mu_ptErr_global;   //!
  //mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Float_t         Reco_QQ_ctau3D[maxBranchSize];
  Float_t         Reco_QQ_ctauErr3D[maxBranchSize];
  TBranch        *b_Reco_QQ_ctau3D;
  TBranch        *b_Reco_QQ_ctauErr3D;
  mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  mytree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  
  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if(isMC){
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
  }


  double mass;
  double pt, y, pt1, pt2, eta1, eta2, eta;
  int nTrkHits1, nTrkHits2;
  double normChi2_global1, normChi2_global2;
  int nMuValHits1, nMuValHits2;
  int StationsMatched1, StationsMatched2;
  double dxy1, dxy2, dxyErr1, dxyErr2, dz1, dz2, dzErr1, dzErr2;
  bool TMOneStaTight1, TMOneStaTight2;
  int nPixWMea1, nPixWMea2;
  int nTrkWMea1, nTrkWMea2;
  double ctau, ctau3D;
  double ctauErr, ctauErr3D;
  int nPixValHits1, nPixValHits2;
  double ptErr_global1, ptErr_global2;
  bool highPurity1, highPurity2;
  double QQVtxProb;
  double QQMassErr;
  double normChi2_inner1, normChi2_inner2;
  double ptErr_inner1, ptErr_inner2;
  int cBin;

  bool isCowboy;
  double cosAlpha;
  double cosAlpha3D;
  double QQdca;

  double phi1, phi2;
  double ptimb;
  double dcosMuTheta;
  
  double weight = 1;
  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double pt_weight = 1;

  double Ncoll_weight = 1;
  double Gen_weight_ = 1;

  int kL2filter = 16;
  int kL3filter = 17;

  string ftnpId = "Nom";
  string ftnpTrk = "Nom";
  string ftnpTrg = "Nom";

  if(tnp_ids==-1) ftnpId = "sysUp";
  else if(tnp_ids==-2) ftnpId = "sysDo";
  else if(tnp_ids==1) ftnpId = "statUp";
  else if(tnp_ids==2) ftnpId = "statDo";
  else if(tnp_ids==99) ftnpId = "TagChange";
  if(tnp_trks==-1) ftnpTrk = "sysUp";
  else if(tnp_trks==-2) ftnpTrk = "sysDo";
  else if(tnp_trks==1) ftnpTrk = "statUp";
  else if(tnp_trks==2) ftnpTrk = "statDo";
  else if(tnp_trks==99) ftnpTrk = "TagChange";
  if(tnp_trgs==-1) ftnpTrg = "sysUp";
  else if(tnp_trgs==-2) ftnpTrg = "sysDo";
  else if(tnp_trgs==1) ftnpTrg = "statUp";
  else if(tnp_trgs==2) ftnpTrg = "statDo";
  else if(tnp_trgs==99) ftnpTrg = "TagChange";
 
  string flavorTag = (isPR) ? "PR" : "NP";
  const char* outFileName = Form("TnPSkim/ThisTEST_OutputSkim_isMC%d_%s_fitVar_noFit%s_tnp_id%s_trk%s_trg%s.root",isMC, flavorTag.c_str(),fitVarStr.c_str(),ftnpId.c_str(), ftnpTrk.c_str(), ftnpTrg.c_str());
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TTree *wtree = new TTree("tree","");
  wtree->Branch("cBin",&cBin);
  wtree->Branch("mass",&mass);
  wtree->Branch("pt",&pt);
  wtree->Branch("y",&y);
  wtree->Branch("pt1",&pt1);
  wtree->Branch("pt2",&pt2);
  wtree->Branch("ptimb",&ptimb);
  wtree->Branch("eta1",&eta1);
  wtree->Branch("eta2",&eta2);
  wtree->Branch("phi1",&phi1);
  wtree->Branch("phi2",&phi2);
  wtree->Branch("dcosMuTheta",&dcosMuTheta);
  wtree->Branch("eta",&eta);
  wtree->Branch("nTrkHits1",&nTrkHits1);
  wtree->Branch("nTrkHits2",&nTrkHits2);
  wtree->Branch("normChi2_global1",&normChi2_global1);
  wtree->Branch("normChi2_global2",&normChi2_global2);
  wtree->Branch("nMuValHits1",&nMuValHits1);
  wtree->Branch("nMuValHits2",&nMuValHits2);
  wtree->Branch("StationsMatched1",&StationsMatched1);
  wtree->Branch("StationsMatched2",&StationsMatched2);
  wtree->Branch("dxy1",&dxy1);
  wtree->Branch("dxy2",&dxy2);
  wtree->Branch("dxyErr1",&dxyErr1);
  wtree->Branch("dxyErr2",&dxyErr2);
  wtree->Branch("dz1",&dz1);
  wtree->Branch("dz2",&dz2);
  wtree->Branch("dzErr1",&dzErr1);
  wtree->Branch("dzErr2",&dzErr2);
  wtree->Branch("nPixWMea1",&nPixWMea1);
  wtree->Branch("nPixWMea2",&nPixWMea2);
  wtree->Branch("nTrkWMea1",&nTrkWMea1);
  wtree->Branch("nTrkWMea2",&nTrkWMea2);
  wtree->Branch("ctau",&ctau);
  wtree->Branch("ctau3D",&ctau3D);
  wtree->Branch("ctauErr",&ctauErr);
  wtree->Branch("ctauErr3D",&ctauErr3D);
  wtree->Branch("cosAlpha",&cosAlpha);
  wtree->Branch("cosAlpha3D",&cosAlpha3D);
  wtree->Branch("QQdca",&QQdca);
  wtree->Branch("isCowboy",&isCowboy);
  wtree->Branch("nPixValHits1",&nPixValHits1);
  wtree->Branch("nPixValHits2",&nPixValHits2);
  wtree->Branch("highPurity1",&highPurity1);
  wtree->Branch("highPurity2",&highPurity2);
  wtree->Branch("QQVtxProb",&QQVtxProb);
  wtree->Branch("QQMassErr",&QQMassErr);
  wtree->Branch("normChi2_inner1",&normChi2_inner1);
  wtree->Branch("normChi2_inner2",&normChi2_inner2);
  wtree->Branch("ptErr_inner1",&ptErr_inner1);
  wtree->Branch("ptErr_inner2",&ptErr_inner2);
  wtree->Branch("weight",&weight);
  wtree->Branch("Ncoll_weight",&Ncoll_weight);
  wtree->Branch("Gen_weight",&Gen_weight_);
  wtree->Branch("tnp_weight",&tnp_weight);
  wtree->Branch("pt_weight",&pt_weight);
  
  if(!isMC){
    wtree->Branch("ptErr_global1",&ptErr_global1);
    wtree->Branch("ptErr_global2",&ptErr_global2);
    wtree->Branch("TMOneStaTight1",&TMOneStaTight1);
    wtree->Branch("TMOneStaTight2",&TMOneStaTight2);
  }


  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  
  double tnp_trig_dimu=-1;

  int count =0;
  int counttnp =0;
  const int nevt =mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    if(isMC) cBin = Centrality;
    else if(!isMC) cBin = getHiBinFromhiHF(SumET_HF);

    if(cBin > cHigh || cBin < cLow) continue;
    
    if(isMC){
      weight = findNcoll(Centrality) * Gen_weight;
      Ncoll_weight = findNcoll(Centrality);
      Gen_weight_ = Gen_weight;
    }
    else if(!isMC){weight = 1;}
    
    bool HLTPass = false;
    if(isSwitch){ 
      if(Centrality>=CentSwitch){kTrigSel = kTrigL1DBOS40100;}
      else if(Centrality<CentSwitch){kTrigSel = kTrigUps;}
    }

    if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTPass=true;
    if(HLTPass==false) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; irqq++) 
    {
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      bool HLTFilterPass=false;
      if( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTFilterPass=true;
      if(HLTFilterPass==false) continue;

      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_QQ_whichGen[irqq] == -1) continue;
      }
      if(Reco_QQ_sign[irqq]!=0) continue;  

//      if ( Reco_QQ_VtxProb[irqq]  < 0.005 ) continue;
      
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
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));
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

      if(isPtWeight){
        pt_weight = f1->Eval(JP_Reco->Pt());
        weight = weight * pt_weight;
      }

      if(isTnP){
        tnp_weight = 1;
        tnp_trig_weight_mupl = -1;
        tnp_trig_weight_mumi = -1;
        tnp_weight = tnp_weight * tnp_weight_muid_pbpb(pt1, eta1, tnp_ids) * tnp_weight_muid_pbpb(pt2, eta2, tnp_ids); //mu id
        tnp_weight = tnp_weight * tnp_weight_trk_pbpb(eta1, tnp_trks) * tnp_weight_trk_pbpb(eta2, tnp_trks); //inner tracker

        if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ) continue;

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

        if( mupl_isL2 && mumi_isL3){
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, tnp_trgs); 
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, tnp_trgs); 
          SelDone = true;
        }   
        else if( mupl_isL3 && mumi_isL2){
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, tnp_trgs); 
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, tnp_trgs); 
          SelDone = true;
        }   
        else if( mupl_isL3 && mumi_isL3){
          double T1_ = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, tnp_trgs, true);
          double T2_ = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, tnp_trgs, true);
          double T1 = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, tnp_trgs, true);
          double T2 = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, tnp_trgs, true);
          double den_ = T1_*T2 + (T1-T1_)*T2_;
          double num_ = T1_*tnp_weight_trg_pbpb(pt1, eta1, 3, tnp_trgs)*T2*tnp_weight_trg_pbpb(pt2, eta2, 2, tnp_trgs) + (T1*tnp_weight_trg_pbpb(pt1, eta1, 2, tnp_trgs)-T1_*tnp_weight_trg_pbpb(pt1, eta1, 3, tnp_trgs))*T2_*tnp_weight_trg_pbpb(pt2, eta2, 3, tnp_trgs);

          tnp_trig_weight_mupl = num_/den_; tnp_trig_weight_mumi = 1;
          if(den_<=0 || num_<=0){cout << "ERROR wrong calculation" << endl; continue;}

          SelDone = true;
        }
        if(SelDone == false){continue;}//cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;
        if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
        tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
      }

      weight = weight * tnp_weight;
      mass = JP_Reco->M();
      pt = JP_Reco->Pt();
      y = JP_Reco->Rapidity();
      eta = JP_Reco->Eta();
      nTrkHits1 = Reco_mu_nTrkHits[Reco_QQ_mupl_idx[irqq]];  
      nTrkHits2 = Reco_mu_nTrkHits[Reco_QQ_mumi_idx[irqq]];  
      normChi2_global1 = Reco_mu_normChi2_global[Reco_QQ_mupl_idx[irqq]];
      normChi2_global2 = Reco_mu_normChi2_global[Reco_QQ_mumi_idx[irqq]];
      nMuValHits1 = Reco_mu_nMuValHits[Reco_QQ_mupl_idx[irqq]];
      nMuValHits2 = Reco_mu_nMuValHits[Reco_QQ_mumi_idx[irqq]];
      nPixWMea1 = Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]];
      nPixWMea2 = Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]];
      nTrkWMea1 = Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]];
      nTrkWMea2 = Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]];
      StationsMatched1 = Reco_mu_StationsMatched[Reco_QQ_mupl_idx[irqq]];
      StationsMatched2 = Reco_mu_StationsMatched[Reco_QQ_mumi_idx[irqq]];
      ctau = 0;
      ctau3D = Reco_QQ_ctau3D[irqq];
      ctauErr = 0;
      ctauErr3D = Reco_QQ_ctauErr3D[irqq];
      cosAlpha = 0;
      cosAlpha3D = 0;
      QQdca = 0;
      isCowboy = 0;
      nPixValHits1 = Reco_mu_nPixValHits[Reco_QQ_mupl_idx[irqq]];
      nPixValHits2 = Reco_mu_nPixValHits[Reco_QQ_mumi_idx[irqq]];
      highPurity1 = Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]];
      highPurity2 = Reco_mu_highPurity[Reco_QQ_mumi_idx[irqq]];
      dxy1 = Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]];
      dxy2 = Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]];
      dxyErr1 = Reco_mu_dxyErr[Reco_QQ_mupl_idx[irqq]];
      dxyErr2 = Reco_mu_dxyErr[Reco_QQ_mumi_idx[irqq]];
      dz1 = Reco_mu_dz[Reco_QQ_mupl_idx[irqq]];
      dz2 = Reco_mu_dz[Reco_QQ_mumi_idx[irqq]];
      dzErr1 = Reco_mu_dzErr[Reco_QQ_mupl_idx[irqq]];
      dzErr2 = Reco_mu_dzErr[Reco_QQ_mumi_idx[irqq]];
      normChi2_inner1 = 0;
      normChi2_inner2 = 0;
      ptErr_inner1 = 0;
      ptErr_inner2 = 0;
      QQVtxProb = Reco_QQ_VtxProb[irqq];
      QQMassErr = 0;
      wtree->Fill();
    }
  }

  wtree->Write();
  outFile->Close();
}
