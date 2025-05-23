#include <ctime>
#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBin.h"
#include "tnp_weight_pp.h"

static const long MAXTREESIZE = 100000000000000;
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);


void SkimTree_Event_psi2S_pp_MC(int nevt=1000000, bool isMC = true, int kTrigSel = kTrigJpsipp, int hiHFBinEdge = -1, int PDtype = 1) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  //TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/MinBias/HIMinimumBias_Run2018_Upsilon_PromptReco_v1.root";
  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameDataReReco = "/work2/Oniatree/JPsi/OniaTree_DoubleMuPD_pp2017_5p02TeV_EOY_merged.root";
  TString fnameDataReRecoPeri = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuonPsiPeri/ReReco_Oniatree_addvn_part*.root";
  TString fnameMC = "/work2/Oniatree/Psi2S/OniatreeMC_Psi2SMM_TuneCUETP8M1_5p02TeV_pythia8_RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV-v2.root";

  TString fPD;
  if(PDtype==1) fPD = "DoubleMuon";
  else if(PDtype==2) fPD = "DBPeri";

  TChain *mytree = new TChain("hionia/myTree");
  if(!isMC){
    if(PDtype==1) mytree->Add(fnameDataReReco.Data());
    else if(PDtype==2) mytree->Add(fnameDataReRecoPeri.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  const long int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Float_t         SumET_HF;
  //int           Reco_QQ_size;
  Short_t           Reco_QQ_size;
  Short_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
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

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  //mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  //mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
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
  //Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  //TBranch        *b_Reco_mu_TMOneStaTight;   //!

  //mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
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


  //TChain *eptree = new TChain("tree");
  //if(!isMC){
  //  eptree->Add(fnameDataReReco.Data());
  //}
  //else if(isMC){
  //  eptree->Add(fnameMC.Data());
  //}
  
  //const int nEP = 29;  // number of event planes in the tree
  //double qx[nEP]; 
  //double qy[nEP]; 
  //TBranch *b_qx;
  //TBranch *b_qy;
  //eptree->SetBranchAddress("qx",qx, &b_qx);
  //eptree->SetBranchAddress("qy",qy, &b_qy);
  

  int trigIndx=0;
  if(kTrigSel == kTrigJpsipp) trigIndx=0; //PbPb : 0
  else if(kTrigSel == kTrigUps) trigIndx=1;
  else if(kTrigSel == kTrigL1DBOS40100) trigIndx=2;
  else if(kTrigSel == kTrigL1DB50100) trigIndx=3;
  
  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;

  //int kL2filter = 13;
  //int kL3filter = 13;

  int count =0;
  int counttnp =0;

  //TString fCentSelHF = "HFNom";
  //if(hiHFBinEdge==1) fCentSelHF = "HFUp";
  //else if(hiHFBinEdge==-1) fCentSelHF = "HFDown";
  TFile* newfile;
  newfile = new TFile(Form("OniaFlowSkim_%sTrig_%sPD_pp_psi2S_isMC%d_230126.root",fTrigName[trigIndx].Data(),fPD.Data(),isMC),"recreate");

  const static long int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nDimu;
  float vz;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu];
  float phi[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float phi1[nMaxDimu];
  float phi2[nMaxDimu];
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float weight0[nMaxDimu]; 
  float weight1[nMaxDimu]; 
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu];
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
  float qxmupl[nMaxDimu];
  float qxmumi[nMaxDimu];
  float qymupl[nMaxDimu];
  float qymumi[nMaxDimu];
  int recoQQsign[nMaxDimu];
  float ctau3D[nMaxDimu];
  float ctau3DErr[nMaxDimu];
  float ctau3DRes[nMaxDimu];
  double TnPweight[nMaxDimu] = {1.};
  double weight = 1;

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
  mmevttree->Branch("runN",&runN,"runN/I");
  mmevttree->Branch("lumi",&lumi,"lumi/I");
  //mmevttree->Branch("cBin",&cBin,"cBin/I");
  mmevttree->Branch("vz",&vz,"vz/F");
  mmevttree->Branch("nDimu",&nDimu,"nDimu/I");
  mmevttree->Branch("mass",mass,"mass[nDimu]/F");
  mmevttree->Branch("y",y,"y[nDimu]/F");
  mmevttree->Branch("pt",pt,"pt[nDimu]/F");
  mmevttree->Branch("pt1",pt1,"pt1[nDimu]/F");
  mmevttree->Branch("pt2",pt2,"pt2[nDimu]/F");
  mmevttree->Branch("eta",eta,"eta[nDimu]/F");
  mmevttree->Branch("eta1",eta1,"eta1[nDimu]/F");
  mmevttree->Branch("eta2",eta2,"eta2[nDimu]/F");
  //mmevttree->Branch("qxa",qxa,"qxa[nDimu]/F");
  //mmevttree->Branch("qxb",qxb,"qxb[nDimu]/F");
  //mmevttree->Branch("qxc",qxc,"qxc[nDimu]/F");
  //mmevttree->Branch("qya",qya,"qya[nDimu]/F");
  //mmevttree->Branch("qyb",qyb,"qyb[nDimu]/F");
  //mmevttree->Branch("qyc",qyc,"qyc[nDimu]/F");
  //mmevttree->Branch("qxdimu",qxdimu,"qxdimu[nDimu]/F");
  //mmevttree->Branch("qydimu",qydimu,"qydimu[nDimu]/F");
  //mmevttree->Branch("qxmupl",qxmupl,"qxmupl[nDimu]/F");
  //mmevttree->Branch("qxmumi",qxmumi,"qxmumi[nDimu]/F");
  //mmevttree->Branch("qymupl",qymupl,"qymupl[nDimu]/F");
  //mmevttree->Branch("qymumi",qymumi,"qymumi[nDimu]/F");
  mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nDimu]/I");
  mmevttree->Branch("ctau3D",ctau3D,"ctau3D[nDimu]/F");
  mmevttree->Branch("ctau3DErr",ctau3DErr,"ctau3DErr[nDimu]/F");
  mmevttree->Branch("ctau3DRes",ctau3DRes,"ctau3DRes[nDimu]/F");
  mmevttree->Branch("weight",&weight,"weight/D");
  mmevttree->Branch("TnPweight",TnPweight,"TnPweight[nDimu]/D");
  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;
//	cout << "Event : " << iev << endl;

    mytree->GetEntry(iev);
    //eptree->GetEntry(iev);
  
    nDimu = 0;
    
	//if( runNb >=327123 ) continue;
    //if(( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) {
	//	cout << "HLTriggers : " << HLTriggers << " kTrigSel : " << kTrigSel << endl;
	//	cout << "HLT PASS" << endl; }

    //if(!( (HLTriggers&((int)pow(2, kTrigSel))) == ((int)pow(2, kTrigSel)) ) ) {
		//cout << "HLT Error!!" << endl;
		//cout << "HLTriggers : " << HLTriggers << " kTrigSel : " << kTrigSel << endl;
		//break;
	//	continue;
	//}
	
	//if (Reco_QQ_size<0) continue;
	//cout << "Reco_QQ_size : " << Reco_QQ_size << endl;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
		//cout << "irqq : " << irqq << endl;
      runN = runNb;
      evt = eventNb;
      lumi = LS;
      //cBin = -999;
      //if(hiHFBinEdge ==0) cBin = getHiBinFromhiHF(SumET_HF);
      //else if(hiHFBinEdge == 1) cBin = getHiBinFromhiHF_Up(SumET_HF);
      //else if(hiHFBinEdge == -1) cBin = getHiBinFromhiHF_Down(SumET_HF);
      //if(cBin==-999){ cout << "ERROR!!! No HF Centrality Matching!!" << endl; return;}
      vz = zVtx;


      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
	  //cout << "HERE1" << endl;
	  //cout << "QQ size : " << Reco_QQ_size << endl;
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
	  //cout << "HERE2" << endl;
	  //cout << "QQ size : " << Reco_QQ_size << endl;
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
	  //cout << "HERE3" << endl;
	  //cout << "QQ size : " << Reco_QQ_size << endl;
      
      weight = 1.;
      //if(isMC) weight = findNcoll(Centrality) * Gen_weight;

      //if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) {
	//	  cout << "Trigger ERROR!!" << endl;
	//	  cout << "Reco_QQ_trig : " << Reco_QQ_trig[irqq] << " kTrigSel : " << kTrigSel << endl;
	//	  continue;
	  //}
      //if(( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) { cout << "Trigger PASS" << endl; }
     
      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      
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
          passMuonTypePl        //                       &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //                        &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      if ( !(muplSoft && mumiSoft) ) { 
	  //cout << "Muon Cut ERORR!!" << endl;
        continue;   
	  }
      //if ( (muplSoft && mumiSoft) ) { cout << "Muon Cut PASS" << endl; }
      
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) {
		  //cout << "VtxProb ERROR!!" << endl;
        continue;
	  }
      //if ( Reco_QQ_VtxProb[irqq] >  0.01 ) { cout << "VtxProb PASS" << endl; }
   
      recoQQsign[irqq] = Reco_QQ_sign[irqq];     
 
      count++;     
      if(isMC){
       tnp_weight = 1;
       tnp_trig_weight_mupl = -1;
       tnp_trig_weight_mumi = -1;
       tnp_weight = tnp_weight * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mupl_Reco->Pt(), mupl_Reco->Eta())) * std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mumi_Reco->Pt(), mumi_Reco->Eta())); //mu id
//       tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0); //inner tracker
       counttnp++;
      }


      // Fill the output tree
      //if ( JP_Reco->Eta() < 0 )  {  
      //  qxa[nDimu] = qx[HFp2];
      //  qya[nDimu] = qy[HFp2];
      //  qxb[nDimu] = qx[HFm2];
      //  qyb[nDimu] = qy[HFm2];

      //}
      //else {
      //  qxa[nDimu] = qx[HFm2];
      //  qya[nDimu] = qy[HFm2];
      //  qxb[nDimu] = qx[HFp2];
      //  qyb[nDimu] = qy[HFp2];
      //}
      //
      //qxc[nDimu] = qx[trackmid2];
      //qyc[nDimu] = qy[trackmid2];

      if(isMC) TnPweight[nDimu] = tnp_weight;
      mass[nDimu] = JP_Reco->M();
      phi[nDimu] = JP_Reco->Phi();
      phi1[nDimu] = mupl_Reco->Phi();
      phi2[nDimu] = mumi_Reco->Phi();
      eta[nDimu] = JP_Reco->Eta();
      y[nDimu] = JP_Reco->Rapidity();
      pt[nDimu] = JP_Reco->Pt();
      pt1[nDimu] = mupl_Reco->Pt();
      pt2[nDimu] = mumi_Reco->Pt();
      eta1[nDimu] = mupl_Reco->Eta();
      eta2[nDimu] = mumi_Reco->Eta();
      //qxdimu[nDimu] = TMath::Cos(2*phi[nDimu]);
      //qydimu[nDimu] = TMath::Sin(2*phi[nDimu]);
      //qxmupl[nDimu] = TMath::Cos(2*phi1[nDimu]);
      //qxmumi[nDimu] = TMath::Cos(2*phi2[nDimu]);
      //qymupl[nDimu] = TMath::Sin(2*phi1[nDimu]);
      //qymumi[nDimu] = TMath::Sin(2*phi2[nDimu]);
      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      ctau3DRes[nDimu] = (Reco_QQ_ctau3D[irqq])/(Reco_QQ_ctauErr3D[irqq]);
      nDimu++;

    } // end of dimuon loop
    
    if(nDimu>0) { mmevttree->Fill(); } //cout << "Fill" << endl; }
    
  } //end of event loop
//  mmtree->Write();  // Don't need to call Write() for trees
  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  newfile->cd();
  mmevttree->Write();
  newfile->Close();
} 
