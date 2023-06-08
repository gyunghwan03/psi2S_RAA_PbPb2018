#include <iostream>

#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tnp_weight_lowptPbPb.h"
#include "TreeSetting.h"

using namespace std;

void getEfficiency_psi_pbpb_SYSTNP(
  bool isPR,
  float ptLow = 0.0, float ptHigh = 30.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 181, bool isTnP = true, bool isPtWeight = true
  ) {

  gStyle->SetOptStat(0);
//  int kTrigSel_=0;
//  if(kTrigSel == kTrigUps) kTrigSel_ = 1;
//  else if(kTrigSel == kTrigL1DBOS40100) kTrigSel_ = 2; 

  float massLow ;
  float massHigh;
  	massLow =2.0;
  	massHigh = 4.8;

  double xmin = ptLow;
  double xmax = ptHigh;

  int CentSwitch = 87;

//  string ftrigSel = (isSwitch) ? "SwitchOn" : "SwitchOff";
//  if(isSwitch) ftrigSel += "_UpsAndL1OS40100";
//  else if(!isSwitch) ftrigSel += Form("_%s",fTrigName[kTrigSel_].Data());

  //input files
  TString inputMC_base, inputMC_ref, inputMC, idstr, trkstr, trgstr;

  string flavorTag = (isPR) ? "PR" : "NP";
  inputMC_base = Form("../TnPSkim/OutputSkim_isMC1_%s_fitVar_noFitNom_tnp_", flavorTag.c_str() ); 
  inputMC_ref = Form("../TnPSkim/OutputSkim_isMC1_%s_fitVar_noFitNom_tnp_idNom_trkNom_trgNom.root", flavorTag.c_str()); 
  auto get_input_config = [&] (int cvr, int mod)
  {
    TString _input_config, _idstr, _trkstr, _trgstr;
    std::cout << "Mode str : " <<mode_str[mod-1].c_str() << std::endl;
    if( cvr == 1 ){ _trkstr = "Nom"; _trgstr = "Nom"; _idstr = Form("%s",mode_str[mod-1].c_str()); }
    if( cvr == 2 ){ _idstr = "Nom"; _trgstr = "Nom"; _trkstr = Form("%s",mode_str[mod-1].c_str()); }
    if( cvr == 3 ){ _idstr = "Nom"; _trkstr = "Nom"; _trgstr = Form("%s",mode_str[mod-1].c_str()); }
    _input_config = Form("id%s_trk%s_trg%s",_idstr.Data(), _trkstr.Data(), _trgstr.Data()); 
    return _input_config;
  };

  TString inputMCs[15];
  int counter =0;
  for( int _cv : {1,2,3}){
    for( int _md : {1,2,3,4,5}){
      inputMCs[counter] = inputMC_base + get_input_config(_cv, _md) + ".root"; 
//      std::cout << "Input MC file["<< _cv*_md -1 << "] = " << inputMCs[_cv*_md -1].Data() << std::endl;
      counter++;
    }
  }

  string rapDef = (yHigh == 2.4) ? "Forward_y_211218" : "y0_1p6_211201";
  TFile *fPtW = new TFile(Form("ratioDataMC_AA_Jpsi_DATA_%s.root", rapDef.c_str()),"read");
  TF1* f1 = (TF1*) fPtW->Get("dataMC_Ratio1");

  TChain* mytree[15];
  for (int i = 0 ; i < 15; i ++){
    mytree[i]= new TChain("tree"); 
    std::cout << "Expecting input MC file["<< i << "] = " << inputMCs[i].Data() << std::endl;
    mytree[i]->Add(inputMCs[i].Data());
  }
  TChain* reftree = new TChain("tree");

  reftree->Add(inputMC_ref.Data());

  //SetBranchAddress
  SetTree_SYSREF settree_;
  settree_.TreeSetting(reftree);
  SetTreeSYS systree[15];
  for (int i = 0 ; i < 15; i ++){
    systree[i].TreeSetting(mytree[i]);
  }
  
  TString histName = Form("%s_pt%.1f_%.1f_y%.1f_%.1f_accYes_mass%.1f_%.1f_cent%d_%d_isTnP%d_isPtWeight%d_%s", flavorTag.c_str(), ptLow,ptHigh,yLow,yHigh,massLow,massHigh,cLow,cHigh,isTnP,isPtWeight, "SYSTNP");
  TH1D *hreco[15], *hreco_tnp[15], *hreco_xtnp[15];
  counter =0;
  for( int _cv : {1,2,3}){
    for( int _md : {1,2,3,4,5}){
      TString typestr = "_" + get_input_config(_cv, _md); 
      hreco[counter] = new TH1D(Form("hreco%s",typestr.Data()),"hreco",1,xmin,xmax);
      hreco_tnp[counter] = new TH1D(Form("hreco_tnp%s"  ,typestr.Data() ),"hreco_tnp",(int) ((xmax-xmin)),xmin,xmax);
      hreco_xtnp[counter] = new TH1D(Form("hreco_xtnp%s",typestr.Data() ),"hreco_xtnp",(int) ((xmax-xmin)),xmin,xmax);
      hreco[counter]->Sumw2();
      hreco_tnp[counter]->Sumw2();
      hreco_xtnp[counter]->Sumw2();

      hreco[counter]->SetTitle("Reco");
      hreco_xtnp[counter]->SetTitle("Reco no tnp");
      hreco_tnp[counter]->SetTitle("Reco with tnp");
      counter++;
    }
  }



// -------NO GEN info in DS--------

  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

//  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
//  double pt_weight = 1;
  
  double tnp_trig_dimu=-1;

  int kL2filter = 38;
  int kL3filter = 39;

  int count =0;
  int counttnp =0;
  const int nevt =reftree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << reftree->GetEntries() <<  " ("<<(int)(100.*iev/reftree->GetEntries()) << "%)" << endl;

    for (int i = 0 ; i < 15; i ++){
      mytree[i]->GetEntry(iev);
    }
    reftree->GetEntry(iev);

    if(!( fabs(y) < yHigh && fabs(y) > yLow && pt < ptHigh && pt >ptLow )) continue;
    if(!( IsAcceptanceQQ(pt1, eta1) && IsAcceptanceQQ(pt2, eta2 )) ) continue;
    if(!( cBin < cHigh && cBin >= cLow)) continue;
    if(!( mass < massHigh && mass > massLow)) continue;
    bool checkID = true; 
    if (checkID) {
      if(!( nTrkWMea1 >5 && nTrkWMea2 >5 && nPixWMea1 > 0 && nPixWMea2 > 0 && fabs(dxy1) < 0.3 && fabs(dxy2) < 0.3 && fabs(dz1) < 20. && fabs(dz2) < 20.) ) continue;

    }
    double ptW =1;
    if( isPtWeight) ptW = f1->Eval(pt);
    for (int i = 0 ; i < 15; i ++){
      hreco[i]->Fill(pt,systree[i].i_weight);
      hreco_tnp[i]->Fill(pt,systree[i].i_weight);
      hreco_xtnp[i]->Fill(pt, (systree[i].i_weight)/(systree[i].i_tnp_weight));
    }

    count++;

  }
  cout << "count " << count << endl;

  TString outFileName = Form("roots/TnP/mc_eff_%s.root",histName.Data());
  TFile* outFile = new TFile(outFileName,"RECREATE");
//  heff->Write();
  for (int i = 0 ; i < 15; i ++){
    hreco[i]->Write();
    hreco_tnp[i]->Write();
    hreco_xtnp[i]->Write();
  }
//  hgen->Write();
  //creco->Close();
  //creco_tnp->Close();
  outFile->Close();

}

