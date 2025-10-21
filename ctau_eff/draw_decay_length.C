#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_v2mass.C"
#include "../rootFitHeaders.h"

using namespace std;

void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void draw_decay_length(double ptLow =  3, double ptHigh = 30, 
    double yLow = 1.6, double yHigh = 2.4,
    double SiMuPtCut = 0, double massLow = 2.6, double massHigh = 3.5, int PRMC=0, bool dimusign=true, bool fAccW = false, bool fEffW = false)
{
  // Jpsi mass range: 2.6 - 3.5
  // 2S mass range: 3.3 - 4.1

  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, SiMuPtCut) ;
  kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh);

  //Basic Setting
  gStyle->SetOptStat(0);
  //  setTDRStyle();

  TFile *fData;
  TFile *fPRMC;
  TFile *fNPMC;
  TTree *treeData;
  TTree *treePRMC;
  TTree *treeNPMC;

  fData = new TFile("../skimmedFiles/OniaFlowSkim_JpsiTrig_DoubleMuonPD_pp_isMC0_221226.root");
  fPRMC = new TFile("../skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root");
  fNPMC = new TFile("../skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_241011.root");

  TCut mCut;
  TCut ptCut;
  TCut yCut;
  TCut cCut;
  TCut sglMuCut;
  TCut accCut;
  TCut dimuSign;
  TCut quality;
  quality = ("fabs(vz)<15");
  dimuSign = ("recoQQsign==0");
  mCut = Form("mass>=%.1f&&mass<=%.1f",massLow, massHigh);
  ptCut = Form("pt>%.1f&&pt<%.1f", ptLow, ptHigh);
  yCut  = Form("abs(y)>%.1f&&abs(y)<%.1f", yLow, yHigh);
  //cCut  = Form("cBin>=%d&&cBin<=%d",cLow,cHigh);
  sglMuCut = Form("abs(eta1)>%.1f&&abs(eta1)<%.1f&&abs(eta2)>%.1f&&abs(eta2)<%.1f", yLow, yHigh, yLow, yHigh);
  accCut= ("( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || \
      ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || \
      ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )");

  cout<<kineLabel<<endl;
  TCut totalCut ="";
  totalCut = quality&&dimuSign&&mCut&&ptCut&&yCut&&sglMuCut&&accCut;
  //totalCut = quality&&dimuSign&&mCut&&ptCut&&yCut&&cCut&&sglMuCut;
  cout<<totalCut<<endl;
  //if(forPR==0)consider="PR_eff";
  //else if(forPR==1)consider="NP_eff";

  float nbins=4000;
  float xmin =  -1;
  float xmax =   3;
  TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
  TH1D* h_massCut = new TH1D("h_massCut",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
  TH1D* h_decayData = new TH1D("h_decayData",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_decayPRMC = new TH1D("h_decayPRMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_decayNPMC = new TH1D("h_decayNPMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_deffData = new TH1D("h_deffData",";l_{J/#psi};Efficiency",nbins,xmin,xmax);
  TH1D* h_deffPRMC = new TH1D("h_deffPRMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);
  TH1D* h_deffNPMC = new TH1D("h_deffNPMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);

  treeData = (TTree*)fData->Get("mmepevt");
  treePRMC = (TTree*)fPRMC->Get("mmepevt");
  treeNPMC = (TTree*)fNPMC->Get("mmepevt");

  treeData->Draw("ctau3D>>h_decayData",totalCut,"l");
  treePRMC->Draw("ctau3D>>h_decayPRMC",totalCut,"l");
  treeNPMC->Draw("ctau3D>>h_decayNPMC",totalCut,"l");

  cout<<"Data Tree: "<<treeData->GetEntries()<<", by total cut: "<<h_decayData->GetEntries()<<", "<<h_decayData->GetEntries()/treeData->GetEntries()*100<<"%"<<endl;
  cout<<"PRMC Tree: "<<treePRMC->GetEntries()<<", by total cut: "<<h_decayPRMC->GetEntries()<<", "<<h_decayData->GetEntries()/treePRMC->GetEntries()*100<<"%"<<endl;
  cout<<"NPMC Tree: "<<treeNPMC->GetEntries()<<", by total cut: "<<h_decayNPMC->GetEntries()<<", "<<h_decayData->GetEntries()/treeNPMC->GetEntries()*100<<"%"<<endl;

  // 단순히 컷 없이 히스토그램만 그림
  int totalPRMC = (int)h_decayPRMC->GetEntries();
  int totalNPMC = (int)h_decayNPMC->GetEntries();

  // mass 히스토그램: 컷 없이 동일 totalCut 만 사용 (ctau 컷 제거)
  treeData->Draw("mass>>h_mass", totalCut, "EP");
  treeData->Draw("mass>>h_massCut", totalCut, "EP"); // h_massCut: now same selection as h_mass (no ctau cut)
  cout<<"How many Jpsi (mass histogram) : "<<h_mass->GetEntries()<<endl;

  // 누적 효율 곡선 계산 (컷 기준 없이 PR / NP 누적 분포만 그림)
  int nBins = h_decayPRMC->GetNbinsX();
  if(totalPRMC>0 && totalNPMC>0){
    for(int ib=1; ib<=nBins; ++ib){
      double cumPR = h_decayPRMC->Integral(1, ib);
      double cumNP = h_decayNPMC->Integral(1, ib);
      h_deffPRMC->SetBinContent(ib, cumPR / (double)totalPRMC);
      h_deffNPMC->SetBinContent(ib, 1.0 - cumNP / (double)totalNPMC);
    }
  } else {
    cout<<"Warning: PRMC or NPMC empty, efficiency histograms will be empty."<<endl;
  }

  TCanvas* c_mass = new TCanvas("c_mass","",600,600);
  h_mass->Draw();
  h_mass->SetMinimum(0);
  h_mass->SetLineColor(kBlack);
  h_massCut->Draw("same");
  h_massCut->SetLineColor(kRed);
  c_mass->Update();
  if(PRMC==0){
  c_mass->SaveAs(Form("PR_pp_mass_%s.pdf", kineLabel.Data()));}
  else if(PRMC==1){
  c_mass->SaveAs(Form("NP_pp_mass_%s.pdf", kineLabel.Data()));}

  float pos_x = 0.23;
  float pos_x_mass = 0.55;
  float pos_y = 0.65;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 16;
  TString perc = "%";

  TCanvas* c_decayL = new TCanvas("c_decayL","",600,600);
  h_deffPRMC->Draw("l");
  h_deffNPMC->Draw("same");
  h_deffNPMC->SetLineColor(kRed+2);
  jumSun(xmin,1,xmax,1,1,1);

  if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff*0.5,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*1.5,text_color,text_size);

  // 저장: PRMC / NPMC 별 폴더에 결과 저장
  c_decayL->Update();
  if(PRMC==0)  c_decayL->SaveAs(Form("PR_pp_decay_%s.pdf",kineLabel.Data()));
  else         c_decayL->SaveAs(Form("NP_pp_decay_%s.pdf",kineLabel.Data()));
  TFile *wf = new TFile(Form("roots_1S_pp/decayL/%s/decay_hist_%s.root", (PRMC==0?"PRMC":"NPMC"), kineLabel.Data()),"recreate");
  wf->cd();

  h_deffPRMC->Write();
  h_deffNPMC->Write();
  h_mass->Write();
  h_massCut->Write();

}

void GetHistSqrt(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
  double content;
  double err;
  for(int i=1; i<=h1->GetNbinsX(); i++){
    content=0;err=0;
    content = h1->GetBinContent(i);
    err = h1->GetBinError(i);
    err = 0.5*err*TMath::Power(content,-0.5);
    h2->SetBinContent(i,TMath::Sqrt(content));
    h2->SetBinError(i,err);
  }
}
