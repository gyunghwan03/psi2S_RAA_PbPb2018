#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "../../../Style.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

void compare_pT_PR(bool isSys=false, double ptHigh=40)
{
  gStyle->SetOptStat(0);
  setTDRStyle();

  TFile *f_mid = new TFile(Form("roots/RAA_psi2S_midRap_pT_%.f.root",ptHigh)); 
  TFile *f_fwd = new TFile(Form("roots/RAA_psi2S_forRap_pT_%.f.root",ptHigh)); 
  TFile *f_mid1S = new TFile("roots/RAA_JPsi_midRap_pT.root"); 
  TFile *f_fwd1S = new TFile("roots/RAA_JPsi_forRap_pT.root"); 
  TFile *fSys1S = new TFile("../../syst_summary_Jpsi/syst_roots/total_syst.root");
  TFile *fSys2S = new TFile("../../syst_summary_psi2S/syst_roots/total_syst.root");
  TFile *f_midJpsi_Old = new TFile("../roots/RAA_PR_Jpsi_HIN_16_025_mid_pT.root");
  TFile *f_fwdJpsi_Old = new TFile("../roots/RAA_PR_Jpsi_HIN_16_025_fwd_pT.root");
  TFile *f_ATLAS_PR = new TFile("../roots/RAA_PR_Jpsi_ATLAS_pT.root");
  TFile *f_ATLAS_NP = new TFile("../roots/RAA_NP_Jpsi_ATLAS_pT.root");

  TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
  TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");
  TH1D *h_midPR1S = (TH1D*) f_mid1S->Get("hRAA_PR");
  TH1D *h_fwdPR1S = (TH1D*) f_fwd1S->Get("hRAA_PR");
  TH1D *hSys1S_mid_PR = (TH1D*) fSys1S->Get("mid_pt_PR");
  TH1D *hSys1S_fwd_PR = (TH1D*) fSys1S->Get("fwd_pt_PR");

  TH1D *hSys2S_mid_PR = (TH1D*) fSys2S->Get("mid_pt_PR");
  TH1D *hSys2S_fwd_PR = (TH1D*) fSys2S->Get("fwd_pt_PR");

  TH1F *h_midPRold_Jpsi =    (TH1F*) f_midJpsi_Old->Get("Table 18/Hist1D_y1");
  TH1F *h_midPRold_JpsiErr = (TH1F*) f_midJpsi_Old->Get("Table 18/Hist1D_y1_e1");
  TH1F *h_midPRold_JpsiSys = (TH1F*) f_midJpsi_Old->Get("Table 18/Hist1D_y1_e2");

  TH1F *h_fwdPRold_Jpsi =    (TH1F*) f_fwdJpsi_Old->Get("Table 22/Hist1D_y1");
  TH1F *h_fwdPRold_JpsiErr = (TH1F*) f_fwdJpsi_Old->Get("Table 22/Hist1D_y1_e1");
  TH1F *h_fwdPRold_JpsiSys = (TH1F*) f_fwdJpsi_Old->Get("Table 22/Hist1D_y1_e2");


  const int nPtBins=6;
  const int nPtBins_fwd=4;
  double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,ptHigh};
  double ptBin_fwd[nPtBins_fwd+1] = {3.5,6.5,9,12,ptHigh};
  double ptBin_old[nPtBins] = {6.5,9,12,15,20,30};
  double ptBin_fwd_old[nPtBins_fwd-1] = {6.5,12,30};
  double ptBin_alice[5] = {0,2,4,6,12};
  double x_alice[4] = {1,3,5,9};
  double binWidth_alice[4] = {1,1,1,3};
  double x[nPtBins]; double binWidth[nPtBins]; 
  double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];
  double x_old[nPtBins-1]; double binWidth_old[nPtBins-1];
  double x_fwd_old[2]={(6.5+12)/2,(12+30)/2}; double binWidth_fwd_old[2]={(12-6.5)/2,(30-12)/2};
  double x_fwd_Jpsi[3]={(3+6.5)/2,(6.5+12)/2,(12+30)/2}; double binWidth_fwd_Jpsi[3]={(6.5-3)/2,(12-6.5)/2,(30-12)/2};

  double midPR_new[nPtBins]; 
  double midPR_new_Err[nPtBins]; 
  double midPR_new_Sys[nPtBins]; 
  double fwdPR_new[nPtBins_fwd]; 
  double fwdPR_new_Err[nPtBins_fwd]; 
  double fwdPR_new_Sys[nPtBins_fwd]; 

  double midPR1S_new[nPtBins]; 
  double midPR1S_new_Err[nPtBins]; 
  double fwdPR1S_new[nPtBins_fwd]; 
  double fwdPR1S_new_Err[nPtBins_fwd]; 

  double SysPR1S_mid[nPtBins]; 
  double SysPR1S_fwd[nPtBins_fwd];
  double SysPR2S_mid[nPtBins]; 
  double SysPR2S_fwd[nPtBins_fwd];

  double PRJpsi_midOld[nPtBins-1];
  double PRJpsi_midOld_Err[nPtBins-1];
  double PRJpsi_midOld_Sys[nPtBins-1];

  double PRJpsi_fwdOld[nPtBins-1];
  double PRJpsi_fwdOld_Err[nPtBins-1];
  double PRJpsi_fwdOld_Sys[nPtBins-1];

  double midPR_old[nPtBins-1] = {0.109,0.128,0.172,0.142,0.123};
  double midPR_old_Err[nPtBins-1] = {0.057,0.034,0.048,0.061,0.109};
  double midPR_old_Sys[nPtBins-1] = {0.018,0.013,0.018,0.020,0.034};
  double fwdPR_old[2] = {0.135,0.247};
  double fwdPR_old_Err[2] = {0.079,0.098}; 
  double fwdPR_old_Sys[2] = {0.038,0.042};

  double alice_var[4]={0.41,0.35,0.22,0.15};
  double alice_err[4]={0.1,0.07,0.07,0.06};
  double alice_sys[4]={0.11,0.08,0.06,0.02};

  for (int i=0; i<nPtBins; i++){
    midPR_new[i] = h_midPR->GetBinContent(i+1);
    midPR_new_Err[i] = h_midPR->GetBinError(i+1);
    midPR1S_new[i] = h_midPR1S->GetBinContent(i+1);
    midPR1S_new_Err[i] = h_midPR1S->GetBinError(i+1);
    x[i] = (ptBin[i+1]+ptBin[i])/2;
    binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
    SysPR1S_mid[i] = midPR1S_new[i]*(hSys1S_mid_PR->GetBinContent(i+1));
    SysPR2S_mid[i] = midPR_new[i]*(hSys2S_mid_PR->GetBinContent(i+1));
  }
  for (int i=0; i<nPtBins-1; i++){
    PRJpsi_midOld[i] = h_midPRold_Jpsi->GetBinContent(i+1);
    PRJpsi_midOld_Err[i] = h_midPRold_JpsiErr->GetBinContent(i+1);
    PRJpsi_midOld_Sys[i] = h_midPRold_JpsiSys->GetBinContent(i+1);
    x_old[i] = (ptBin_old[i+1]+ptBin_old[i])/2;
    binWidth_old[i] = (ptBin_old[i+1]-ptBin_old[i])/2;
  }
  for (int i=0; i<3; i++){
    PRJpsi_fwdOld[i]     = h_fwdPRold_Jpsi->GetBinContent(i+1);
    PRJpsi_fwdOld_Err[i] = h_fwdPRold_JpsiErr->GetBinContent(i+1);
    PRJpsi_fwdOld_Sys[i] = h_fwdPRold_JpsiSys->GetBinContent(i+1);
  }
  for (int i=0; i<nPtBins_fwd; i++){
    fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
    fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
    fwdPR1S_new[i] = h_fwdPR1S->GetBinContent(i+1);
    fwdPR1S_new_Err[i] = h_fwdPR1S->GetBinError(i+1);
    x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
    binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
    SysPR1S_fwd[i] = fwdPR1S_new[i]*(hSys1S_fwd_PR->GetBinContent(i+1));
    SysPR2S_fwd[i] = fwdPR_new[i]*(hSys2S_fwd_PR->GetBinContent(i+1));
  }

  TGraphErrors *g_midPR;
  TGraphErrors *g_midPR1S;
  TGraphErrors *g_midPR_old;
  TGraphErrors *g_midPRold_Jpsi;
  TGraphErrors *g_fwdPRold_Jpsi;

  if(isSys==1){
    g_midPR = new TGraphErrors(nPtBins, x, midPR_new, 0, midPR_new_Err);
    g_midPR1S = new TGraphErrors(nPtBins, x, midPR1S_new, 0, midPR1S_new_Err);
    g_midPR_old = new TGraphErrors(nPtBins - 1, x_old, midPR_old, 0, midPR_old_Err);
    g_midPRold_Jpsi = new TGraphErrors(nPtBins - 1, x_old, PRJpsi_midOld, 0, PRJpsi_midOld_Err);
    g_fwdPRold_Jpsi = new TGraphErrors(3, x_fwd_Jpsi, PRJpsi_fwdOld, 0, PRJpsi_fwdOld_Err);
  }
  else {
    g_midPR = new TGraphErrors(nPtBins, x, midPR_new, binWidth, midPR_new_Err);
    g_midPR1S = new TGraphErrors(nPtBins, x, midPR1S_new, binWidth, midPR1S_new_Err);
    g_midPR_old = new TGraphErrors(nPtBins - 1, x_old, midPR_old, binWidth_old, midPR_old_Err);
    g_midPRold_Jpsi = new TGraphErrors(nPtBins - 1, x_old, PRJpsi_midOld, binWidth_old, PRJpsi_midOld_Err);
    g_fwdPRold_Jpsi = new TGraphErrors(3, x_fwd_Jpsi, PRJpsi_fwdOld, binWidth_fwd_Jpsi, PRJpsi_fwdOld_Err);
  }
  TGraphErrors *g_midPRSys1S = new TGraphErrors(nPtBins,x,midPR1S_new,binWidth,SysPR1S_mid);
  TGraphErrors *g_midPRSys = new TGraphErrors(nPtBins,x,midPR_new,binWidth,SysPR2S_mid);
  TGraphErrors *g_midPR_oldSys = new TGraphErrors(nPtBins-1,x_old,midPR_old,binWidth_old,midPR_old_Sys);
  TGraphErrors *g_midPRold_JpsiSys = new TGraphErrors(nPtBins-1,x_old,PRJpsi_midOld,binWidth_old,PRJpsi_midOld_Sys);
  TGraphErrors *g_fwdPRold_JpsiSys = new TGraphErrors(3,x_fwd_Jpsi,PRJpsi_fwdOld,binWidth_fwd_Jpsi,PRJpsi_fwdOld_Sys);

  TGraphErrors *g_fwdPR;
  TGraphErrors *g_fwdPR1S;
  TGraphErrors *g_fwdPR_old;
  TGraphErrors *g_alice;

  if(isSys==1){
    g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR_new, 0, fwdPR_new_Err);
    g_fwdPR1S = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR1S_new, 0, fwdPR1S_new_Err);
    g_fwdPR_old = new TGraphErrors(2, x_fwd_old, fwdPR_old, 0, fwdPR_old_Err);
    g_alice = new TGraphErrors(4, x_alice, alice_var, 0, alice_err);
  }
  else{
    g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR_new, binWidth_fwd, fwdPR_new_Err);
    g_fwdPR1S = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR1S_new, binWidth_fwd, fwdPR1S_new_Err);
    g_fwdPR_old = new TGraphErrors(2, x_fwd_old, fwdPR_old, binWidth_fwd_old, fwdPR_old_Err);
    g_alice = new TGraphErrors(4, x_alice, alice_var, binWidth_alice, alice_err);
  }
  TGraphErrors *g_fwdPRSys1S = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR1S_new,binWidth_fwd,SysPR1S_fwd);
  TGraphErrors *g_fwdPRSys = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,binWidth_fwd,SysPR2S_fwd);
  TGraphErrors *g_fwdPR_oldSys = new TGraphErrors(2,x_fwd_old,fwdPR_old,binWidth_fwd_old,fwdPR_old_Sys);
  TGraphErrors *g_aliceSys = new TGraphErrors(4,x_alice,alice_var,binWidth_alice,alice_sys);

  TGraphErrors *g_ATLAS_PR = (TGraphErrors*)f_ATLAS_PR->Get("Table 5/Graph1D_y1");
  TGraphAsymmErrors *g_ATLAS_PR_Sys = (TGraphAsymmErrors*)f_ATLAS_PR->Get("Table 5/Graph1D_y1"); 

  TH1F *h_ATLAS_PR_Syspl = (TH1F*) f_ATLAS_PR->Get("Table 5/Hist1D_y1_e2plus");
  TH1F *h_ATLAS_PR_Sysmi = (TH1F*) f_ATLAS_PR->Get("Table 5/Hist1D_y1_e2minus");

  for(int i=0; i<g_ATLAS_PR->GetN(); i++){
    double PR_Syspl, PR_Sysmi, NP_Syspl, NP_Sysmi;
    double x_PR, x_NP, y_PR, y_NP;

    g_ATLAS_PR->GetPoint(i,x_PR,y_PR);

    PR_Syspl = h_ATLAS_PR_Syspl->GetBinContent(i*2+1);
    PR_Sysmi = h_ATLAS_PR_Sysmi->GetBinContent(i*2+1);

    g_ATLAS_PR_Sys->SetPointEYhigh(i,fabs(PR_Syspl));
    g_ATLAS_PR_Sys->SetPointEYlow(i,fabs(PR_Sysmi));
    }

  cout <<"HERE" << endl;
  double x_1st = (6.5+3)/2;
  TArrow *fwdPR_upper = new TArrow(x_1st,0,x_1st,0.395,0.027,"<-|");

  int iPeriod = 101;
  int iPos = 33;

  double text_x = 0.18;
  double text_y = 0.8;
  double y_diff = 0.08;
  float pos_x = 0.21;
  float pos_y = 0.87;
  float pos_y_diff = 0.051;
  int text_color = 1;
  float text_size = 25;

  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();

  g_midPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPR->GetXaxis()->CenterTitle();
  g_midPR->GetYaxis()->SetTitle("R_{AA}");
  g_midPR->GetYaxis()->CenterTitle();
  g_midPR->SetTitle();
  g_midPR->GetXaxis()->SetLimits(0.,ptHigh);
  g_midPR->SetMinimum(0.);
  g_midPR->SetMaximum(1.44);

  g_midPR->SetMarkerColor(kBlue+2);
  g_midPR->SetLineColor(kBlue+2);
  g_midPR->SetMarkerStyle(20);
  g_midPR->SetMarkerSize(1.4);
  g_midPRSys->SetLineColor(kBlue-4);
  g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);

  g_midPR1S->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPR1S->GetXaxis()->CenterTitle();
  g_midPR1S->GetYaxis()->SetTitle("R_{AA}");
  g_midPR1S->GetYaxis()->CenterTitle();
  g_midPR1S->SetTitle();
  g_midPR1S->GetXaxis()->SetLimits(0.,ptHigh);
  g_midPR1S->SetMinimum(0.);
  g_midPR1S->SetMaximum(1.44);
  g_midPR1S->SetMarkerStyle(33);
  g_midPR1S->SetMarkerSize(1.9);
  g_midPR1S->SetMarkerColor(kViolet+2);
  g_midPR1S->SetLineColor(kViolet+2);
  g_midPRSys1S->SetLineColor(kViolet-4);
  g_midPRSys1S->SetFillColorAlpha(kViolet-9,0.40);


  g_midPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPR_old->GetXaxis()->CenterTitle();
  g_midPR_old->GetYaxis()->SetTitle("R_{AA}");
  g_midPR_old->GetYaxis()->CenterTitle();
  g_midPR_old->SetTitle();
  g_midPR_old->GetXaxis()->SetLimits(0.,ptHigh);
  g_midPR_old->SetMinimum(0.);
  g_midPR_old->SetMaximum(1.44);
  g_midPR_old->SetMarkerStyle(25);
  g_midPR_old->SetMarkerSize(1.4);
  g_midPR_old->SetMarkerColor(15);
  g_midPR_old->SetLineColor(15);

  g_midPR_oldSys->SetLineColor(15);
  g_midPR_oldSys->SetFillColorAlpha(15,0.4);

  g_midPRold_Jpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPRold_Jpsi->GetXaxis()->CenterTitle();
  g_midPRold_Jpsi->GetYaxis()->SetTitle("R_{AA}");
  g_midPRold_Jpsi->GetYaxis()->CenterTitle();
  g_midPRold_Jpsi->SetTitle();
  g_midPRold_Jpsi->GetXaxis()->SetLimits(0.,ptHigh);
  g_midPRold_Jpsi->SetMinimum(0.);
  g_midPRold_Jpsi->SetMaximum(1.44);
  g_midPRold_Jpsi->SetMarkerStyle(27);
  g_midPRold_Jpsi->SetMarkerSize(1.4);
  g_midPRold_Jpsi->SetMarkerColor(40);
  g_midPRold_Jpsi->SetLineColor(40);

  g_ATLAS_PR->SetMarkerStyle(27);
  g_ATLAS_PR->SetMarkerSize(1.4);
  g_ATLAS_PR->SetMarkerColor(45);
  g_ATLAS_PR->SetLineColor(45);
  g_ATLAS_PR_Sys->SetLineColor(45);
  g_ATLAS_PR_Sys->SetFillColorAlpha(45,0.4);

  g_midPRold_JpsiSys->SetLineColor(40);
  g_midPRold_JpsiSys->SetFillColorAlpha(40,0.4);


  g_midPR->Draw("AP");
  g_midPR1S->Draw("P");
  g_midPR_old->Draw("P");
  g_midPRold_Jpsi->Draw("P");
  g_ATLAS_PR->Draw("P");
  if(isSys==1){
	g_midPRSys1S->Draw("5");
    g_midPRSys->Draw("5");
    g_midPR_oldSys->Draw("5");
    g_midPRold_JpsiSys->Draw("5");
    g_ATLAS_PR_Sys->Draw("5");
  }

  TLegend *leg1 = new TLegend(0.38,0.69,0.65,0.54);
  leg1->SetTextSize(text_size);
  leg1->SetTextFont(43);
  leg1->SetBorderSize(0);
  leg1->AddEntry(g_midPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%", "pe");
  leg1->AddEntry(g_midPR1S,"#bf{Prompt J/#psi} New Data, Cent. 0-90%", "pe");
  leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, Cent. 0-100%", "pe");
  leg1->AddEntry(g_midPRold_Jpsi,"#bf{Prompt J/#psi} HIN-16-025, Cent. 0-100%", "pe");
  leg1->AddEntry(g_ATLAS_PR,"#bf{Prompt J/#psi} ATLAS, |y| < 2, Cent. 0-80%", "pe");
  leg1->Draw("SAME");
  jumSun(0,1,ptHigh,1);

  drawText(Form("6.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c1, iPeriod, iPos);

  c1->SaveAs(Form("./figs/compare_PR_mid_pT_%.f_Sys%d.pdf",ptHigh,isSys));

  TCanvas *c2 = new TCanvas("c2","",900,800);
  c2->cd();
  g_fwdPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_fwdPR->GetXaxis()->CenterTitle();
  g_fwdPR->GetYaxis()->SetTitle("R_{AA}");
  g_fwdPR->GetYaxis()->CenterTitle();
  g_fwdPR->SetTitle();
  g_fwdPR->GetXaxis()->SetLimits(0.,ptHigh);
  g_fwdPR->SetMinimum(0.);
  g_fwdPR->SetMaximum(1.44);

  g_fwdPR->SetMarkerColor(kBlue+2);
  g_fwdPR->SetLineColor(kBlue+2);
  g_fwdPR->SetMarkerStyle(20);
  g_fwdPR->SetMarkerSize(1.4);
  g_fwdPRSys->SetLineColor(kBlue-4);
  g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

  g_fwdPR1S->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_fwdPR1S->GetXaxis()->CenterTitle();
  g_fwdPR1S->GetYaxis()->SetTitle("R_{AA}");
  g_fwdPR1S->GetYaxis()->CenterTitle();
  g_fwdPR1S->SetTitle();
  g_fwdPR1S->GetXaxis()->SetLimits(0.,ptHigh);
  g_fwdPR1S->SetMinimum(0.);
  g_fwdPR1S->SetMaximum(1.44);
  g_fwdPR1S->SetMarkerStyle(33);
  g_fwdPR1S->SetMarkerSize(2.0);
  g_fwdPR1S->SetMarkerColor(kViolet+2);
  g_fwdPR1S->SetLineColor(kViolet+2);
  g_fwdPRSys1S->SetLineColor(kViolet-4);
  g_fwdPRSys1S->SetFillColorAlpha(kViolet-9,0.40);


  g_fwdPR_old->SetMarkerStyle(25);
  g_fwdPR_old->SetMarkerSize(1.4);
  g_fwdPR_old->SetMarkerColor(15);
  g_fwdPR_old->SetLineColor(15);

  g_alice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_alice->GetXaxis()->CenterTitle();
  g_alice->GetYaxis()->SetTitle("R_{AA}");
  g_alice->GetYaxis()->CenterTitle();
  g_alice->SetTitle();
  g_alice->GetXaxis()->SetLimits(0.,ptHigh);
  g_alice->SetMinimum(0.);
  g_alice->SetMaximum(1.44);
  g_alice->SetMarkerStyle(33);
  g_alice->SetMarkerSize(1.8);
  g_alice->SetMarkerColor(45);
  g_alice->SetLineColor(45);
  
  g_aliceSys->SetLineColor(42);
  g_aliceSys->SetFillColorAlpha(42,0.4);

  g_fwdPR_oldSys->SetLineColor(15);
  g_fwdPR_oldSys->SetFillColorAlpha(15,0.4);

  g_fwdPRold_Jpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_fwdPRold_Jpsi->GetXaxis()->CenterTitle();
  g_fwdPRold_Jpsi->GetYaxis()->SetTitle("R_{AA}");
  g_fwdPRold_Jpsi->GetYaxis()->CenterTitle();
  g_fwdPRold_Jpsi->SetTitle();
  g_fwdPRold_Jpsi->GetXaxis()->SetLimits(0.,ptHigh);
  g_fwdPRold_Jpsi->SetMinimum(0.);
  g_fwdPRold_Jpsi->SetMaximum(1.44);
  g_fwdPRold_Jpsi->SetMarkerStyle(27);
  g_fwdPRold_Jpsi->SetMarkerSize(1.4);
  g_fwdPRold_Jpsi->SetMarkerColor(40);
  g_fwdPRold_Jpsi->SetLineColor(40);

  g_fwdPRold_JpsiSys->SetLineColor(40);
  g_fwdPRold_JpsiSys->SetFillColorAlpha(40,0.4);

  fwdPR_upper->SetLineColor(15);
  fwdPR_upper->SetLineWidth(2);

  g_fwdPR->Draw("AP");
  g_fwdPR1S->Draw("P");
  g_fwdPR_old->Draw("P");
  g_fwdPRold_Jpsi->Draw("P");
  if(isSys==1){
    g_fwdPRSys1S->Draw("5");
    g_fwdPRSys->Draw("5");
    g_fwdPR_oldSys->Draw("5");
    g_aliceSys->Draw("5");
    g_fwdPRold_JpsiSys->Draw("5");
  }
  fwdPR_upper->Draw("");
  g_alice->Draw("P");

  TLegend *leg2 = new TLegend(0.29,0.69,0.63,0.5);
  leg2->SetTextSize(text_size);
  leg2->SetTextFont(43);
  leg2->SetBorderSize(0);
  leg2->AddEntry(g_fwdPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%", "pe");
  leg2->AddEntry(g_fwdPR1S,"#bf{Prompt J/#psi} New Data, Cent. 0-90%", "pe");
  leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, Cent. 0-100%", "pe");
  leg2->AddEntry(g_fwdPRold_Jpsi,"#bf{Prompt J/#psi} HIN-16-025, Cent. 0-100%", "pe");
  leg2->AddEntry(g_alice,"#bf{Inclusive #psi(2S)} ALICE, 2.5 < y < 4, Cent.0-90%", "pe");
  leg2->Draw("SAME");
  jumSun(0,1,ptHigh,1);

  drawText(Form("3.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c2, iPeriod, iPos);

  c2->SaveAs(Form("./figs/compare_PR_fwd_pT_%.f_Sys%d.pdf",ptHigh,isSys));

  TCanvas *c3 = new TCanvas("c3","",900,800);
  c3->cd();
  g_midPR->Draw("AP");
  g_midPR_old->Draw("P");
  if(isSys==1){
    g_midPRSys->Draw("5");
    g_midPR_oldSys->Draw("5");
  }

  TLegend *leg3 = new TLegend(0.38,0.69,0.65,0.56);
  leg3->SetTextSize(text_size);
  leg3->SetTextFont(43);
  leg3->SetBorderSize(0);
  leg3->AddEntry(g_midPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%", "pe");
  leg3->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, Cent. 0-100%", "pe");
  leg3->Draw("SAME");
  jumSun(0,1,ptHigh,1);

  drawText(Form("6.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
  drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c3, iPeriod, iPos);

  c3->SaveAs(Form("./figs/compare_PR_mid_pT_%.f_Sys%d_psi2S.pdf",ptHigh,isSys));

  TCanvas *c4 = new TCanvas("c4","",900,800);
  c4->cd();
  g_fwdPR->Draw("AP");
  g_fwdPR_old->Draw("P");
  g_alice->Draw("P");
  if(isSys==1){
    g_fwdPRSys->Draw("5");
    g_fwdPR_oldSys->Draw("5"); 
    g_aliceSys->Draw("5");
  }
  fwdPR_upper->Draw("");

  TLegend *leg4 = new TLegend(0.29,0.69,0.63,0.56);
  leg4->SetTextSize(text_size);
  leg4->SetTextFont(43);
  leg4->SetBorderSize(0);
  leg4->AddEntry(g_fwdPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%", "pe");
  leg4->AddEntry(g_fwdPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, Cent. 0-100%", "pe");
  leg4->AddEntry(g_alice,"#bf{Inclusive #psi(2S)} ALICE, 2.5 < y < 4, Cent. 0-90%", "pe");
  leg4->Draw("SAME");
  jumSun(0,1,ptHigh,1);

  drawText(Form("3.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
  drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c4, iPeriod, iPos);

  c4->SaveAs(Form("./figs/compare_PR_fwd_pT_%.f_Sys%d_psi2S.pdf",ptHigh,isSys));

}
