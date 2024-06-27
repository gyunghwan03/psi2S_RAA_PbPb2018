#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"
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

void compare_pT_PR(bool isSys=false)
{
  gStyle->SetOptStat(0);
  setTDRStyle();

  TFile *f_mid = new TFile("roots/RAA_psi2S_midRap_pT.root"); 
  TFile *f_fwd = new TFile("roots/RAA_psi2S_forRap_pT.root"); 
  TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst.root");
  TFile *fOld = new TFile("roots/RAA_PR_Jpsi_HIN_16_025_pT.root");

  TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
  TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

  TH1D *hSys_PR = (TH1D*) fSys->Get("mid_pt_PR");
  TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_pt_PR");

  TH1F *h_PR_Jpsi = (TH1F*) fOld->Get("Table 18/Hist1D_y1");
  TH1F *h_PR_JpsiErr = (TH1F*) fOld->Get("Table 18/Hist1D_y1_e1");
  TH1F *h_PR_JpsiSys = (TH1F*) fOld->Get("Table 18/Hist1D_y1_e2");


  const int nPtBins=6;
  const int nPtBins_fwd=4;
  double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,50};
  double ptBin_fwd[nPtBins_fwd+1] = {3.5,6.5,9,12,50};
  double ptBin_old[nPtBins] = {6.5,9,12,15,20,30};
  double ptBin_fwd_old[nPtBins_fwd-1] = {6.5,12,30};
  double ptBin_alice[5] = {0,2,4,6,12};
  double x_alice[4] = {1,3,5,9};
  double binWidth_alice[4] = {1,1,1,3};
  double x[nPtBins]; double binWidth[nPtBins]; 
  double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];
  double x_old[nPtBins-1]; double binWidth_old[nPtBins-1];
  double x_fwd_old[2]={(6.5+12)/2,(12+30)/2}; double binWidth_fwd_old[2]={(12-6.5)/2,(30-12)/2};

  double midPR_new[nPtBins]; 
  double midPR_new_Err[nPtBins]; 
  double midPR_new_Sys[nPtBins]; 
  double fwdPR_new[nPtBins_fwd]; 
  double fwdPR_new_Err[nPtBins_fwd]; 
  double fwdPR_new_Sys[nPtBins_fwd]; 

  double SysPR[nPtBins]; 
  double SysPR_fwd[nPtBins_fwd];

  double PRJpsi[nPtBins-1];
  double PRJpsi_Err[nPtBins-1];
  double PRJpsi_Sys[nPtBins-1];

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
    x[i] = (ptBin[i+1]+ptBin[i])/2;
    binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
    SysPR[i] = midPR_new[i]*(hSys_PR->GetBinContent(i+1));
    PRJpsi[i] = h_PR_Jpsi->GetBinContent(i+1);
  }
  for (int i=0; i<nPtBins-1; i++){
    PRJpsi_Err[i] = h_PR_JpsiErr->GetBinContent(i+1);
    PRJpsi_Sys[i] = h_PR_JpsiSys->GetBinContent(i+1);
    x_old[i] = (ptBin_old[i+1]+ptBin_old[i])/2;
    binWidth_old[i] = (ptBin_old[i+1]-ptBin_old[i])/2;
  }
  for (int i=0; i<nPtBins_fwd; i++){
    fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
    fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
    x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
    binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
    SysPR_fwd[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
  }

  TGraphErrors *g_midPR;
  TGraphErrors *g_midPR_old;
  TGraphErrors *g_midPR_Jpsi;

  if(isSys==1){
    g_midPR = new TGraphErrors(nPtBins, x, midPR_new, 0, midPR_new_Err);
    g_midPR_old = new TGraphErrors(nPtBins - 1, x_old, midPR_old, 0, midPR_old_Err);
    g_midPR_Jpsi = new TGraphErrors(nPtBins - 1, x_old, PRJpsi, 0, PRJpsi_Err);
  }
  else {
    g_midPR = new TGraphErrors(nPtBins, x, midPR_new, binWidth, midPR_new_Err);
    g_midPR_old = new TGraphErrors(nPtBins - 1, x_old, midPR_old, binWidth_old, midPR_old_Err);
    g_midPR_Jpsi = new TGraphErrors(nPtBins - 1, x_old, PRJpsi, binWidth_old, PRJpsi_Err);
  }
  TGraphErrors *g_midPRSys = new TGraphErrors(nPtBins,x,midPR_new,binWidth,SysPR);
  TGraphErrors *g_midPR_oldSys = new TGraphErrors(nPtBins-1,x_old,midPR_old,binWidth_old,midPR_old_Sys);
  TGraphErrors *g_midPR_JpsiSys = new TGraphErrors(nPtBins-1,x_old,PRJpsi,binWidth_old,PRJpsi_Sys);

  TGraphErrors *g_fwdPR;
  TGraphErrors *g_fwdPR_old;
  TGraphErrors *g_alice;

  if(isSys==1){
    g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR_new, 0, fwdPR_new_Err);
    g_fwdPR_old = new TGraphErrors(2, x_fwd_old, fwdPR_old, 0, fwdPR_old_Err);
	g_alice = new TGraphErrors(4,x_alice,alice_var,0,alice_err);
  }
  else{
    g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR_new, binWidth_fwd, fwdPR_new_Err);
    g_fwdPR_old = new TGraphErrors(2, x_fwd_old, fwdPR_old, binWidth_fwd_old, fwdPR_old_Err);
	g_alice = new TGraphErrors(4,x_alice,alice_var,binWidth_alice,alice_err);
  }
  TGraphErrors *g_fwdPRSys = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,binWidth_fwd,SysPR_fwd);
  TGraphErrors *g_fwdPR_oldSys = new TGraphErrors(2,x_fwd_old,fwdPR_old,binWidth_fwd_old,fwdPR_old_Sys);
  TGraphErrors *g_aliceSys = new TGraphErrors(4,x_alice,alice_var,binWidth_alice,alice_sys);


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
  g_midPR->GetXaxis()->SetLimits(0.,50.);
  g_midPR->SetMinimum(0.);
  g_midPR->SetMaximum(1.44);

  g_midPR->SetMarkerColor(kBlue+2);
  g_midPR->SetLineColor(kBlue+2);
  g_midPR->SetMarkerStyle(20);
  g_midPR->SetMarkerSize(1.4);
  g_midPRSys->SetLineColor(kBlue-4);
  g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);
  g_midPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPR_old->GetXaxis()->CenterTitle();
  g_midPR_old->GetYaxis()->SetTitle("R_{AA}");
  g_midPR_old->GetYaxis()->CenterTitle();
  g_midPR_old->SetTitle();
  g_midPR_old->GetXaxis()->SetLimits(0.,50.);
  g_midPR_old->SetMinimum(0.);
  g_midPR_old->SetMaximum(1.44);
  g_midPR_old->SetMarkerStyle(25);
  g_midPR_old->SetMarkerSize(1.4);
  g_midPR_old->SetMarkerColor(15);
  g_midPR_old->SetLineColor(15);

  g_midPR_oldSys->SetLineColor(15);
  g_midPR_oldSys->SetFillColorAlpha(15,0.4);

  g_midPR_Jpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midPR_Jpsi->GetXaxis()->CenterTitle();
  g_midPR_Jpsi->GetYaxis()->SetTitle("R_{AA}");
  g_midPR_Jpsi->GetYaxis()->CenterTitle();
  g_midPR_Jpsi->SetTitle();
  g_midPR_Jpsi->GetXaxis()->SetLimits(0.,50.);
  g_midPR_Jpsi->SetMinimum(0.);
  g_midPR_Jpsi->SetMaximum(1.44);
  g_midPR_Jpsi->SetMarkerStyle(33);
  g_midPR_Jpsi->SetMarkerSize(1.8);
  g_midPR_Jpsi->SetMarkerColor(38);
  g_midPR_Jpsi->SetLineColor(38);

  g_midPR_JpsiSys->SetLineColor(38);
  g_midPR_JpsiSys->SetFillColorAlpha(38,0.4);


  g_midPR->Draw("AP");
  g_midPR_old->Draw("P");
  g_midPR_Jpsi->Draw("P");
  if(isSys==1){
    g_midPRSys->Draw("5");
    g_midPR_oldSys->Draw("5");
    g_midPR_JpsiSys->Draw("5");
  }

  TLegend *leg1 = new TLegend(0.42,0.69,0.69,0.56);
	leg1->SetTextSize(text_size);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
  leg1->AddEntry(g_midPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%");
  leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, Cent. 0-100%");
  leg1->AddEntry(g_midPR_Jpsi,"#bf{Prompt J/#psi} HIN-16-025, Cent. 0-100%");
  leg1->Draw("SAME");
  jumSun(0,1,50,1);

  drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c1, iPeriod, iPos);

  c1->SaveAs(Form("./figs/compare_PR_mid_pT_Sys%d.pdf",isSys));

  TCanvas *c2 = new TCanvas("c2","",900,800);
  c2->cd();
  g_fwdPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_fwdPR->GetXaxis()->CenterTitle();
  g_fwdPR->GetYaxis()->SetTitle("R_{AA}");
  g_fwdPR->GetYaxis()->CenterTitle();
  g_fwdPR->SetTitle();
  g_fwdPR->GetXaxis()->SetLimits(0.,50.);
  g_fwdPR->SetMinimum(0.);
  g_fwdPR->SetMaximum(1.44);

  g_fwdPR->SetMarkerColor(kBlue+2);
  g_fwdPR->SetLineColor(kBlue+2);
  g_fwdPR->SetMarkerStyle(20);
  g_fwdPR->SetMarkerSize(1.4);
  g_fwdPRSys->SetLineColor(kBlue-4);
  g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

  g_fwdPR_old->SetMarkerStyle(25);
  g_fwdPR_old->SetMarkerSize(1.4);
  g_fwdPR_old->SetMarkerColor(15);
  g_fwdPR_old->SetLineColor(15);

  g_alice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_alice->GetXaxis()->CenterTitle();
  g_alice->GetYaxis()->SetTitle("R_{AA}");
  g_alice->GetYaxis()->CenterTitle();
  g_alice->SetTitle();
  g_alice->GetXaxis()->SetLimits(0.,50.);
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

  fwdPR_upper->SetLineColor(15);
  fwdPR_upper->SetLineWidth(2);

  g_fwdPR->Draw("AP");
  g_fwdPR_old->Draw("P");
  if(isSys==1){
    g_fwdPRSys->Draw("5");
    g_fwdPR_oldSys->Draw("5");
	g_aliceSys->Draw("5");
  }
  fwdPR_upper->Draw("");
  g_alice->Draw("P");

  TLegend *leg2 = new TLegend(0.22,0.69,0.49,0.56);
  leg2->SetTextSize(text_size);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
  leg2->AddEntry(g_fwdPR,"#bf{Prompt #psi(2S)} New Data, Cent. 0-90%");
  leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, Cent. 0-100%");
  leg2->AddEntry(g_alice,"#bf{Inclusive #psi(2S)} ALICE, 2.5 < y < 4, Cent.0-90%");
  leg2->Draw("SAME");
  jumSun(0,1,50,1);

  drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6< |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
  CMS_lumi_v2mass(c2, iPeriod, iPos);

  c2->SaveAs(Form("./figs/compare_PR_fwd_pT_Sys%d.pdf",isSys));
}
