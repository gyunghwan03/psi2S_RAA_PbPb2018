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

void compare_pT_NP(bool isSys=true)
{
  gStyle->SetOptStat(0);
  setTDRStyle();

  TFile *f_mid = new TFile("roots/RAA_psi2S_midRap_pT.root"); 
  TFile *f_fwd = new TFile("roots/RAA_psi2S_forRap_pT.root"); 
  TFile *fSys = new TFile("../syst_summary_psi2S/syst_roots/total_syst.root");
  TFile *f_old = new TFile("./roots/Raa_NP_Jpsi_HIN_16_025.root");

  TH1D *h_midNP = (TH1D*) f_mid->Get("hRAA_NP");
  TH1D *h_fwdNP = (TH1D*) f_fwd->Get("hRAA_NP");

  TH1D *hSys_NP = (TH1D*) fSys->Get("mid_pt_NP");
  TH1D *hSys_fwd_NP = (TH1D*) fSys->Get("fwd_pt_NP");

  TH1F *hOld = (TH1F*) f_old->Get("Table 27/Hist1D_y1");
  TH1F *hOld_Err = (TH1F*) f_old->Get("Table 27/Hist1D_y1_e1");
  TH1F *hOld_Sys = (TH1F*) f_old->Get("Table 27/Hist1D_y1_e2");


  const int nPtBins=6;
  const int nPtBins_old=12;
  const int nPtBins_fwd=5;
  double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,50};
  double ptBin_fwd[nPtBins_fwd+1] = {3.5,5,6.5,9,12,50};
  double ptBin_old[nPtBins_old+1] = {6.5,7.5,8.5,9.5,11,13,15,17.5,20,25,30,35,50};
  double ptBin_fwd_old[nPtBins_fwd-1] = {6.5,12,50};
  double x[nPtBins]; double binWidth[nPtBins]; 
  double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];
  double x_old[nPtBins_old]; double binWidth_old[nPtBins_old];
  double x_fwd_old[2]={(6.5+12)/2,(12+30)/2}; double binWidth_fwd_old[2]={(12-6.5)/2,(30-12)/2};

  double midNP_new[nPtBins]; 
  double midNP_new_Err[nPtBins]; 
  double midNP_new_Sys[nPtBins]; 
  double fwdNP_new[nPtBins_fwd]; 
  double fwdNP_new_Err[nPtBins_fwd]; 
  double fwdNP_new_Sys[nPtBins_fwd]; 

  double SysNP[nPtBins]; 
  double SysNP_fwd[nPtBins_fwd];


  double midNP_old[nPtBins_old];
  double midNP_old_Err[nPtBins_old];
  double midNP_old_Sys[nPtBins_old];

  for (int i=0; i<nPtBins_old; i++){
    midNP_old[i] = hOld->GetBinContent(i+1);
    midNP_old_Err[i] = hOld_Err->GetBinContent(i+1);
    midNP_old_Sys[i] = hOld_Sys->GetBinContent(i+1);
    x_old[i] = (ptBin_old[i+1]+ptBin_old[i])/2;
    binWidth_old[i] = (ptBin_old[i+1]-ptBin_old[i])/2;
    cout << "Val : " << midNP_old[i] << " Err : " << midNP_old_Err[i] << " Sys : " << midNP_old_Sys[i] << " x : " << x_old[i] << " bin width : " << binWidth_old[i] << endl;
  }

  //double fwdPR_old[2] = {0.135,0.247};
  //double fwdPR_old_Err[2] = {0.079,0.098}; 
  //double fwdPR_old_Sys[2] = {0.038,0.042};

  for (int i=0; i<nPtBins; i++){
    midNP_new[i] = h_midNP->GetBinContent(i+1);
    midNP_new_Err[i] = h_midNP->GetBinError(i+1);
    x[i] = (ptBin[i+1]+ptBin[i])/2;
    binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
    SysNP[i] = midNP_new[i]*(hSys_NP->GetBinContent(i+1));
  }
  //for (int i=0; i<nPtBins_fwd; i++){
  //  fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
  //  fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
  //  x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
  //  binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
  //  SysPR_fwd[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
  //}

  TGraphErrors *g_midNP;
  TGraphErrors *g_midNP_old;

  if (isSys == 1)
  {
    g_midNP = new TGraphErrors(nPtBins, x, midNP_new, 0, midNP_new_Err);
    g_midNP_old = new TGraphErrors(nPtBins_old, x_old, midNP_old, 0, midNP_old_Err);
  }
  else
  {
    g_midNP = new TGraphErrors(nPtBins, x, midNP_new, binWidth, midNP_new_Err);
    g_midNP_old = new TGraphErrors(nPtBins_old, x_old, midNP_old, binWidth_old, midNP_old_Err);
  }
  TGraphErrors *g_midNPSys = new TGraphErrors(nPtBins,x,midNP_new,binWidth,SysNP);
  TGraphErrors *g_midNP_oldSys = new TGraphErrors(nPtBins_old,x_old,midNP_old,binWidth_old,midNP_old_Sys);

  //TGraphErrors *g_fwdPR = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,0,fwdPR_new_Err);
  //TGraphErrors *g_fwdPRSys = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,binWidth_fwd,SysPR_fwd);
  //TGraphErrors *g_fwdPR_old = new TGraphErrors(2,x_fwd_old,fwdPR_old,0,fwdPR_old_Err);
  //TGraphErrors *g_fwdPR_oldSys = new TGraphErrors(2,x_fwd_old,fwdPR_old,binWidth_fwd_old,fwdPR_old_Sys);
  //double x_1st = (6.5+3)/2;
  //TArrow *fwdPR_upper = new TArrow(x_1st,0,x_1st,0.395,0.027,"<-|");

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

  g_midNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midNP->GetXaxis()->CenterTitle();
  g_midNP->GetYaxis()->SetTitle("R_{AA}");
  g_midNP->GetYaxis()->CenterTitle();
  g_midNP->SetTitle();
  g_midNP->GetXaxis()->SetLimits(0.,50.);
  g_midNP->SetMinimum(0.);
  g_midNP->SetMaximum(1.44);

  g_midNP->SetMarkerColor(kRed+3);
  g_midNP->SetLineColor(kRed+3);
  g_midNP->SetMarkerStyle(20);
  g_midNP->SetMarkerSize(1.4);
  g_midNPSys->SetLineColor(kRed-4);
  g_midNPSys->SetFillColorAlpha(kRed-9,0.40);
  g_midNP_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_midNP_old->GetXaxis()->CenterTitle();
  g_midNP_old->GetYaxis()->SetTitle("R_{AA}");
  g_midNP_old->GetYaxis()->CenterTitle();
  g_midNP_old->SetTitle();
  g_midNP_old->GetXaxis()->SetLimits(0.,50.);
  g_midNP_old->SetMinimum(0.);
  g_midNP_old->SetMaximum(1.44);
  g_midNP_old->SetMarkerStyle(25);
  g_midNP_old->SetMarkerSize(1.4);
  g_midNP_old->SetMarkerColor(15);
  g_midNP_old->SetLineColor(15);

  g_midNP_oldSys->SetLineColor(15);
  g_midNP_oldSys->SetFillColorAlpha(15,0.4);


  g_midNP->Draw("AP");
  g_midNP_old->Draw("P");
  if(isSys==1){
    g_midNPSys->Draw("5");
    g_midNP_oldSys->Draw("5");
  }

  TLegend *leg1 = new TLegend(0.25,0.69,0.69,0.56);
	leg1->SetTextSize(text_size);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
  leg1->AddEntry(g_midNP,"#bf{NonPrompt #psi(2S)} New Data, |y| < 1.6, Cent. 0-90%", "pe");
  leg1->AddEntry(g_midNP_old,"#bf{NonPrompt J/#psi} HIN-16-025, |y| < 2.4, Cent. 0-100%", "pe");
  leg1->Draw("SAME");
  jumSun(0,1,50,1);

  drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff, text_color, text_size);
	//drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
  CMS_lumi_v2mass(c1, iPeriod, iPos);

  c1->SaveAs(Form("./figs/compare_NP_mid_pT_Sys%d.pdf",isSys));
/*
  TCanvas *c2 = new TCanvas("c2","",900,800);
  c2->cd();
  g_fwdNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_fwdNP->GetXaxis()->CenterTitle();
  g_fwdNP->GetYaxis()->SetTitle("R_{AA}");
  g_fwdNP->GetYaxis()->CenterTitle();
  g_fwdNP->SetTitle();
  g_fwdNP->GetXaxis()->SetLimits(0.,50.);
  g_fwdNP->SetMinimum(0.);
  g_fwdNP->SetMaximum(1.44);

  g_fwdNP->SetMarkerColor(kBlue+2);
  g_fwdNP->SetLineColor(kBlue+2);
  g_fwdNP->SetMarkerStyle(20);
  g_fwdNP->SetMarkerSize(1.4);
  g_fwdNPSys->SetLineColor(kBlue-4);
  g_fwdNPSys->SetFillColorAlpha(kBlue-9,0.40);

  g_fwdNP_old->SetMarkerStyle(25);
  g_fwdNP_old->SetMarkerSize(1.4);
  g_fwdNP_old->SetMarkerColor(15);
  g_fwdNP_old->SetLineColor(15);

  g_fwdNP_oldSys->SetLineColor(15);
  g_fwdNP_oldSys->SetFillColorAlpha(15,0.4);

  fwdNP_upper->SetLineColor(15);
  fwdNP_upper->SetLineWidth(2);

  g_fwdNP->Draw("AP");
  g_fwdNPSys->Draw("5");
  g_fwdNP_old->Draw("P");
  g_fwdNP_oldSys->Draw("5");
  fwdNP_upper->Draw("");

  TLegend *leg2 = new TLegend(0.62,0.835,0.72,0.735);
  leg2->SetTextSize(text_size);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
  leg2->AddEntry(g_fwdNP,"Prompt #psi(2S) New Data");
  leg2->AddEntry(g_fwdNP_old,"Propmt #psi(2S) HIN-16-025");
  leg2->Draw("SAME");
  jumSun(0,1,50,1);

  drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("1.6< |y| < 2.4", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
  CMS_lumi_v2mass(c2, iPeriod, iPos);

  c2->SaveAs("./figs/compare_fwd_pT.pdf");
*/
}
