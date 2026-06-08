#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
#include <iostream>

using namespace std;
using namespace RooFit;

valErr getFrac_psi2S_pbpb(int i=0);
void draw_Bfraction_pbpb_y0_2p4()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	int iPeriod = 101;
	int iPos = 33;

	const int nPtBins=11;
	TFile *fpp[nPtBins+1];

	double ptBin[nPtBins+1] = {6.5,7.5,8.5,9.5,11,13,15,17.5,20,25,30,50};
	double fracPP[nPtBins]; double fracPbPb[nPtBins];
	double fracErrPP[nPtBins]; double fracErrPbPb[nPtBins];
	double x[nPtBins]; double binWidth[nPtBins];

	double ptBin_HIN[12]={6.5,7.5,8.5,9.5,11,13,15,17.5,20,25,30,50};
	double var_HIN[11] = {0.272,0.288,0.329,0.356,0.396,0.450,0.483,0.537,0.525,0.547,0.588};
	double err_HIN[11] = {0.040,0.028,0.026,0.026,0.028,0.030,0.032,0.047,0.038,0.054,0.081};
	double x_HIN[11]; double binWidth_HIN[11];

	double BfracPP[nPtBins]; double BfracPbPb[nPtBins];

	TString kineLabel[nPtBins+1];

	TH1D *hfracPP = new TH1D("hfracPP",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);
	TH1D *hfracPbPb = new TH1D("hfracPbPb",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);


	for(int i=0;i<nPtBins;i++)
	{
		x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		valErr fracPP; 
		fracPP = getFrac_psi2S_pbpb(i);

		BfracPP[i] = fracPP.val; 

		fracErrPP[i] = fracPP.err;
		
		hfracPP->SetBinContent(i+1,fracPP.val);
		hfracPP->SetBinError(i+1,fracPP.err);

		cout << "pT	" << ptBin[i] << "	-	" << ptBin[i+1] << "	fraction PP :	" << fracPP.val << " +/- " << fracPP.err << endl;
	}

	for(int i=0; i<12; i++)
	{
		x_HIN[i] = (ptBin_HIN[i+1]+ptBin_HIN[i])/2;
		binWidth_HIN[i] = (ptBin_HIN[i+1]-ptBin_HIN[i])/2;
	}


	TGraphErrors* gfracPP = new TGraphErrors(nPtBins,x,BfracPP,binWidth,fracErrPP); 
	TGraphErrors* gfracHIN = new TGraphErrors(11,x_HIN,var_HIN,binWidth_HIN,err_HIN);

	int text_size=17;
	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	TCanvas *c1 = new TCanvas("c1","",550,520);
	c1->cd();
	gfracPP->SetMarkerColor(kBlue+2);
	gfracPP->SetLineColor(kBlue+2);
	gfracPP->SetMarkerStyle(20);
	gfracHIN->SetMarkerColor(kRed+2);
	gfracHIN->SetLineColor(kRed+2);
	gfracHIN->SetMarkerStyle(33);
	gfracHIN->SetMarkerSize(1.4);
	
	gfracPP->SetTitle("");
	gfracPP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gfracPP->GetYaxis()->SetTitle("non-propmt fraction");
	gfracPP->GetXaxis()->CenterTitle();
	gfracPP->GetYaxis()->CenterTitle();
	gfracPP->GetYaxis()->SetRangeUser(0,0.8);
	gfracPP->GetXaxis()->SetLimits(0,50);
	gfracHIN->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gfracHIN->GetYaxis()->SetTitle("non-propmt fraction");
	gfracHIN->GetXaxis()->CenterTitle();
	gfracHIN->GetYaxis()->CenterTitle();
	gfracHIN->GetXaxis()->SetLimits(0,50);
	gfracHIN->GetYaxis()->SetRangeUser(0,0.75);

	gfracPP->Draw("AP");
	gfracHIN->Draw("P");

	TLegend *leg = new TLegend(0.2,0.9,0.3,0.8);
	leg->SetTextSize(text_size);
	leg->SetTextFont(43);
	leg->SetBorderSize(0);
	leg->AddEntry(gfracPP,"PbPb 5.02 TeV (2018 DATA), |y|<2.4","PE");
	leg->AddEntry(gfracHIN,"PbPb 5.02 TeV (HIN-16-025), |y| < 2.4","P");
	leg->Draw("SAME");

	//CMS_lumi_v2mass(c1,iPeriod,iPos);

	c1->SaveAs("figs/Bfraction_pbpb_y0_2p4.pdf");


}
valErr getFrac_psi2S_pbpb(int i) {
	double ptBin[12] = {6.5,7.5,8.5,9.5,11,13,15,17.5,20,25,30,50};
	TString kineLabel[12];
	kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],0,2.4,0.0);
	TFile* inf = new TFile(Form("roots/2DFit_No_Weight/Final/2DFitResult_%s_centrality0-200_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
