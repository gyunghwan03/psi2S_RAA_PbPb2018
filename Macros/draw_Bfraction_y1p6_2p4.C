#include <iostream>
#include "../rootFitHeaders.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
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

using namespace std;
using namespace RooFit;

valErr getFrac_psi2S_PbPb(int i=0);
valErr getFrac_psi2S_pp(int i=0);
void draw_Bfraction_y1p6_2p4()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	int iPeriod = 101;
	int iPos = 33;

	const int nPtBins=4;
	TFile *fPbPb[nPtBins+1];
	TFile *fpp[nPtBins+1];

	double ptBin[nPtBins+1] = {3,4.5,6.5,12,50};
	double fracPP[nPtBins]; double fracPbPb[nPtBins];
	double fracErrPP[nPtBins]; double fracErrPbPb[nPtBins];
	double x[nPtBins]; double binWidth[nPtBins];

	double BfracPP[nPtBins]; double BfracPbPb[nPtBins];

	TString kineLabel[nPtBins+1];

	TH1D *hfracPP = new TH1D("hfracPP",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);
	TH1D *hfracPbPb = new TH1D("hfracPbPb",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);


	for(int i=0;i<nPtBins;i++)
	{
		x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		valErr fracPP; valErr fracPbPb;
		fracPP = getFrac_psi2S_pp(i);
		fracPbPb = getFrac_psi2S_PbPb(i);
		
		BfracPP[i] = fracPP.val; BfracPbPb[i] = fracPbPb.val;

		fracErrPP[i] = fracPP.err;
		fracErrPbPb[i] = fracPbPb.err;

		hfracPP->SetBinContent(i+1,fracPP.val);
		hfracPP->SetBinError(i+1,fracPP.err);
		hfracPbPb->SetBinContent(i+1,fracPbPb.val);
		hfracPbPb->SetBinError(i+1,fracPbPb.err);

		//gfracPP->SetPoint(i,(ptBin[i]+ptBin[i+1])/2, fracPP.val);
		//gfracPP->SetPointError(i,i,fracPP.err);
		//gfracPbPb->SetPoint(i,(ptBin[i]+ptBin[i+1])/2, fracPbPb.val);
		//gfracPbPb->SetPointError(i,i,fracPbPb.err);

		cout << "pT	" << ptBin[i] << "	-	" << ptBin[i+1] << "	fraction PP :	" << fracPP.val << " +/- " << fracPP.err << endl;
		cout << "				fraction PbPb :	" << fracPbPb.val << " +/- " << fracPbPb.err << endl;
	}

	TGraphErrors* gfracPP = new TGraphErrors(nPtBins,x,BfracPP,binWidth,fracErrPP); 
	TGraphErrors* gfracPbPb = new TGraphErrors(nPtBins,x,BfracPbPb, binWidth,fracErrPbPb);

	int text_size=17;
	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	TCanvas *c1 = new TCanvas("c1","",550,520);
	c1->cd();
	gfracPP->SetMarkerColor(kBlue+2);
	gfracPP->SetLineColor(kBlue+2);
	gfracPP->SetMarkerStyle(20);
	gfracPbPb->SetMarkerColor(kRed+2);
	gfracPbPb->SetLineColor(kRed+2);
	gfracPbPb->SetMarkerStyle(21);
	
	gfracPP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gfracPP->GetYaxis()->SetTitle("non-propmt fraction");
	gfracPP->GetXaxis()->CenterTitle();
	gfracPP->GetYaxis()->CenterTitle();
	gfracPbPb->GetXaxis()->CenterTitle();
	gfracPbPb->GetYaxis()->CenterTitle();
	gfracPP->GetYaxis()->SetRangeUser(0,0.75);
	gfracPbPb->GetYaxis()->SetRangeUser(0,0.75);
	gfracPP->GetXaxis()->SetLimits(0,50);
	gfracPP->Draw("AP");
	gfracPbPb->Draw("P");

	TLegend *leg = new TLegend(0.2,0.9,0.3,0.8);
	leg->SetTextSize(text_size);
	leg->SetTextFont(43);
	leg->SetBorderSize(0);
	leg->AddEntry(gfracPP,"pp 5.02 TeV, 1.6<|y|<2.4","P");
	leg->AddEntry(gfracPbPb,"PbPb 5.02 TeV, 1.6<|y|<2.4","P");
	leg->Draw("SAME");

	CMS_lumi_v2mass(c1,iPeriod,iPos);

	c1->SaveAs("Bfraction_y1p6_2p4.pdf");


}
valErr getFrac_psi2S_PbPb(int i) {
	double ptBin[5] = {3,4.5,6.5,12,50};
	TString kineLabel[5];
	kineLabel[i] = getKineLabel(ptBin[i],ptBin[i+1],1.6,2.4,0.0,0,200);
	TFile* inf = new TFile(Form("./psi2S_Corr/roots/2DFit_230324/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
valErr getFrac_psi2S_pp(int i) {
	double ptBin[5] = {3,4.5,6.5,12,50};
	TString kineLabel[5];
	kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],1.6,2.4,0.0);
	TFile* inf = new TFile(Form("./pp_psi2S_Corr/roots/2DFit_230323/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
