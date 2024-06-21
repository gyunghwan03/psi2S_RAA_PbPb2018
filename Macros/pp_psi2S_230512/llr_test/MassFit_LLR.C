#include <iostream>
#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
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
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void MassFit_LLR(
		float ptLow=3, float ptHigh=4.5,
		float yLow=1.6, float yHigh=2.4,
		int ordC=1,
		int PR=2, //0=PR, 1=NP, 2=Inc.
		float ctau3DErrCut=1.
		)
{
	TString DATE="LLR";
	gStyle->SetEndErrorSize(0);
	gSystem->mkdir(Form("roots"));
	gSystem->mkdir(Form("figs"));
	gSystem->mkdir(Form("roots/2DFit_%s/",DATE.Data()));
	gSystem->mkdir(Form("figs/2DFit_%s",DATE.Data()));
	gSystem->mkdir(Form("figs/2DFit_%s/Mass/",DATE.Data()));

	TString bCont;
	if(PR==0) bCont="Prompt";
	else if(PR==1) bCont="NonPrompt";
	else if(PR==2) bCont="Inclusive";

	TString fOrd;
	if (ordC==1) fOrd="1stCheby";
	else if (ordC==2) fOrd="2ndCheby";
	else if (ordC==3) fOrd="3rdCheby";
	else if (ordC==4) fOrd="4thCheby";
	else if (ordC==5) fOrd="5thCheby";
	else if (ordC==6) fOrd="6thCheby";
	else if (ordC==0) fOrd="Expo";

	cout << "Order : " << fOrd << endl;

	RooMsgService::instance().getStream(0).removeTopic(Caching);
	RooMsgService::instance().getStream(1).removeTopic(Caching);
	RooMsgService::instance().getStream(0).removeTopic(Plotting);
	RooMsgService::instance().getStream(1).removeTopic(Plotting);
	RooMsgService::instance().getStream(0).removeTopic(Integration);
	RooMsgService::instance().getStream(1).removeTopic(Integration);
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

	TFile* f1; TFile* f2; TFile* f3;
	TString kineCut;
	massLow=3.3;
	massHigh=4.1;

	f1 = new TFile("../../../skimmedFiles/OniaRooDataSet_isMC0_Psi2S_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230515.root");
	//f1 = new TFile("../../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20210111.root");
	//f1 = new TFile("../../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20201127.root");
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>3.3 && mass<4.1",ptLow, ptHigh, yLow, yHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	TString nan_cut = "&& !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes)";
	kineCut = OS+accCut+kineCut + nan_cut;
	TString kineLabel = getKineLabelpp (ptLow, ptHigh,yLow, yHigh, 0.0);

	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<endl;
	cout << "####################################" << endl;
	//RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
	RooDataSet *dsAB = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->SetName("dsAB");
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();
	//***********************************************************************
	//****************************** MASS FIT *******************************
	//***********************************************************************

	//         The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
	double paramsupper[8] = {0.4,    1.0,     4.9, 5.9, 1.0,     25.};
	double paramslower[8] = {0.01,   0.0,     0.1, 0.7, 0.0,    -25};
	//SIGNAL: initial params
	//double sigma_1_init = 0.04;
	//double x_init = 0.23;
	//double alpha_1_init = 2.6;
	//double n_1_init = 3.6;
	//double f_init = 0.4;
	double m_lambda_init = 2;

	//Bkg : initial params
	double sl1_init = 0.01;
	double sl2_init = 0.01;
	double sl3_init = 0.01;
	double sl4_init = 0.01;
	double sl5_init = 0.01;
	double sl6_init = 0.01;

	double sig_high = 1000000;
	double bkg_high = 1000000;

	if (ptLow==3.5&&ptHigh==6.5) {
		sig_high = 100000;
		bkg_high = 5000000;
	} else if(ptLow==6.5&&ptHigh==9&&yLow==1.6) {
		sig_high = 1000000;
		bkg_high = 25000000;
	} else if(ptLow==6.5&&ptHigh==9&&yLow==0) {
		sig_high = 50000;
		bkg_high = 3000000;
	} else if(ptLow==9&&ptHigh==12&&yLow==1.6) {
		sig_high = 370000;
		bkg_high = 800000;
	} else if(ptLow==12&&ptHigh==50) {
		sig_high = 20000;
		bkg_high = 1000000;
	} else if(ptLow==3.5&&ptHigh==50) {
		sig_high = 200000;
		bkg_high = 10000000;
	} else if(ptLow==9&&ptHigh==12&&yLow==0) {
		sig_high = 40000;
		bkg_high = 1000000;
	} else if(ptLow==12&&ptHigh==15) {
		sig_high = 20000;
		bkg_high = 10000000;
	} else if(ptLow==15&&ptHigh==20) {
		sig_high = 20000;
		bkg_high = 8000000;
	} else if(ptLow==20&&ptHigh==25) {
		sig_high = 10000;
		bkg_high = 5000000;
	} else if(ptLow==25&&ptHigh==50) {
		sig_high = 100000;
		bkg_high = 5000000;
	} else if(ptLow==6.5&&ptHigh==50) {
		sig_high = 3000000;
		bkg_high = 10000000;
	}

	//SIGNAL
	TFile * f_fit = new TFile(Form("../roots_MC/Mass/mc_MassFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",kineLabel.Data()));
	RooDataSet *dataset_fit = (RooDataSet*)f_fit->Get("datasetMass");
	RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
	ws_fit->import(*dataset_fit);

	Double_t alpha_MC_value = ws_fit->var("alpha_1_A")->getVal();
	Double_t alpha_MC_value_err = ws_fit->var("alpha_1_A")->getError();
	Double_t n_MC_value = ws_fit->var("n_1_A")->getVal();
	Double_t n_MC_value_err = ws_fit->var("n_1_A")->getError();
	Double_t xA_MC_value = ws_fit->var("x_A")->getVal();
	Double_t xA_MC_value_err =  ws_fit->var("x_A")->getError();
	Double_t f_MC_value = ws_fit->var("f")->getVal();
	Double_t f_MC_value_err = ws_fit->var("f")->getError();
	Double_t sigma_MC_value = ws_fit->var("sigma_1_A")->getVal();
	Double_t sigma_MC_value_err = ws_fit->var("sigma_1_A")->getError();

	double sigma_index = 5;

	Double_t alpha_lower = alpha_MC_value-(sigma_index*alpha_MC_value_err);
	Double_t alpha_higher = alpha_MC_value+(sigma_index*alpha_MC_value_err);
	Double_t xA_lower = xA_MC_value-(sigma_index*xA_MC_value_err);
	Double_t xA_higher = xA_MC_value+(sigma_index*xA_MC_value_err);
	Double_t n_lower = n_MC_value-(sigma_index*n_MC_value_err);
	Double_t n_higher = n_MC_value+(sigma_index*n_MC_value_err);
	Double_t f_lower = f_MC_value-(sigma_index*f_MC_value_err);
	Double_t f_higher = f_MC_value+(sigma_index*f_MC_value_err);
	Double_t sigma_lower = sigma_MC_value-(sigma_index*sigma_MC_value_err);
	Double_t sigma_higher = sigma_MC_value+(sigma_index*sigma_MC_value_err);
	
	double alpha_1_init = alpha_MC_value; double n_1_init = n_MC_value;
	double sigma_1_init = sigma_MC_value; double x_init = xA_MC_value; double f_init = f_MC_value;

	//SIGNAL
	RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.Psi2S, pdgMass.Psi2S-0.01, pdgMass.Psi2S+0.01) ;
	//RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.Psi2S, pdgMass.Psi2S -0.05, pdgMass.Psi2S + 0.05) ;
	//RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
	//RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
	//RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
	//RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
	//RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
	//RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
	//RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
	//RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);

	RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init);
	RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
	RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
	RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init);
	RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
	RooRealVar    n_1_A("n_1_A","power order", n_1_init);
	RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
	RooRealVar   *f = new RooRealVar("f","cb fraction", f_init);

	//Set up crystal ball shapes
	RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
	RooAddPdf*  pdfMASS_Jpsi;
	//DOUBLE CRYSTAL BALL
	RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	//BACKGROUND
	RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[5], paramsupper[5]);
	RooRealVar *sl1 = new RooRealVar("sl1","sl1", sl1_init, -1., 1.);
	RooRealVar *sl2 = new RooRealVar("sl2","sl2", sl2_init, -1., 1.);
	RooRealVar *sl3 = new RooRealVar("sl3","sl3", sl3_init, -1., 1.);
	RooRealVar *sl4 = new RooRealVar("sl4","sl4", sl4_init, -1., 1.);
	RooRealVar *sl5 = new RooRealVar("sl5","sl5", sl5_init, -1., 1.);
	RooRealVar *sl6 = new RooRealVar("sl6","sl6", sl6_init, -1., 1.);
	//RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",10000,0,1000000);
	//THIS IS THE BACKGROUND FUNCTION
	//RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
	//RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)*@2+@3",RooArgList(*(ws->var("mass")), m_lambda_A, *sl1, *sl2));
	//RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
	RooChebychev *pdfMASS_bkg;
	RooGenericPdf *pdfMASS_bkg_expo;
	if (ordC==1) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1)); }
	else if (ordC==2) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2)); }
	else if (ordC==3) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3)); }
	else if (ordC==4) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3, *sl4)); }
	else if (ordC==5) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3, *sl4, *sl5)); }
	else if (ordC==6) { pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3, *sl4, *sl5, *sl6)); }
	else if (ordC==0) { pdfMASS_bkg_expo = new RooGenericPdf("pdfMASS_bkg_expo","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));  }
	//Build the model
	//RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals", 0, 5000000);
	//RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg", 0, 1000000);
	RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals", 0, sig_high);
	RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg", 0, bkg_high);
	RooAddPdf* pdfMASS_Tot;
	if (ordC==0) { pdfMASS_Tot	= new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg_expo),RooArgList(*N_Jpsi,*N_Bkg)); }
	else { pdfMASS_Tot  = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg)); }
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *bkg_1order),RooArgList(*N_Jpsi,*N_Bkg));
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*N_JpsiPR,*N_JpsiNP,*N_Bkg));
	ws->import(*pdfMASS_Tot);
	TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
	c_A->cd();
	TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
	pad_A_1->SetTicks(1,1);
	pad_A_1->Draw(); pad_A_1->cd();
	RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
	myPlot_A->SetTitle("");
	ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

	pad_A_1->cd();
	gPad->SetLogy();
	RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
	dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));
	cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
	RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), NumCPU(nCPU) );
	cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
	ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack));
	//ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
	if (ordC==0) ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Bkg_expo"),Components(RooArgSet(*pdfMASS_bkg_expo)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
	else ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
	//make a pretty plot
	myPlot2_A->SetFillStyle(4000);
	myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
	//myPlot2_A->GetYaxis()->CenterTitle();
	//myPlot2_A->GetYaxis()->SetTitleSize(0.058);
	//myPlot2_A->GetYaxis()->SetLabelSize(0.054);
	//myPlot2_A->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());
	TH1* h = ws->data("dsAB")->createHistogram("hist", *ws->var("mass"), Binning(myPlot_A->GetNbinsX(),myPlot_A->GetXaxis()->GetXmin(),myPlot_A->GetXaxis()->GetXmax()));
	Double_t YMax = h->GetBinContent(h->GetMaximumBin());
	// Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
	Double_t YMin = 1e99;
	for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
	double Ydown;
	double Yup;
	Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
	Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
	myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);

	//myPlot2_A->SetMinimum(2*10);
	myPlot2_A->GetXaxis()->SetLabelSize(0);
	myPlot2_A->GetXaxis()->SetTitleSize(0);
	myPlot2_A->GetXaxis()->CenterTitle();
	myPlot2_A->GetXaxis()->SetRangeUser(massLow, massHigh);
	myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	myPlot2_A->Draw();

	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("n_{#psi (2S)} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
	drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);

	TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
	c_A->cd();
	pad_A_2->Draw();
	pad_A_2->cd();
	pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
	pad_A_2->SetBottomMargin(0.67);
	pad_A_2->SetBottomMargin(0.4);
	pad_A_2->SetFillStyle(4000);
	pad_A_2->SetFrameFillStyle(4000);
	pad_A_2->SetTicks(1,1);

	RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");
	RooHist* hpull_A = frameTMP->pullHist("dataOS","pdfMASS_tot", true);
	hpull_A->SetMarkerSize(0.8);
	RooPlot* pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution")) ;
	pullFrame_A->addPlotable(hpull_A,"P") ;
	pullFrame_A->SetTitle("");
	pullFrame_A->SetTitleSize(0);
	pullFrame_A->GetYaxis()->SetTitleOffset(0.3) ;
	pullFrame_A->GetYaxis()->SetTitle("Pull") ;
	pullFrame_A->GetYaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetYaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetYaxis()->SetRangeUser(-3.8,3.8);
	pullFrame_A->GetYaxis()->CenterTitle();

	pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	pullFrame_A->GetXaxis()->SetTitleOffset(1.05) ;
	pullFrame_A->GetXaxis()->SetLabelOffset(0.04) ;
	pullFrame_A->GetXaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetXaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetXaxis()->CenterTitle();

	pullFrame_A->GetYaxis()->SetTickSize(0.04);
	pullFrame_A->GetYaxis()->SetNdivisions(404);
	pullFrame_A->GetXaxis()->SetTickSize(0.03);
	pullFrame_A->Draw() ;

	TLine *l1 = new TLine(massLow,0,massHigh,0);
	l1->SetLineStyle(1);
	l1->Draw("same");

	printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_tot", nMassBin, false);

	fitMass->Print();
	Double_t theNLL = fitMass->minNll();
	cout.precision(15);
	cout << " *** NLL : " << theNLL << endl;
	//RooRealVar *Par1 = ws->var("sl1");
	//RooRealVar *Par2 = ws->var("sl2");
	//cout << "Chebychev Par1 : " << Par1->getVal() << " +/- " << Par1->getError() << endl;
	//cout << "Chebychev Par2 : " << Par2->getVal() << " +/- " << Par2->getError() << endl;

	//ws->Print();
	RooArgSet* fitargs = new RooArgSet();
	fitargs->add(fitMass->floatParsFinal());
	//RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", RooArgList(*ws->var("N_Jpsi"),*ws->var("N_Bkg")) );
	RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", *fitargs );
	datasetMass->add(*fitargs);
	RooWorkspace *wsmass = new RooWorkspace("workspaceMass");

	c_A->Update();
	

	TTree* tree = new TTree("tree","");
	tree->Branch("theNLL", &theNLL, "theNLL/D");
	tree->Fill();

	TH1D* hLLR = new TH1D("hLLR", "", 1,0,1);
	hLLR->SetBinContent(1,theNLL);

	TFile* outFile;
	if(PR==2){
		outFile = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%s_%s.root",DATE.Data(), bCont.Data(), kineLabel.Data(), fOrd.Data()),"recreate");
		c_A->SaveAs(Form("figs/2DFit_%s/Mass/Mass_%s_%s_%s.pdf",DATE.Data(), bCont.Data(), kineLabel.Data(),fOrd.Data())); }
	//  outFile->cd();
	pdfMASS_Tot->Write();
	datasetMass->Write();
	fitMass->Write();
	outFile->Close();

	TFile* LLRFile;
	LLRFile = new TFile(Form("roots/2DFit_%s/NLL_MassFit_%s_%s_%s.root",DATE.Data(), bCont.Data(), kineLabel.Data(), fOrd.Data()),"recreate");
	LLRFile->cd();
	hLLR->Write();
	LLRFile->Close();

	TString txtName = Form("figs/2DFit_%s/Mass/Mass_%s_%s_%s.txt",DATE.Data(), bCont.Data(), kineLabel.Data(), fOrd.Data());
	ofstream txtFile;
	txtFile.open(txtName.Data());
	txtFile.precision(15);
	//fitMass->floatParsFinal().printMultiline(txtFile, 1111);
	//txtFile << "\n" ;
	//if (ordC==1) { txtFile << "sl1_init : " <<    sl1_init; }
	//else if (ordC==2) {
	//	txtFile << "sl1_init : " <<   sl1_init << "\n";
	//	txtFile << "sl2_init : " <<   sl2_init << "\n";
	//}
	//else if (ordC==3) {
	//	txtFile << "sl1_init : " <<   sl1_init << "\n";
	//	txtFile << "sl2_init : " <<   sl2_init << "\n";
	//	txtFile << "sl3_init : " <<   sl3_init << "\n";
	//}
	//else if (ordC==4) {
	//	txtFile << "sl1_init : " <<   sl1_init << "\n";
	//	txtFile << "sl2_init : " <<   sl2_init << "\n";
	//	txtFile << "sl3_init : " <<   sl3_init << "\n";
	//	txtFile << "sl4_init : " <<   sl4_init << "\n";
	//}
	////else if (ordC==0) {
	////	txtFile << "sl4_init : " <<   sl4_init << "\n";
	////}
	//txtFile << "\n";
	//txtFile << " *** NLL : " <<  theNLL << "\n";
	txtFile << theNLL;
	txtFile.close();
}

