#include <iostream>
#include "../../../../rootFitHeaders.h"
#include "../../../../commonUtility.h"
#include "../../../../JpsiUtility.h"
#include "../../../../cutsAndBin.h"
#include "../../../../CMS_lumi_v2mass.C"
#include "../../../../tdrstyle.C"
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

void mc_MassFit_CBGauss(
		double ptLow=6.5, double ptHigh=9.0,
		double yLow=0, double yHigh=1.6,
		int PR=0, //0=PR, 1=NP, 2=Inc.
		int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
		)
{
    TString DATE = "No_Weight";
    gStyle->SetEndErrorSize(0);
    gSystem->mkdir(Form("roots_MC"));
    gSystem->mkdir(Form("roots_MC/Mass_CBGauss"));
    gSystem->mkdir(Form("figs"));
    gSystem->mkdir(Form("figs/2DFit_%s",DATE.Data()),kTRUE);
	gSystem->mkdir(Form("figs/2DFit_%s/mc_Mass_CBGauss",DATE.Data()),kTRUE);

	TString bCont;
	if(PR==0) bCont="Prompt";
	else if(PR==1) bCont="NonPrompt";
	else if(PR==2) bCont="Inclusive";

	TString fname;
	if (PRw==1) fname="PR";
	else if (PRw==2) fname="NP";
	
	RooMsgService::instance().getStream(0).removeTopic(Caching);
	RooMsgService::instance().getStream(1).removeTopic(Caching);
	RooMsgService::instance().getStream(0).removeTopic(Plotting);
	RooMsgService::instance().getStream(1).removeTopic(Plotting);
	RooMsgService::instance().getStream(0).removeTopic(Integration);
	RooMsgService::instance().getStream(1).removeTopic(Integration);
	RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
	RooMsgService::instance().getStream(1).removeTopic(Fitting);
	RooMsgService::instance().getStream(1).removeTopic(Minimization);
	RooMsgService::instance().getStream(1).removeTopic(InputArguments);
	RooMsgService::instance().getStream(1).removeTopic(Eval);
	RooMsgService::instance().getStream(1).removeTopic(DataHandling);
	// RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooMsgService::instance().setSilentMode(true);

    // MC
	TFile* f1 = new TFile("../../../../skimmedFiles/OniaRooDataSet_isMC1_Psi2S_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230517.root", "read");

	// cout << "Input file: "
	// << Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",
	// fname.Data(),fEffW,fAccW,isPtW,isTnP) << endl;
	// exit(0);

	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>3.4 && mass<4.0",ptLow, ptHigh, yLow, yHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;

	massLow=3.4;
	massHigh=4.0;

	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<endl;
	cout << "####################################" << endl;
	RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset),WeightVar(*ws->var("weight")));
	RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->SetName("dsAB");
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();
	//***********************************************************************
	//****************************** MASS FIT *******************************
	//***********************************************************************

	//         The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
	//pt3-4.5
	//double paramsupper[8] = {0.4,    1.0,     4.9, 3.9, 1.0,     25.0};
	//double paramslower[8] = {0.01,   0.0,     1.1, 1.1, 0.0,    -25.0};//pt3-4.5 m_lambda==-25.0
	//double paramsupper[8] = {0.4,    1.0,     4.9, 2.9, 1.0,     25.0};
	//double paramslower[8] = {0.01,   0.0,     1., 1., 0.0,      0.0};//pt3-4.5 m_lambda==-25.0
	//Cent.10-20
    
	Double_t sigma_up, x_up, alpha_up, n_up, f_up;
	Double_t sigma_lo, x_lo, alpha_lo, n_lo, f_lo;

	sigma_up=1; x_up=2; alpha_up=10; n_up=10; f_up=1;
	sigma_lo=1e-3; x_lo=1e-3; alpha_lo=1e-3; n_lo=1e-3; f_lo=1e-3;

    
	//SIGNAL: initial params
	double sigma_1_init = 0.0555;
    double x_init = 0.366;
	double alpha_1_init = 2.9581;
	double n_1_init = 1.8425;
	double f_init = 0.6886;
    double sl1_mean = 0.01, sl2_mean = 0.04, sl3_mean = 0.06;
    double N_Jpsi_high = 2e+5; // 2500000
	double N_Bkg_high = 200000;
    double fit_limit = 3.86;
	double m_lambda_init = 5;
	double psi_2S_mass = pdgMass.Psi2S;

;   if(ptLow==3.5&&ptHigh==6.5) {
        N_Jpsi_high = 100000;
        sigma_up = 0.8;
    }
;   if(ptLow==3.5&&ptHigh==50) {
        N_Jpsi_high = 200000;
        sigma_up = 0.8;
    }
	if(ptLow==8&&ptHigh==12) { n_1_init = 5.3;}
	if(ptLow==12&&ptHigh==20) { n_1_init = 4.6;}

    if(ptLow==6.5&&ptHigh==9&&yLow==1.6) {
        N_Jpsi_high = 50000; // 2500000
        sigma_up = 0.8;
        n_1_init=10;
        n_up = 16;
        n_lo = 0;
    }
    if(ptLow==9&&ptHigh==12&&yLow==1.6) {
        N_Jpsi_high = 20000; // 2500000
        sigma_up = 0.8;
        n_1_init=8;
        n_up = 20;
        n_lo = 1e-6;
    }
    if(ptLow==6.5&&ptHigh==9&&yLow==0) {
        N_Jpsi_high = 50000; // 2500000
    }
     if(ptLow==5&&ptHigh==6.5) {
               N_Jpsi_high = 50000;
        f_init=0.8; sigma_1_init = 0.04; n_1_init = 1;
        x_init = 0.7; alpha_1_init = 3;
        x_up = 2;
        alpha_up=20;
        // m_lambda_init=2;
        n_up = 10;
    }
       if(ptLow==6.5&&ptHigh==12) {
        N_Jpsi_high = 200000;
        x_up = 2;
        f_init=0.3; sigma_1_init = 0.04; n_1_init = 2;
        x_init = 0.71; alpha_1_init = 2;
        m_lambda_init=2;
        n_up = 1000;
    }
    
     
     if(ptLow==20&&ptHigh==25) {
               N_Jpsi_high = 50000;
        f_init=0.8; sigma_1_init = 0.04; n_1_init = 1;
        x_init = 0.7; alpha_1_init = 3;
        x_up = 3;
        alpha_up=20;
        // m_lambda_init=2;
        n_up = 10;
     }
     if(ptLow==25&&ptHigh==30) {
               N_Jpsi_high = 50000;
        f_init=0.8; sigma_1_init = 0.04; n_1_init = 1;
        x_init = 0.7; alpha_1_init = 3;
        x_up = 3;
        alpha_up=20;
        // m_lambda_init=2;
        n_up = 10;
    }
    
     if(ptLow==30&&ptHigh==50) {
               N_Jpsi_high = 50000;
        f_init=0.8; sigma_1_init = 0.04; n_1_init = 1;
        x_init = 1; alpha_1_init = 3;
        x_up = 3;
        alpha_up=20;
        // m_lambda_init=2;
        n_up = 20;
    }


    //  if(ptLow==25&&ptHigh==30) {
    //     n_up =20; // 2500000
    //     n_lo = 1e-5;
    //     n_1_init = 2;
    // }
    //  if(ptLow==30&&ptHigh==50) {
    //     n_up =20; // 2500000
    //     n_lo = 1e-5;
    //     n_1_init = 2;
    // }

     if(ptLow==6.5&&ptHigh==50) {
        N_Jpsi_high = 400000; // 2500000
    }
    if(ptLow==9&&ptHigh==12&&yLow==0) {
        N_Jpsi_high = 35000;
    }
    if(ptLow==12&&ptHigh==15&&yLow==0) {
        N_Jpsi_high = 20000; // 2500000
    }


    if(ptLow==12&&ptHigh==50) {
        N_Jpsi_high = 10000;
        f_init=0.4; sigma_1_init = 0.03; n_1_init = 2;
        x_init = 0.4; alpha_1_init = 2;
    }
    if(ptLow==15&&ptHigh==20) {
        N_Jpsi_high = 500000;
        sigma_up=0.07;
        f_init=0.2766; sigma_1_init = 0.0527; n_1_init = 4.116;
        x_init = 0.8; alpha_1_init = 1.57;
        m_lambda_init=2;
    }
    if(ptLow==20&&ptHigh==50) {
        N_Jpsi_high = 100000;
        f_init=0.866; sigma_1_init = 0.0327; n_1_init = 1;
        x_init = 0.8; alpha_1_init = 1;
    }
    if(ptLow==20&&ptHigh==25) {
        N_Jpsi_high = 100000;
        sigma_up = 0.0513; sigma_lo = 0.0513;
        x_up = 0.5045; x_lo = 0.5045;
        alpha_up = 1.5956; alpha_lo = 1.5956;
        f_up = 0.3431; f_lo = 0.3431;
    }
    if(ptLow==25&&ptHigh==30) {
        N_Jpsi_high = 100000;
        sigma_up = 0.0513; sigma_lo = 0.0513;
        x_up = 0.5045; x_lo = 0.5045;
        alpha_up = 1.5956; alpha_lo = 1.5956;
        f_up = 0.3431; f_lo = 0.3431;
    }
    if(ptLow==30&&ptHigh==50) {
        N_Jpsi_high = 100000;
        sigma_up = 0.0513; sigma_lo = 0.0513;
        x_up = 0.5045; x_lo = 0.5045;
        alpha_up = 1.5956; alpha_lo = 1.5956;
        f_up = 0.3431; f_lo = 0.3431;
    }
    double paramsupper[8] = {sigma_up, x_up, alpha_up, n_up, f_up,  25.0};
    double paramslower[8] = {sigma_lo, x_lo, alpha_lo, n_lo, f_lo,  -5.0};

    // cout << psi_2S_mass << endl;
    // exit(1);

    //SIGNAL
    RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",psi_2S_mass, psi_2S_mass, psi_2S_mass) ;
    RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
    RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
    RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
    RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
    RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
    RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
    RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
    RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);
    //RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, f_init, f_init);
    //f->setConstant(kTRUE);

    //Set up crystal ball shapes
    RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
    RooAddPdf*  pdfMASS_Jpsi;

    //DOUBLE CRYSTAL BALL
    //RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
    //pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	
	RooGaussian* gauss = new RooGaussian("gauss","gaussian PDF",*(ws->var("mass")),mean,sigma_2_A);
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal",RooArgList(*cb_1_A,*gauss), RooArgList(*f) );

    //BACKGROUND
    RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[5], paramsupper[5]);
    RooRealVar *sl1 = new RooRealVar("sl1","sl1", sl1_mean, -10., 10.);
    RooRealVar *sl2 = new RooRealVar("sl2","sl2", sl2_mean, -10., 10.);
    RooRealVar *sl3 = new RooRealVar("sl3","sl3", sl3_mean, -10., 10.);

    //THIS IS THE BACKGROUND FUNCTION
    RooChebychev *pdfMASS_bkg;
    pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));

    //Build the model
    RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0, N_Jpsi_high);
    RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0, N_Bkg_high);
    //RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
    RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi),RooArgList(*N_Jpsi));
    //pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *bkg_1order),RooArgList(*N_Jpsi,*N_Bkg));
    //pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*N_JpsiPR,*N_JpsiNP,*N_Bkg));
    ws->import(*pdfMASS_Tot);

    //nMassBin = 25;
    TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
    c_A->cd();
    TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
    pad_A_1->SetTicks(1,1);
    pad_A_1->Draw(); pad_A_1->cd();
    RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
    myPlot_A->SetTitle("");
    ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

    pad_A_1->cd();
    bool logY_flag = true;
    if (logY_flag == true) {
        gPad->SetLogy();
    }
    RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
    dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));
    bool isWeighted = ws->data("dsAB")->isWeighted();
    cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
    RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), Range(3.4, fit_limit), NumCPU(nCPU));
    cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
	fitMass->Print("V");
    
    // Check and get fitted parameters
    const RooArgList & fitParams = fitMass->floatParsFinal();
    
    auto & fitN_Jpsi = (RooRealVar &) fitParams[0];
    auto & fitAlpha = (RooRealVar &) fitParams[1];
    auto & fitFraction = (RooRealVar &) fitParams[2];
    auto & fitMean = (RooRealVar &) fitParams[3];
    auto & fitN_1 = (RooRealVar &) fitParams[4];
    auto & fitSigma_1 = (RooRealVar &) fitParams[5];
    auto & fitX = (RooRealVar &) fitParams[6];
    // cout << fitN_Jpsi.getVal() << endl;
    
    /*
    for ( int i = 0; i < fitParams.getSize(); ++i)
    {
      auto & fitPar = (RooRealVar &) fitParams[i];
      std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
    }
    */
    
    double f_factor = (double) fitFraction.getVal();
    ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack), Range(3.4, fit_limit));
    // ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Range(3.4, 3.87), Normalization(fitFraction.getVal()));
    // ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Range(3.4, 3.87), Normalization(1-fitFraction.getVal()));
    ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Normalization(fitFraction.getVal()), Range(3.4, fit_limit));
    ws->pdf("gauss")->plotOn(myPlot2_A,Name("gauss"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(3.4, fit_limit));
    //ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(3.4, fit_limit));
    //ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));

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
    //myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);


    myPlot2_A->GetXaxis()->SetLabelSize(0);
    myPlot2_A->GetXaxis()->SetTitleSize(0);
    myPlot2_A->GetXaxis()->CenterTitle();
    //if (logY_flag == true) {
    //    myPlot2_A->SetMinimum(2*10);
    //    myPlot2_A->SetMaximum(10000000000);
    //}
    //else {
    //    myPlot2_A->GetYaxis()->SetRangeUser(0,YMax);
    //}
	myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);
    
    myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    myPlot2_A->Draw();
    
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c;",ptLow, ptHigh),text_x,text_y,text_color,text_size);
    if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
    else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
    drawText(Form("N_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
    // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
    drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(),ws->var("alpha_1_A")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    drawText(Form("f = %.4f #pm %.4f",fitFraction.getVal(),fitFraction.getError()),text_x,text_y-y_diff*4,text_color,text_size);
    drawText(Form("n_{1} = %.4f #pm %.4f",fitN_1.getVal(),fitN_1.getError()),text_x,text_y-y_diff*5,text_color,text_size);
    drawText(Form("#sigma_{1} = %.4f #pm %.4f",fitSigma_1.getVal(),fitSigma_1.getError()),text_x,text_y-y_diff*6,text_color,text_size);
    drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f",fitX.getVal(),fitX.getError()),text_x,text_y-y_diff*7,text_color,text_size);
    

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

    TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
    outh->GetXaxis()->SetBinLabel(1,"Jpsi");

    float temp1 = ws->var("N_Jpsi")->getVal();
    float temp1err=ws->var("N_Jpsi")->getError();

    outh->SetBinContent(1,temp1);
    outh->SetBinError(1,temp1err);

    fitMass->Print();
    Double_t theNLL = fitMass->minNll();
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
    TString kineLabel = getKineLabelpp (ptLow, ptHigh,yLow, yHigh, 0.0);
    
    TFile* outFile;
    outFile = new TFile(Form("roots_MC/Mass_CBGauss/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
    c_A->SaveAs(Form("figs/2DFit_%s/mc_Mass_CBGauss/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));


    pdfMASS_Tot->Write();
    datasetMass->Write();
    outh->Write();
    // ws->Write();
    fitMass->Write();
    outFile->Close();
        fitMass->Print("v");
    // Get # of entries for every bins.
    // To compare events number point by point
    // Due to there was events number issue
    // ofstream fout(Form("entry_check_mass_fit_of_2D_fit_Psi2S_Inclusive_%s.txt", kineLabel.Data()));
    // for (int i = 1; i <= h->GetNbinsX(); i++)
    // fout << h->GetBinContent(i) << endl;
}
