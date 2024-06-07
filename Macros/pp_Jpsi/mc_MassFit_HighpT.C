#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
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
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
using namespace std;
using namespace RooFit;


void check_convergence(RooFitResult *fit_result);

void mc_MassFit_HighpT(
		float ptLow=6.5, float ptHigh=9.0,
		float yLow=0, float yHigh=1.6,
		int PR=0, //0=PR, 1=NP, 2=Inc.
		int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
		)
{
    TString DATE = "2DFit_No_Weight";
    gStyle->SetEndErrorSize(0);
    gSystem->mkdir("roots");
    gSystem->mkdir("roots_MC");
    gSystem->mkdir("roots_MC/Mass");
    gSystem->mkdir("figs");
    gSystem->mkdir(Form("figs/%s",DATE.Data()),kTRUE);
    gSystem->mkdir(Form("figs/%s/mc_Mass",DATE.Data()),kTRUE);

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
    TFile* f1 = new TFile("../../skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root", "read");
    //TFile* f1 = new TFile("../../skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230810.root", "read"); // mass range (2.2 - 4)
    


	// cout << "Input file: "
	// << Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",
	// fname.Data(),fEffW,fAccW,isPtW,isTnP) << endl;
	// exit(0);

	massLow=2.6;
	massHigh=3.5;
    
    //double fit_limit = 3.26;     //1st value
    double fit_limit = 3.22;     

	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f",ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;


	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<endl;
	cout << "####################################" << endl;
	RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset),WeightVar(*ws->var("weight")));
	RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );

    cout << "======= dataset ======\n";
    dataset->Print("V");
    cout << endl << endl;

    cout << "======= datasetW ======\n";
    datasetW->Print("V");
    cout << endl << endl;

    cout << "======= dsAB ======\n";
    dsAB->Print("V");
    cout << endl << endl;

	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->SetName("dsAB");
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();
	//***********************************************************************
	//****************************** MASS FIT *******************************
	//***********************************************************************

	//         The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
    double paramsupper[8] = {0.1,  1.0  ,     10,     10.,   1.0,     25.0};
    double paramslower[8] = {1e-6,  0.0,   1e-6,   1e-06,   1e-6,      -25.0};
    
	//SIGNAL: initial params
	//1st value
	//double sigma_1_init = 0.15;
    //double x_init = 0.07;
	//double alpha_1_init = 1.5;
	//double n_1_init = 1.0;
	//double f_init = 0.58;
    //double sl1_mean = 0.05, sl2_mean = 0.04, sl3_mean = 0.06;
    //double N_Jpsi_high = 800000; // 2500000
	//double N_Bkg_high = 800000;

	double sigma_1_init = 0.015;
    double x_init = 0.07;
	double alpha_1_init = 1.;
	double n_1_init = 1.;
	double f_init = 0.58;
    double sl1_mean = 0.05, sl2_mean = 0.04, sl3_mean = 0.06;
    double N_Jpsi_high = 800000; // 2500000
	double N_Bkg_high = 800000;

    //double fit_limit = 3.3;
    //if(ptLow==6.5&&ptHigh==50) {
    //    N_Jpsi_high = 5000000;
    //    sigma_1_init = 0.03;
    //    x_init = 0.2;
	//    alpha_1_init = 1.97;
	//    n_1_init = 1.48;
	//    f_init = 0.1;
    //}
    //else if(ptLow==12&&ptHigh==15) {
    //    N_Jpsi_high = 50000;
    //    sigma_1_init = 0.03;
    //    x_init = 0.2;
	//    alpha_1_init = 1.97;
	//    n_1_init = 1.48;
	//    f_init = 0.38;
    //}
    //else if(ptLow==25&&ptHigh==50) {
    //    N_Jpsi_high = 5000;
    //    sigma_1_init = 0.03;
    //    x_init = 0.2;
	//    alpha_1_init = 1.97;
	//    n_1_init = 1.48;
	//    f_init = 0.1;
    //}

	double m_lambda_init = 5;
	double mass_jpsi = pdgMass.JPsi;
    // cout << psi_2S_mass << endl;
    // exit(1);

    //SIGNAL
    RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi-0.1, pdgMass.JPsi+0.1);
    RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
    RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
    RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
    RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
    RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
    RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
    RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
    RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);

    //Set up crystal ball shapes
    RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
    RooAddPdf*  pdfMASS_Jpsi;

    //DOUBLE CRYSTAL BALL
    RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
    pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );

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

    //nMassBin = 46;
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
    RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,fit_limit), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(4));
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
    ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack), Range(massLow, fit_limit));
    // ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Range(3.4, 3.87), Normalization(fitFraction.getVal()));
    // ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Range(3.4, 3.87), Normalization(1-fitFraction.getVal()));
    ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Normalization(fitFraction.getVal()), Range(massLow, fit_limit));
    ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(massLow, fit_limit));
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
    myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);


    myPlot2_A->GetXaxis()->SetLabelSize(0);
    myPlot2_A->GetXaxis()->SetTitleSize(0);
    myPlot2_A->GetXaxis()->CenterTitle();
    //if (logY_flag == true) {
    //    myPlot2_A->SetMinimum(2*10);
    //    myPlot2_A->SetMaximum(10000000000);
    //}
    //else {
    //    myPlot2_A->GetYaxis()->SetRangeUser(0,280000);
    //}
	myPlot2_A->GetYaxis()->SetRangeUser(Ydown,Yup);
    
    myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    myPlot2_A->Draw();
    
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh),text_x,text_y,text_color,text_size);
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
	outFile = new TFile(Form("roots_MC/Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
    c_A->SaveAs(Form("figs/%s/mc_Mass/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));


    pdfMASS_Tot->Write();
    datasetMass->Write();
    outh->Write();
    // ws->Write();

    outFile->Close();
    fitMass->Print("v");
    check_convergence(fitMass);
}

void check_convergence(RooFitResult *fit_result)
{
    int hesse_code = fit_result->status();
    double edm = fit_result->edm();
    double mll = fit_result->minNll();
    //RooRealVar* par_fitresult = (RooRealVar*)fit_result->floatParsFinal().find("N_Jpsi");
    RooRealVar* fit_para = (RooRealVar*)fit_result->floatParsFinal().at(0);
    double val_ = fit_para->getVal();
    double err_ = fit_para->getError();
    double min_ = fit_para->getMin();
    double max_ = fit_para->getMax();

    cout << "\n###### Fit Convergence Check ######\n";
    cout << "Caveat: mean was fixed so it's not a matter" << endl;
    cout << "Hesse: " << hesse_code << "\t edm: " << edm << "\t mll: " << mll << endl;
    int cnt_ = 0;
    for (int idx = 0; idx < fit_result->floatParsFinal().getSize(); idx++) {
        RooRealVar *fit_para = (RooRealVar *)fit_result->floatParsFinal().at(idx);
        double val_ = fit_para->getVal();
        double err_ = fit_para->getError();
        double min_ = fit_para->getMin();
        double max_ = fit_para->getMax();
        if ((val_ - err_ > min_) && (val_ + err_ < max_)) {
            // No work is intended
        }
        else {
            cout << "\033[31m" << "[Stuck] " << fit_para->GetName() << " \033[0m // Final Value : " << val_ << " (" << min_ << " ~ " << max_ << ")" << endl << endl;
            cnt_++;
        }
    }
    if (cnt_ == 0) cout << "[Fit Converged]" << endl;
}
