#include <iostream>
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
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
using namespace std;
using namespace RooFit;


void check_convergence(RooFitResult *fit_result);

void mc_MassFit_DCBGauss(
		double ptLow=6.5, double ptHigh=9.0,
		double yLow=0, double yHigh=1.6,
		int PRw=1, //1=PR, 2=NP, 0=Inc.
		int PR=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
		)
{
    TString DATE = "2DFit_No_Weight";
	TString fname;
	if (PRw==1) fname="PR";
	else if (PRw==2) fname="NP";

    nCPU=28;

    gStyle->SetEndErrorSize(0);
    gSystem->mkdir("roots");
    gSystem->mkdir("roots_MC");
    gSystem->mkdir("roots_MC/DCBGauss");
    gSystem->mkdir("figs");
    gSystem->mkdir("figs/DCBGauss",kTRUE);
    gSystem->mkdir("figs/DCBGauss/mc_Mass",kTRUE);

	TString bCont;
	if(PR==0) bCont="Prompt";
	else if(PR==1) bCont="NonPrompt";
	else if(PR==2) bCont="Inclusive";

	
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
    TFile* f1 = new TFile("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root", "read");
    //TFile* f1;
    //if(PRw==1) f1= new TFile("../../skimmedFiles/OniaRooDataSet_miniAOD_isMC1_JPsi_Prompt_cent0_200_Effw0_Accw0_PtW0_TnP0_240530.root", "read");
    //else if(PRw==2) f1 = new TFile("../../skimmedFiles/OniaRooDataSet_miniAOD_isMC1_JPsi_NonPrompt_cent0_200_Effw0_Accw0_PtW0_TnP0_250421.root","read");


	// cout << "Input file: "
	// << Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",
	// fname.Data(),fEffW,fAccW,isPtW,isTnP) << endl;
	// exit(0);

	massLow=2.6;
	massHigh=3.5;

	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5",ptLow, ptHigh, yLow, yHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;


	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<< endl;
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

    double sigma_low  = 0.001;
    double sigma_high = 0.08;

    double x_low      = 1.;
    double x_high     = 5.0;

    double alpha_low  = 1.0;
    double alpha_high = 6.0;

    double n_low      = 1.0;
    double n_high     = 8.0;

    double f_low      = 0.0;
    double f_high     = 1.0;
	//         The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
    
	//SIGNAL: initial params
    double alpha_1_init = 2.9801; double n_1_init = 1.7473;
    double sigma_1_init = 0.0427; double f_init = 0.0726;
    double x_init = 3.5297;
    double N_Jpsi_high = 1e+8;
    double fit_limit = 3.36;
    double paramsupper[8] = {sigma_high, x_high, alpha_high, n_high, f_high,  25.0};
    double paramslower[8] = {sigma_low,  x_low,  alpha_low,   n_low,  f_low, -25.0};

	double m_lambda_init = 5;
	double psi_2S_mass = pdgMass.Psi2S;
    // cout << psi_2S_mass << endl;
    // exit(1);

    // ----------------------------
    // SIGNAL: mass model
    // ----------------------------

    // 공통 mean
    RooRealVar mean("m_{J/#Psi}",
                    "mean of the signal gaussian mass PDF",
                    pdgMass.JPsi,
                    pdgMass.JPsi-0.1,
                    pdgMass.JPsi+0.1);

    // sigma1, sigma2 (sigma2 = x_A * sigma1)
    RooRealVar   *x_A      = new RooRealVar("x_A","sigma ratio",
                                            x_init, paramslower[1], paramsupper[1]);
    RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",
                            sigma_1_init, paramslower[0], paramsupper[0]);
    RooFormulaVar sigma_2_A("sigma_2_A","@0*@1", RooArgList(sigma_1_A, *x_A));

    // tails (alpha, n 공유)
    RooRealVar    alpha_1_A("alpha_1_A","tail shift",
                            alpha_1_init, paramslower[2], paramsupper[2]);
    RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0", RooArgList(alpha_1_A));

    RooRealVar    n_1_A("n_1_A","power order",
                        n_1_init , paramslower[3], paramsupper[3]);
    RooFormulaVar n_2_A("n_2_A","1.0*@0", RooArgList(n_1_A));

    // ----------------------------
    // 1) Double Crystal Ball (CB1 + CB2)
    // ----------------------------

    // CB1
    RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "Crystal Ball 1",
                                        *(ws->var("mass")), mean,
                                        sigma_1_A, alpha_1_A, n_1_A);

    // CB2 (sigma2_A 사용)
    RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "Crystal Ball 2",
                                        *(ws->var("mass")), mean,
                                        sigma_2_A, alpha_2_A, n_2_A);

    // DoubleCB 안에서 CB2의 비율
    RooRealVar *f_CB2 = new RooRealVar("f_CB2","fraction of CB2 in DoubleCB",
                                       f_init, paramslower[4], paramsupper[4]);

    // DoubleCB = f_CB2 * CB2 + (1-f_CB2) * CB1
    RooAddPdf* pdfDoubleCB =
      new RooAddPdf("pdfDoubleCB","Double Crystal Ball",
                    RooArgList(*cb_2_A, *cb_1_A),  // 첫 번째가 f_CB2에 해당
                    RooArgList(*f_CB2));

    // ----------------------------
    // 2) CB sigma와 연동된 Gaussian
    //     sigmaG = kG * sigma_1_A
    // ----------------------------
    RooRealVar   kG("kG","Gauss / CB sigma ratio",
                8, 5, 10.);
    RooFormulaVar sigmaG("sigmaG","@0*@1",
                     RooArgList(sigma_1_A, kG));

    RooGaussian sigG("sigG","Gaussian core linked to CB sigma",
                 *(ws->var("mass")), mean, sigmaG);


    // ----------------------------
    // 3) Signal = DoubleCB + Gauss
    //     pdfMASS_Jpsi = (1-fG)*DoubleCB + fG*Gauss
    // ----------------------------
    RooRealVar *fG = new RooRealVar("fG","Gaussian fraction in signal",
                                .03, 1e-3, 1.);

    RooAddPdf* pdfMASS_Jpsi =
      new RooAddPdf("pdfMASS_Jpsi","Signal (DoubleCB + Gauss)",
                    RooArgList(sigG, *pdfDoubleCB), // 첫 번째가 fG
                    RooArgList(*fG));

    // ----------------------------
    // 4) Extended yield
    // ----------------------------
    RooRealVar *N_Jpsi = new RooRealVar("N_Jpsi","inclusive Jpsi signals",
                                        0, N_Jpsi_high);

    // (A) RooExtendPdf 버전 (추천)
    RooExtendPdf* pdfMASS_Tot =
      new RooExtendPdf("pdfMASS_Tot","Jpsi only",
                       *pdfMASS_Jpsi, *N_Jpsi);

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
    RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,fit_limit), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
    //RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
    cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
	fitMass->Print("V");
    
    // Check and get fitted parameters
    const RooArgList & fitParams = fitMass->floatParsFinal();
    
    auto & fitN_Jpsi = (RooRealVar &) fitParams[0];
    auto & fitAlpha = (RooRealVar &) fitParams[1];
    auto & fitFraction = (RooRealVar &) fitParams[2];
    auto & fitFraction2 = (RooRealVar &) fitParams[3];
    auto & fitMean = (RooRealVar &) fitParams[4];
    auto & fitN_1 = (RooRealVar &) fitParams[5];
    auto & fitSigmaG = (RooRealVar &) fitParams[6];
    auto & fitSigma_1 = (RooRealVar &) fitParams[7];
    auto & fitX = (RooRealVar &) fitParams[8];
    // cout << fitN_Jpsi.getVal() << endl;
    
    
    for ( int i = 0; i < fitParams.getSize(); ++i)
    {
      auto & fitPar = (RooRealVar &) fitParams[i];
      std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
    }
    
    
    double f_factor = (double) fitFraction.getVal();
    ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack), Range(2.6, fit_limit));
    // ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Range(3.4, 3.87), Normalization(fitFraction.getVal()));
    // ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Range(3.4, 3.87), Normalization(1-fitFraction.getVal()));
    ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Normalization(fitFraction.getVal()), Range(2.6, fit_limit));
    ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(2.6, fit_limit));
    ws->pdf("sigG")->plotOn(myPlot2_A,Name("sigG"), LineColor(kMagenta+2), Normalization(1-f_factor), Range(2.6, fit_limit));
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

    TLegend* leg_B = new TLegend(text_x+0.5,text_y-0.2,text_x+0.7,text_y); leg_B->SetTextSize(text_size);
    leg_B->SetTextFont(43);
    leg_B->SetBorderSize(0);
    leg_B->AddEntry(myPlot2_A->findObject("pdfMASS_tot"),"Total Fit","l");
    leg_B->AddEntry(myPlot2_A->findObject("cball_1_A"),"CB1","l");
    leg_B->AddEntry(myPlot2_A->findObject("cball_2_A"),"CB2","l");
    leg_B->AddEntry(myPlot2_A->findObject("sigG"),"Gauss","l");
    leg_B->Draw();  
    
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh),text_x,text_y,text_color,text_size);
    if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
    else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
    drawText(Form("N_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
    // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
    drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(),ws->var("alpha_1_A")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    drawText(Form("f = %.4f #pm %.4f", ws->var("f_CB2")->getVal(), ws->var("f_CB2")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
    drawText(Form("n_{1} = %.4f #pm %.4f",ws->var("n_1_A")->getVal(),ws->var("n_1_A")->getError()),text_x,text_y-y_diff*5,text_color,text_size);
    drawText(Form("#sigma_{1} = %.4f #pm %.4f",ws->var("sigma_1_A")->getVal(),ws->var("sigma_1_A")->getError()),text_x,text_y-y_diff*6,text_color,text_size);
    drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f",ws->var("x_A")->getVal(),ws->var("x_A")->getError()),text_x,text_y-y_diff*7,text_color,text_size);
    drawText(Form("#sigma_{Gauss} / #sigma_{1} = %.4f #pm %.4f",ws->var("kG")->getVal(),ws->var("kG")->getError()),text_x,text_y-y_diff*8,text_color,text_size);
    

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
	outFile = new TFile(Form("roots_MC/DCBGauss/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
    c_A->SaveAs(Form("figs/DCBGauss/mc_Mass/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
//	outFile = new TFile(Form("jpsi_0to1p6/roots/MC/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
//    c_A->SaveAs(Form("jpsi_0to1p6/figs/MC/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
//	outFile = new TFile(Form("jpsi_1p6to2p4/roots/MC/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
//    c_A->SaveAs(Form("jpsi_1p6to2p4/figs/MC/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));


    pdfMASS_Tot->Write();
    datasetMass->Write();
    outh->Write();
    // ws->Write();

    outFile->Close();
    fitMass->Print("v");
    Double_t theNLL = fitMass->minNll();
    cout << " *** NLL : " << theNLL << endl;
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
