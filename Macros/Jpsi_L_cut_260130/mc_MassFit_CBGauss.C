#include <iostream>
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
using namespace std;
using namespace RooFit;

void mc_MassFit_CBGauss(
		double ptLow=6.5, double ptHigh=9.0,
		double yLow=0, double yHigh=1.6,
        int nGauss = 1,
		int PR=0, //0=PR, 1=NP, 2=Inc.
		int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
		)
{
    nCPU=28;
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
	TFile* f1 = new TFile("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root", "read");

	// cout << "Input file: "
	// << Form("../data/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210111.root",
	// fname.Data(),fEffW,fAccW,isPtW,isTnP) << endl;
	// exit(0);

	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5",ptLow, ptHigh, yLow, yHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;

	massLow=2.6;
	massHigh=3.5;

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


    ws->var("mass")->setRange(2.6, massHigh);
    ws->var("mass")->setRange("massRange", 2.6, massHigh);

    // --- build mass signal model ---
    RooRealVar massMean("massMean", "", 3.096, 3.09, 3.102);
    // massMean.setConstant();
    RooRealVar massSigma1("massSigma1", "CB sigma", 0.01, 0.001, 1);
    RooRealVar massSigma12("massSigma12", "ratio of sigma 1 vs 2", 1.5, 1, 20);
    RooFormulaVar massSigma2("massSigma2", "@0*@1", {massSigma1, massSigma12});

    RooRealVar alphaL("alphaL", "", 1.5, 0.1, 3.0);
    RooRealVar nL("nL", "nL", 1.5, 1, 5.0);

	  // RooRealVar alphaR("alphaR", "", 1.0, 0, 10.0);
    // RooRealVar nR("nR", "", 2.0, 1, 10.0);
    
    RooRealVar alphaLR("alphaLR", "", 1.2, 1, 5.0);
    RooRealVar nLR("nLR", "", 1.3, 1, 20.0);
    RooFormulaVar alphaR("alphaR", "-@0*@1", {alphaL, alphaLR});
    RooFormulaVar nR("nR", "@0*@1", {nL, nLR});

    // RooFormulaVar alphaR("alphaR", "@0", {alphaL});
    // RooFormulaVar nR("nR", "@0", {nL});

    // RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2)
    // RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);

    // RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.70, 0.00, 1.00);
	
	RooCBShape CB1("CB1", "left-tail CB", *ws->var("mass"), massMean, massSigma1, alphaL, nL);
    RooCBShape CB2("CB2", "right-tail CB", *ws->var("mass"), massMean, massSigma2, alphaR, nR);
    // RooAddPdf massFitModel("massFitModel", "", RooArgList(CB1, CB2), RooArgList(fCB1));

    // Gauss
    RooRealVar massSigma1G("massSigma1G", "ratio of sigma1 vs sigmaG", 1.1, 1e-6, 10);
    RooFormulaVar massSigmaG("massSigmaG", "@0*@1", {massSigma1, massSigma1G});
    RooGaussian massG("massG", "core Gauss", *ws->var("mass"), massMean, massSigma1G);

    RooRealVar massSigma1G2("massSigma1G2", "ratio of sigma1 vs sigmaG", 1.1, 1e-6, 10);
    RooFormulaVar massSigmaG2("massSigmaG2", "@0*@1", {massSigmaG, massSigma1G2});
    RooGaussian massG2("massG2", "core Gauss", *ws->var("mass"), massMean, massSigmaG2);
  
    // RooCrystalBall DCB("DCB", "", *ws->var("mass"), massMean, massSigma1, alphaL, nL, alphaR, nR);

    RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.70, 1e-6, 1.00);
    RooRealVar fCB2("fCB2", "frac(Gauss in total)", 0.70, 1e-5, 1.00);
    RooRealVar fGaus1("fGaus1", "frac(Gauss in total)", 0.70, 1e-6, 1.00);

    RooAddPdf* massFitModel;
    if(nGauss==2) {
        massFitModel = new RooAddPdf("massFitModel", "", RooArgList(CB1, CB2, massG, massG2), RooArgList(fCB1,fCB2,fGaus1), true); // Recursive fracton
    }
    else if (nGauss==1) {
        massFitModel = new RooAddPdf("massFitModel", "", RooArgList(CB1, CB2, massG), RooArgList(fCB1,fCB2), true); // Recursive fracton
    }
    //RooAddPdf DCB("DCB", "", RooArgList(CB1, CB2), RooArgList(fCB1));

    // --- model ---
    RooRealVar nSigMass("nSigMass", "mass signal yield", 1e4, 1, 1e8);

    double N_Jpsi_high = 1e+8;
    double fit_limit = massHigh;
    //Build the model
    RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",100, N_Jpsi_high);
    RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi",RooArgList(*massFitModel),RooArgList(*N_Jpsi));
    ws->import(*pdfMASS_Tot);

    //nMassBin = 25;
    TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
    c_A->cd();
    
    // Create pad_A_1
    TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.227, 0.98, 1.0);
    pad_A_1->SetTicks(1,1);
    pad_A_1->Draw(); 
    pad_A_1->cd();
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
    RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), Range(2.6, fit_limit), NumCPU(nCPU));
    cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
	fitMass->Print("V");
    
    
    // Check and get fitted parameters
    const RooArgList & fitParams = fitMass->floatParsFinal();
    
    // Print all fit parameters for debugging
    for ( int i = 0; i < fitParams.getSize(); ++i)
    {
      auto & fitPar = (RooRealVar &) fitParams[i];
      std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
    }
    
    // Removed hardcoded array access that caused segmentation fault
    // The number of parameters depends on nGauss value
    //ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack), Range(2.6, fit_limit));
    // ws->pdf("cball_1_A")->plotOn(myPlot2_A,Name("cball_1_A"), LineColor(kBlue+2), Range(2.6, 3.87), Normalization(fitFraction.getVal()));
    // ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Range(2.6, 3.87), Normalization(1-fitFraction.getVal()));
    //ws->pdf("DCB")->plotOn(myPlot2_A,Name("DCB"), LineColor(kBlue+2), Normalization(fitFraction.getVal()), Range(2.6, fit_limit));
    //ws->pdf("massG")->plotOn(myPlot2_A,Name("massG"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(2.6, fit_limit));
    //ws->pdf("massG2")->plotOn(myPlot2_A,Name("massG2"), LineColor(kRed+2), Normalization(1-fitFraction.getVal()), Range(2.6, fit_limit));
    //ws->pdf("cball_2_A")->plotOn(myPlot2_A,Name("cball_2_A"), LineColor(kGreen+2), Normalization(1-fitFraction.getVal()), Range(2.6, fit_limit));
    //ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));

    ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack), Range(2.6, fit_limit));
    massFitModel->plotOn(myPlot2_A, Range("massRange"), Name("model"));
    massFitModel->plotOn(myPlot2_A, Range("massRange"), Components("CB1, CB2"), Name("DCB"), LineStyle(kDotted), LineColor(kRed));
    massFitModel->plotOn(myPlot2_A, Range("massRange"), Components("massG"), Name("G"), LineStyle(kDotted), LineColor(kGreen));
    massFitModel->plotOn(myPlot2_A, Range("massRange"), Components("massG2"), Name("G2"), LineStyle(kDotted), LineColor(kOrange));

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
    
    if(ws->var("alphaL")) drawText(Form("alpha_{L} = %.4f #pm %.4f", ws->var("alphaL")->getVal(),ws->var("alphaL")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    if(ws->var("nL")) drawText(Form("n_{L} = %.4f #pm %.4f", ws->var("nL")->getVal(),ws->var("nL")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
    if(ws->var("alphaLR")) drawText(Form("alpha_{L/R} = %.4f #pm %.4f", ws->var("alphaLR")->getVal(),ws->var("alphaLR")->getError()),text_x,text_y-y_diff*5,text_color,text_size);
    if(ws->var("nLR")) drawText(Form("n_{L/R} = %.4f #pm %.4f", ws->var("nLR")->getVal(),ws->var("nLR")->getError()),text_x,text_y-y_diff*6,text_color,text_size);
    if(ws->var("fCB1")) drawText(Form("f_{CB1} = %.4f #pm %.4f", ws->var("fCB1")->getVal(),ws->var("fCB1")->getError()),text_x,text_y-y_diff*7,text_color,text_size);
    if(ws->var("fCB2")) drawText(Form("f_{CB2} = %.4f #pm %.4f", ws->var("fCB2")->getVal(),ws->var("fCB2")->getError()),text_x,text_y-y_diff*8,text_color,text_size);
    if(ws->var("fGaus1")) drawText(Form("f_{G1} = %.4f #pm %.4f", ws->var("fGaus1")->getVal(),ws->var("fGaus1")->getError()),text_x,text_y-y_diff*9,text_color,text_size);
    if(ws->var("massSigma1")) drawText(Form("#sigma_{1} = %.4f #pm %.4f", ws->var("massSigma1")->getVal(),ws->var("massSigma1")->getError()),text_x,text_y-y_diff*10,text_color,text_size);
    if(ws->var("massSigma12")) drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f", ws->var("massSigma12")->getVal(),ws->var("massSigma12")->getError()),text_x,text_y-y_diff*11,text_color,text_size);
    if(ws->var("massSigma1G")) drawText(Form("#sigma_{G} / #sigma_{1} = %.4f #pm %.4f", ws->var("massSigma1G")->getVal(),ws->var("massSigma1G")->getError()),text_x,text_y-y_diff*12,text_color,text_size);
    if(ws->var("massSigma1G2")) drawText(Form("#sigma_{G2} / #sigma_{G} = %.4f #pm %.4f", ws->var("massSigma1G2")->getVal(),ws->var("massSigma1G2")->getError()),text_x,text_y-y_diff*13,text_color,text_size);
    //drawText(Form("s"))
    // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
    //drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(),ws->var("alpha_1_A")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    //drawText(Form("f = %.4f #pm %.4f",fitFraction.getVal(),fitFraction.getError()),text_x,text_y-y_diff*4,text_color,text_size);
    //drawText(Form("n_{1} = %.4f #pm %.4f",fitN_1.getVal(),fitN_1.getError()),text_x,text_y-y_diff*5,text_color,text_size);
    //drawText(Form("#sigma_{1} = %.4f #pm %.4f",fitSigma_1.getVal(),fitSigma_1.getError()),text_x,text_y-y_diff*6,text_color,text_size);
    //drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f",fitX.getVal(),fitX.getError()),text_x,text_y-y_diff*7,text_color,text_size);

    // Create frameTMP BEFORE adding legend (important!)
    RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");

    // Draw Legend
    TLegend* leg_sig = new TLegend(0.62, 0.55, 0.88, 0.85);
    leg_sig->SetTextSize(0.035);  // Use NDC units (0-1), not pixels
    leg_sig->SetBorderSize(0);
    leg_sig->SetFillStyle(0);
    
    TObject* dataObj = myPlot2_A->findObject("dataOS");
    TObject* modelObj = myPlot2_A->findObject("pdfMASS_tot");
    TObject* dcbObj = myPlot2_A->findObject("DCB");
    TObject* gObj = myPlot2_A->findObject("G");
    TObject* g2Obj = myPlot2_A->findObject("G2");
    
    if(dataObj) leg_sig->AddEntry(dataObj, "Data", "pe");
    if(modelObj) leg_sig->AddEntry(modelObj, "Total Fit", "l");
    if(dcbObj) leg_sig->AddEntry(dcbObj, "CBs", "l");
    if(gObj) leg_sig->AddEntry(gObj, "Gauss", "l");
    if(g2Obj) leg_sig->AddEntry(g2Obj, "Gauss2", "l");
    
    leg_sig->Draw("SAME");
    pad_A_1->GetListOfPrimitives()->Add(leg_sig);
    pad_A_1->Modified();
    pad_A_1->Update();

    // Create and draw pad_A_2
    c_A->cd();
    TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
    pad_A_2->Draw();
    pad_A_2->cd();
    pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
    pad_A_2->SetBottomMargin(0.67);
    pad_A_2->SetBottomMargin(0.4);
    pad_A_2->SetFillStyle(4000);
    pad_A_2->SetFrameFillStyle(4000);
    pad_A_2->SetTicks(1,1);
    
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

    c_A->cd();
    pad_A_1->cd();
    pad_A_1->Modified();
    pad_A_1->Update();
    pad_A_2->cd();
    pad_A_2->Modified();
    pad_A_2->Update();
    c_A->cd();
    c_A->Modified();
    c_A->Update();
    
    TString kineLabel = getKineLabelpp (ptLow, ptHigh,yLow, yHigh, 0.0);
    
    TFile* outFile;
    outFile = new TFile(Form("roots_MC/Mass_CBGauss/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Gauss%d.root", kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, nGauss),"recreate");
    c_A->SaveAs(Form("figs/2DFit_%s/mc_Mass_CBGauss/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Gauss%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, nGauss));
    

    //cout << "\n\n\nHello World\n\n\n" << endl;
    //exit(1);

    pdfMASS_Tot->Write();
    datasetMass->Write();
    outh->Write();
    ws->Write();
    fitMass->Write();
    outFile->Close();
    fitMass->Print("v");
    
}
