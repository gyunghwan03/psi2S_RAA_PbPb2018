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
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void MassFit_FixPar_Data_pp_wide(
    double ptLow=3, double ptHigh=4.5,
    double yLow=1.6, double yHigh=2.4,
    int PRw=1, int nGauss = 2, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
    )
{
  nCPU = 28; // Change to appropriate your cpu
  //TString DATE = "10_60";
  //TString DATE = "20_40";
  //TString DATE = "0_180";
  TString DATE;
  DATE="No_Weight";
  //if(ptLow==6.5&&ptHigh==50) DATE=Form("%i_%i",0,180);
  //else DATE=Form("%i_%i",cLow/2,cHigh/2);
  TString fname;
  if (PRw==1) fname="PR";
  else if (PRw==2) fname="NP";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots_1S_pp/%sMC/Mass",fname.Data()),kTRUE);
  gSystem->mkdir(Form("figs_1S_pp/%sMC/Mass",fname.Data()),kTRUE);


  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TString kineLabel = getKineLabelpp (ptLow, ptHigh,yLow, yHigh, 0.0);

  TFile* f1; TFile* f2; TFile* f3;
  TString kineCut;
  
  massLow=2.6;
  massHigh=3.5;
  //massLow=2.75;

  double l_cut = 10.0;

//  f1 = new TFile(Form("../../skimmedFiles/v2Cut_Nom/OniaRooDataSet_isMC0_Psi2S_%s_m3.3-4.1_OS_Effw%d_Accw%d_PtW%d_TnP%d_221013_root618.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));
  f1 = new TFile(Form("../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230209.root"));

  f2 = new TFile(Form("../../ctau_eff/roots_1S_pp/decayL/%sMC/decay_hist_%s_m2.6-3.5.root", fname.Data(), kineLabel.Data()));
  TH1D *h_Lcut = (TH1D *)f2->Get("h_Lcut");
  l_cut = h_Lcut->GetBinContent(1);
  cout << "l_cut: " << l_cut << endl;
  
  TString ctauCut;
  if(PRw==1) ctauCut = Form("ctau3D<%.5f &&",l_cut);
  else if(PRw==2) ctauCut = Form("ctau3D>%.5f &&",l_cut);

  kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5",ptLow, ptHigh, yLow, yHigh);

  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  TString OS="recoQQsign==0 &&";

  TString nan_cut = "&& !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes) ";
  kineCut = OS+accCut+ctauCut+kineCut + nan_cut;

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<endl;
  cout << "####################################" << endl;
  //RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset),WeightVar(*ws->var("weight")));
  RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset));
  RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  cout << "******** New Combined Dataset ***********" << endl;
  dsAB->SetName("dsAB");
  //ws->import(*dsAB);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
  ws->import(*dsAB);
  dsAB->Print("V");
  //***********************************************************************
  //****************************** MASS FIT *******************************
  //***********************************************************************

  //TFile * f_fit = new TFile(Form("../pp_Jpsi/roots_MC/Mass/mc_MassFitResult_%s_PRw_Effw%d_Accw%d_PtW%d_TnP%d.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));
  //RooDataSet *dataset_fit = (RooDataSet*)f_fit->Get("datasetMass");
  cout << "\n=== import mc fit results for mass pdf parameters ===\n" << endl;
  TFile * f_fit = new TFile(Form("./roots_MC/Mass_CBGauss/mc_MassFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0_Gauss%d.root", nGauss));
  RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
  ws_fit = (RooWorkspace*)f_fit->Get("workspace");

  // Left tail parameters
  Double_t alphaL_MC_value = ws_fit->var("alphaL")->getVal();
  Double_t alphaL_MC_value_err = ws_fit->var("alphaL")->getError();
  Double_t nL_MC_value = ws_fit->var("nL")->getVal();
  Double_t nL_MC_value_err = ws_fit->var("nL")->getError();
  
  // Right tail parameters
  Double_t alphaLR_MC_value = ws_fit->var("alphaLR")->getVal();
  Double_t alphaLR_MC_value_err = ws_fit->var("alphaLR")->getError();
  Double_t nLR_MC_value = ws_fit->var("nLR")->getVal();
  Double_t nLR_MC_value_err = ws_fit->var("nLR")->getError();

  Double_t massSigma1_MC_value = ws_fit->var("massSigma1")->getVal();
  Double_t massSigma1_MC_value_err = ws_fit->var("massSigma1")->getError();
  Double_t massSigma12_MC_value = ws_fit->var("massSigma12")->getVal();
  Double_t massSigma12_MC_value_err = ws_fit->var("massSigma12")->getError();
  Double_t fCB1_MC_value = ws_fit->var("fCB1")->getVal();
  Double_t fCB1_MC_value_err = ws_fit->var("fCB1")->getError();
  Double_t fCB2_MC_value = ws_fit->var("fCB2")->getVal(); 
  Double_t fCB2_MC_value_err = ws_fit->var("fCB2")->getError();
  //Double_t fGaus1_MC_value = ws_fit->var("fGaus1")->getVal();
  //Double_t fGaus1_MC_value_err = ws_fit->var("fGaus1")->getError();
  Double_t massSigma1G_MC_value = ws_fit->var("massSigma1G")->getVal();
  Double_t massSigma1G_MC_value_err = ws_fit->var("massSigma1G")->getError();
    Double_t massSigma1G2_MC_value;
    Double_t massSigma1G2_MC_value_err; 
    Double_t fGaus1_MC_value;
    Double_t fGaus1_MC_value_err;
  if(nGauss==2) { 
    massSigma1G2_MC_value = ws_fit->var("massSigma1G2")->getVal();
    massSigma1G2_MC_value_err = ws_fit->var("massSigma1G2")->getError(); 
    fGaus1_MC_value = ws_fit->var("fGaus1")->getVal();
    fGaus1_MC_value_err = ws_fit->var("fGaus1")->getError();
  }


  double sigma_index = 5;

  cout << "\n=== Mass Fit Model Definition ===\n" << endl;
  RooRealVar    massMean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ;
  RooRealVar massSigma1("massSigma1", "CB sigma", massSigma1_MC_value, 0.01, 0.05);
  //RooRealVar massSigma1("massSigma1", "CB sigma", massSigma1_MC_value, massSigma1_MC_value - (sigma_index*massSigma1_MC_value_err), massSigma1_MC_value + (sigma_index*massSigma1_MC_value_err));
  RooRealVar massSigma12("massSigma12", "ratio of sigma 1 vs 2", massSigma12_MC_value);
  //RooRealVar massSigma12("massSigma12", "ratio of sigma 1 vs 2", massSigma12_MC_value,0.1,10);
  RooFormulaVar massSigma2("massSigma2", "@0*@1", {massSigma1, massSigma12});

  RooRealVar alphaL("alphaL", "",  alphaL_MC_value);
  RooRealVar nL("nL", "nL", nL_MC_value);

  RooRealVar alphaLR("alphaLR", "",  alphaLR_MC_value);
  RooRealVar nLR("nLR", "", nLR_MC_value);
  RooFormulaVar alphaR("alphaR", "-@0*@1", {alphaL, alphaLR});
  RooFormulaVar nR("nR", "@0*@1", {nL, nLR});

	
  RooRealVar fCB1("fCB1", "frac(CB1 in total)", fCB1_MC_value);
  RooRealVar fCB2("fCB2", "frac(CB2 in total)", fCB2_MC_value);
  RooRealVar fGaus1("fGaus1", "frac(Gauss1 in total)", fGaus1_MC_value);
  cout << "HERE" << endl;
	RooCBShape CB1("CB1", "left-tail CB", *ws->var("mass"), massMean, massSigma1, alphaL, nL);
  RooCBShape CB2("CB2", "right-tail CB", *ws->var("mass"), massMean, massSigma2, alphaR, nR);
  RooAddPdf DCB("DCB", "", RooArgList(CB1, CB2), RooArgList(fCB2));
  // RooAddPdf massFitModel("massFitModel", "", RooArgList(CB1, CB2), RooArgList(fCB1));

  // Gauss
  RooRealVar massSigma1G("massSigma1G", "ratio of sigma1 vs sigmaG", massSigma1G_MC_value);
  RooFormulaVar massSigmaG("massSigmaG", "@0*@1", {massSigma1, massSigma1G});
  RooGaussian massG("massG", "core Gauss", *ws->var("mass"), massMean, massSigma1G);

  RooRealVar massSigma1G2("massSigma1G2", "ratio of sigma1 vs sigmaG", massSigma1G2_MC_value);
  RooFormulaVar massSigmaG2("massSigmaG2", "@0*@1", {massSigmaG, massSigma1G2});
  RooGaussian *massG2 = nullptr;
  cout << "HERE2" << endl;
  if(nGauss==2) {
    massG2 = new RooGaussian("massG2", "core Gauss2", *ws->var("mass"), massMean, massSigmaG2);
  }
  cout << "HERE3" << endl;
  //RooRealVar massSigma1G2("massSigma1G2", "ratio of sigma1 vs sigmaG", 1.1, 1, 100);
  //RooFormulaVar massSigmaG2("massSigmaG2", "@0*@1", {massSigmaG, massSigma1G2});
  //RooGaussian massG2("massG2", "core Gauss", *ws->var("mass"), massMean, massSigmaG2);
   
  //RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.70, 0.00, 1.00);
  //RooRealVar fCB2("fCB2", "frac(Gauss in total)", 0.70, 0.00, 1.00);
  //RooRealVar fGaus1("fGaus1", "frac(Gauss in total)", 0.70, 0.00, 1.00);
  //RooRealVar fGaus1("fGaus1", "frac(Gauss1 in total)", fGaus1_MC_value);


    
  RooAddPdf*  pdfMASS_Jpsi;
  //DOUBLE-SIDED CRYSTAL BALL (both left and right tails)
  if(nGauss==2) pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(CB1,CB2,massG,*massG2), RooArgList(fCB1,fCB2,fGaus1) );
  else if(nGauss==1) pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(CB1,CB2,massG), RooArgList(fCB1,fCB2) );
  cout << "HERE4" << endl;
  //pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal",RooArgList(*cb_1_A,*gauss), RooArgList(*f) );

  //pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
  //BACKGROUND
  //RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[5], paramsupper[5]);


  Double_t NBkg_limit = 2.0e+08;
  Double_t NJpsi_limit = 1.0e+08;
  double s1_init = 0.0; double s2_init = 0.0; double s3_init = 0.0; double s4_init = 0.0;

  if(ptLow==6.5&&ptHigh==40&&PRw==1) { NJpsi_limit=1e+8; NBkg_limit = 2e+5; s1_init = 0.; s2_init = 0.; }

  RooRealVar *sl1 = new RooRealVar("sl1","sl1", s1_init, -1., 1.); // 15<pt<50 v2==-1.2 : 0.01
  RooRealVar *sl2 = new RooRealVar("sl2","sl2", s2_init, -1., 1.);
  RooRealVar *sl3 = new RooRealVar("sl3","sl3", s3_init, -1., 1.);
  RooRealVar *sl4 = new RooRealVar("sl4","sl4", s4_init, -1., 1.);

  //THIS IS THE BACKGROUND FUNCTION
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)*@2+@3",RooArgList(*(ws->var("mass")), m_lambda_A, *sl1, *sl2));
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
  RooChebychev *pdfMASS_bkg;
  //RooExponential *pdfMASS_bkg;
  //if(ptLow==3){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));}
  //if(ptLow!=3){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));}
  //  if(ptLow=15){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  /*if(ptLow<=6.5){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1));}
    else if(ptLow<=6.5&&ptHigh==50){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));}
    else pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));*/
  //pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1));
  //*sl1, *sl2, *sl3, *sl4, *sl5, *sl6
  //pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList());
  //pdfMASS_bkg = new RooExponential("pdfMASS_bkg","Background",*(ws->var("mass")),*sl1);
  if (ptLow<6.5) pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));
  else if(ptLow==6.5&&ptHigh==40&&PRw==1) pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));
  else pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1));

  //pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1,*sl2));
  //if(ptLow==9&&ptHigh==12){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(ptLow==20&&ptHigh==50){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(ptLow==12&&ptHigh==15){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(cLow==10&&cHigh==20){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));}
  
  
  //Build the model
  RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,NJpsi_limit);
  RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,NBkg_limit);

  //RooGaussian x_constraint("x_constraint","x_constraint",*x_A,RooConst(xA_MC_value),RooConst(xA_MC_value_err));

  //RooProdPdf model("model","model with constraint",RooArgSet(*pdfMASS_Jpsi,RooArgSet(f_constraint,sigma1_constraint)));
  //RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
  RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
  cout << "HERE5" << endl;
  ws->import(*pdfMASS_Tot);

  cout << "HERE6" << endl;
  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.25, 0.98, 1.0);
  pad_A_1->SetTicks(1,1);
  pad_A_1->Draw(); pad_A_1->cd();
  RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  myPlot_A->SetTitle("");
  ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));
  pad_A_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
  dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));

  bool isWeighted = ws->data("dsAB")->isWeighted();
  //bool isWeighted = true;
  cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
  RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));


  cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
//  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,VisualizeError(*fitMass,1),FillColor(kOrange));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Tot"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("CB1"), Components(RooArgSet(*pdfMASS_bkg, CB1)), LineColor(kMagenta+2),LineWidth(2),LineStyle(kDashed));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("CB2"), Components(RooArgSet(*pdfMASS_bkg, CB2)), LineColor(kMagenta-2),LineWidth(2),LineStyle(kDashed));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("massG"), Components(RooArgSet(*pdfMASS_bkg, massG)), LineColor(kGreen+2),LineWidth(2),LineStyle(kDashed));
  if(nGauss==2) ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("massG2"), Components(RooArgSet(*pdfMASS_bkg, massG2)), LineColor(kGreen-2),LineWidth(2),LineStyle(kDashed));
  //model.plotOn(myPlot2_A,Name("model"), LineColor(kBlack));
  dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));
  //model.plotOn(myPlot2_A,Name("model_Bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
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
  myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
  myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2_A->Draw();

  TLegend* leg_B = new TLegend(text_x+0.5,text_y-0.2,text_x+0.7,text_y); leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot2_A->findObject("dataOS"),"Data","pe");
  leg_B->AddEntry(myPlot2_A->findObject("pdfMASS_Tot"),"Total","l");
  leg_B->AddEntry(myPlot2_A->findObject("pdfMASS_bkg"),"Background","l");
  leg_B->AddEntry(myPlot2_A->findObject("CB1"),"CB1","l");
  leg_B->AddEntry(myPlot2_A->findObject("CB2"),"CB2","l");
  leg_B->AddEntry(myPlot2_A->findObject("massG"),"massG","l");
  if(nGauss==2) leg_B->AddEntry(myPlot2_A->findObject("massG2"),"massG2","l");
  leg_B->Draw("same");  
  /*drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
    if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
    else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
    drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
    drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);*/
  if(yLow==0)drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f",ptLow, ptHigh, yHigh),text_x,text_y,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f;", ptLow, ptHigh, yLow, yHigh), text_x,text_y,text_color,text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f,  N_{Bkg} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError(),ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError())
            ,text_x,text_y-y_diff*1,text_color,text_size);
   // drawText(Form("#alpha = %.4f (fixed)  f = %.4f (fixed)  n_{1} = %.4f (fixed)", ws->var("alpha_1_A")->getVal(), fitFraction.getVal(), fitN_1.getVal()),text_x,text_y-y_diff*2,text_color,text_size);
   // drawText(Form("#sigma_{1} = %.4f #pm %.4f   #sigma_{2} / #sigma_{1} = %.4f (fixed)",fitSigma_1.getVal(),fitSigma_1.getError(),fitX.getVal()),text_x,text_y-y_diff*3,text_color,text_size);

    //drawText(Form("sl1 = %.4f #pm %.4f "
    //            ,fitSl1.getVal(),fitSl1.getError())
	//		,text_x,text_y-y_diff*4,text_color,text_size);
  //if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  //else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  //else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f,  %.1f < v_{2} < %.1f",yLow, yHigh, v2, v2+0.3), text_x,text_y-y_diff,text_color,text_size);
  //drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("N_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("N_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  drawText(Form("m_{J/#psi} = %.4f #pm %.4f",ws->var("m_{J/#Psi}")->getVal(),ws->var("m_{J/#Psi}")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("#alpha_{L} = %.4f (fixed), #alpha_{R} = %.4f (fixed)",ws->var("alpha_L_A")->getVal(),ws->var("alpha_R_A")->getVal()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("f_{J/#psi} = %.4f (fixed)",ws->var("f")->getVal()),text_x,text_y-y_diff*4,text_color,text_size);
  //drawText(Form("n_{L} = %.4f (fixed), n_{R} = %.4f (fixed)",ws->var("n_L_A")->getVal(),ws->var("n_R_A")->getVal()),text_x,text_y-y_diff*5,text_color,text_size);
  //drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)",(ws->var("sigma_1_A")->getVal())*1000,(ws->var("sigma_1_A")->getError())*1000,ws->var("x_A")->getVal()),text_x,text_y-y_diff*6,text_color,text_size);
  if(PRw==1) drawText(Form("l_{J/#psi} < %.4f", l_cut),text_x,text_y-y_diff*7,text_color,text_size);
  else if(PRw==2) drawText(Form("l_{J/#psi} > %.4f", l_cut),text_x,text_y-y_diff*7,text_color,text_size);
  //drawText(Form("(#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)",ws->var("x_A")->getVal()),text_x,text_y-y_diff*7,text_color,text_size);


  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.001, 0.98, 0.32);
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
  RooHist* hpull_A = frameTMP->pullHist("dataOS","pdfMASS_Tot", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot* pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame_A->addPlotable(hpull_A,"P") ;
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_A->GetYaxis()->SetTitle("Pull") ;
  pullFrame_A->GetYaxis()->SetTitleSize(0.08) ;
  pullFrame_A->GetYaxis()->SetLabelSize(0.08) ;
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_A->GetXaxis()->SetLabelSize(0.08) ;
  pullFrame_A->GetXaxis()->SetTitleSize(0.08) ;
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw() ;

  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
  outh->GetXaxis()->SetBinLabel(1,"Jpsi");

  float temp1 = ws->var("N_Jpsi")->getVal();
  float temp1err=ws->var("N_Jpsi")->getError();

  outh->SetBinContent(1,temp1);
  outh->SetBinError(1,temp1err);

  fitMass->Print();
  Double_t theNLL = fitMass->minNll();
  cout << " *** NLL : " << std::setprecision(15) << theNLL << endl;

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

  TFile* outFile;
  outFile = new TFile(Form("roots_1S_pp/%sMC/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Gauss%d.root", fname.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP, nGauss),"recreate");
  c_A->SaveAs(Form("figs_1S_pp/%sMC/Mass_Fixed_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Gauss%d.pdf", fname.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP, nGauss));
  pdfMASS_Tot->Write();
  pdfMASS_bkg->Write();
  datasetMass->Write();
  //ws->Write();
  outh->Write();
  fitMass->Write();
  outFile->Close();
}
