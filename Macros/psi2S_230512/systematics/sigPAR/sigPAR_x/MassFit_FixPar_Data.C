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
#include "../../../../../rootFitHeaders.h"
#include "../../../../../commonUtility.h"
#include "../../../../../JpsiUtility.h"
#include "../../../../../cutsAndBin.h"
#include "../../../../../CMS_lumi_v2mass.C"
#include "../../../../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void check_convergence(RooFitResult *fit_result);
void MassFit_FixPar_Data(
    double ptLow=3, double ptHigh=4.5,
    double yLow=1.6, double yHigh=2.4,
    int cLow=0, int cHigh=200,
    int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
    )
{
  //TString DATE = "10_60";
  //TString DATE = "20_40";
  //TString DATE = "0_180";
  TString DATE;
  DATE="No_Weight";
  //if(ptLow==6.5&&ptHigh==50) DATE=Form("%i_%i",0,180);
  //else DATE=Form("%i_%i",cLow/2,cHigh/2);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/Mass",DATE.Data()),kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/Mass",DATE.Data()),kTRUE);

  TString fname;
  if (PRw==1) fname="PR";
  else if (PRw==2) fname="NP";

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TString kineLabel = getKineLabel (ptLow, ptHigh,yLow, yHigh, 0.0, cLow, cHigh);

  TFile* f1; TFile* f2; TFile* f3;
  TString kineCut;
  
  massLow=3.3;
  massHigh=4.1;
  //massLow=2.75;

//  f1 = new TFile(Form("../../skimmedFiles/v2Cut_Nom/OniaRooDataSet_isMC0_Psi2S_%s_m3.3-4.1_OS_Effw%d_Accw%d_PtW%d_TnP%d_221013_root618.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));
  f1 = new TFile(Form("../../../../../skimmedFiles/OniaRooDataSet_miniAOD_isMC0_Psi2S_cent0_200_Effw0_Accw0_PtW0_TnP0_230514.root"));
  //f1 = new TFile(Form("../../skimmedFiles/v2Cut_Nom/OniaRooDataSet_isMC0_Psi2S_%s_m3.3-4.1_OS_Effw%d_Accw%d_PtW%d_TnP%d_220808.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));
//  f1 = new TFile("../../skimmedFiles/vnCut/OniaRooDataSet_isMC0_JPsi_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS_Effw0_Accw0_PtW1_TnP1_211110.root");
//  f1 = new TFile("/Users/hwan/tools/2019/CMS/JPsi/Jpsi_v2_PbPb2018/skimmedFiles/vnCut/OniaRooDataSet_isMC0_JPsi_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS_Effw0_Accw0_PtW1_TnP1_211110.root");

  kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>3.3 && mass<4.1 && cBin>=%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);

  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  TString OS="recoQQsign==0 &&";

  //TString SglMuPt="pt1>0.5&&pt2>0.5"

  kineCut = OS+accCut+kineCut;

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow<<"-"<<cHigh<<"%"<<endl;
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
  //***********************************************************************
  //****************************** MASS FIT *******************************
  //***********************************************************************

  TFile * f_fit = new TFile(Form("../../../roots_MC/Mass/mc_MassFitResult_%s_PRw_Effw%d_Accw%d_PtW%d_TnP%d.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));
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

  if (f_higher>1.0) f_higher=1.0;
  if (f_lower<0.0) f_lower=0.0;

  double paramslower[6] = {0.,   0., 0.,  0., 0., 0.};
  double paramsupper[6] = {10., 10., 0.2, 2., 1., 25.0};

  double alpha_1_init = alpha_MC_value; double n_1_init = n_MC_value;
  double sigma_1_init = sigma_MC_value; double x_init = xA_MC_value; double f_init = f_MC_value;
  //SIGNAL FIX
  RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.Psi2S, pdgMass.Psi2S -0.05, pdgMass.Psi2S + 0.05 ) ;
  //RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init);
  RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[3], paramsupper[3]);
  RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
  //RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[0], paramsupper[0]);
  RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init);
  RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
  //RooRealVar    n_1_A("n_1_A","power order", n_1_init,n_1_init*0.9,n_1_init*1.1);
  RooRealVar    n_1_A("n_1_A","power order", n_1_init);
  RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
  //RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, f_init*0.99,f_init*1.01);
  RooRealVar   *f = new RooRealVar("f","cb fraction", f_init);

  //SIGNAL FREE
  //RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init);
  //RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init,alpha_1_init*0.9,alpha_1_init*1.1);
  //RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[1], paramsupper[1]);
  //RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);
  //Set up crystal ball shapes
  RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooAddPdf*  pdfMASS_Jpsi;
  //DOUBLE CRYSTAL BALL
  RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
  //RooGaussian* gauss = new RooGaussian("gauss","gaussian PDF",*(ws->var("mass")),mean,sigma_2_A);
  pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
  //pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal",RooArgList(*cb_1_A,*gauss), RooArgList(*f) );

  //pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
  //BACKGROUND
  //RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[5], paramsupper[5]);
  RooRealVar *sl1 = new RooRealVar("sl1","sl1", 0.0, -1., 1.); // 15<pt<50 v2==-1.2 : 0.01
  RooRealVar *sl2 = new RooRealVar("sl2","sl2", 0.0, -1., 1.);
  RooRealVar *sl3 = new RooRealVar("sl3","sl3", 0.0, -1., 1.);
  RooRealVar *sl4 = new RooRealVar("sl4","sl4", 0.0, -1., 1.);
  RooRealVar *sl5 = new RooRealVar("sl5","sl5", 0.0, -1., 1.);
  RooRealVar *sl6 = new RooRealVar("sl6","sl6", 0.0, -1., 1.);
  //RooRealVar *sl4 = new RooRealVar("sl4","sl4", .4, -5., 5.);
  //THIS IS THE BACKGROUND FUNCTION
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)*@2+@3",RooArgList(*(ws->var("mass")), m_lambda_A, *sl1, *sl2));
  //RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
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
  //RooGenericPdf *pdfMASS_bkg;
  RooChebychev *pdfMASS_bkg;
  if (ptLow<6.5) pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));
  else pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1));
  //pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1,*sl2));
  //if(ptLow==9&&ptHigh==12){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(ptLow==20&&ptHigh==50){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(ptLow==12&&ptHigh==15){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2, *sl3));}
  //if(cLow==10&&cHigh==20){pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1, *sl2));}
  //Build the model
  Double_t NBkg_limit = 2.0e+07;
  Double_t NJpsi_limit = 10.0e+06;
  Double_t NJpsi_lower = 0.;
  if (ptLow==12&&ptHigh==15)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==3.5&&ptHigh==5) {
	  NJpsi_lower = 1400;
	  NBkg_limit = 1.7e+5;
	  NJpsi_limit = 1600; }
  else if (ptLow==5&&ptHigh==6.5)  {
	   NBkg_limit = 8e+4;
	   NJpsi_limit = 1700; }
  else if (ptLow==6.5&&ptHigh==9)  {
	   NBkg_limit = 5e+5;
	   NJpsi_limit = 5e+3; }
  else if (ptLow==6.5&&ptHigh==12)  {
	   NBkg_limit = 100000;
	   NJpsi_limit = 10000; }
  else if (ptLow==15&&ptHigh==20)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==20&&ptHigh==25)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==12&&ptHigh==50)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==15&&ptHigh==50)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==25&&ptHigh==50)  {
	   NBkg_limit = 500000;
	   NJpsi_limit = 10000; }
  else if (ptLow==30&&ptHigh==50)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==20&&ptHigh==30)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==30&&ptHigh==50)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==20&&ptHigh==50)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==3.5&&cLow==60&&cHigh==100)  {
	  NBkg_limit = 50000;
	  NJpsi_limit = 2000; }
  else if (ptLow==6.5&&cLow==0&&cHigh==20)  {
	  NBkg_limit = 3e+4;
	  NJpsi_limit = 1000; }
  else if (ptLow==6.5&&cLow==20&&cHigh==40)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==6.5&&cLow==40&&cHigh==60)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (ptLow==6.5&&cLow==60&&cHigh==80) {
	  NJpsi_limit = 700;
	  NBkg_limit = 1e+4;   }
  else if (ptLow==6.5&&cLow==80&&cHigh==100) {
	  NJpsi_limit = 400;
	  NJpsi_lower = 300;
	  NBkg_limit = 4e+3;   }
  else if (ptLow==6.5&&cLow==100&&cHigh==180) {
	  NJpsi_limit = 500;
	  NBkg_limit = 3e+3;   }
  else if (ptLow==3.5&&cLow==0&&cHigh==20)  {
	  NJpsi_lower = 950;
	  NBkg_limit = 1e+6;
	  NJpsi_limit = 1500; }
  else if (ptLow==3.5&&cLow==20&&cHigh==40)  {
	  //NJpsi_lower = 600;
	  NBkg_limit = 8.3e+4;
	  NJpsi_limit = 1580; }
  else if (ptLow==3.5&&cLow==40&&cHigh==60)  {
	  //NJpsi_lower = 600;
	  NBkg_limit = 5.6e+4;
	  NJpsi_limit = 1150; }
  else if (ptLow==3.5&&cLow==60&&cHigh==80)  {
	  NJpsi_lower = 400;
	  NBkg_limit = 500000;
	  NJpsi_limit = 1000; }
  else if (ptLow==3.5&&cLow==80&&cHigh==100)  {
	  NBkg_limit = 1.5e+4;
	  NJpsi_limit = 530; }
  else if (ptLow==3.5&&cLow==100&&cHigh==180)  {
	  NBkg_limit = 7000;
	  NJpsi_limit = 400; }
  else if (cLow==0&&cHigh==40)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (cLow==20&&cHigh==60) {
	  NBkg_limit = 150000;
	  NJpsi_limit = 4000; }
  else if (cLow==40&&cHigh==80)  {
	  NBkg_limit = 500000;
	  NJpsi_limit = 10000; }
  else if (yLow==0.4&&yHigh==0.9)  {
	  NBkg_limit = 50000;
	  NJpsi_limit = 1000; }
  else if (yLow==1.6&&yHigh==2.)  {
	  NBkg_limit = 50000;
	  NJpsi_limit = 5000; }

  RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",NJpsi_lower,NJpsi_limit);
  RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,NBkg_limit);
  //RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,700000);
  //RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,1400000);
  /*RooRealVar *N_Jpsi;
	RooRealVar *N_Bkg;
  //RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,700000);
  //RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,700000);
  if(ptLow==3&&cLow==40){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,1600000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,6400000);}
  else if(ptLow==15&&ptHigh==50&&cLow==20){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,50000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,20000);}
  else if(ptLow==3&&ptHigh==6.5&&cLow==40){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,1574378);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,6283499);}
  else if(ptLow==3.5&&cLow==40){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,1350000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,4900000);}
  else if(ptLow==6.5&&cLow==40){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,450000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,1200000);}
  else if(ptLow==6.5&&ptHigh==50&&cLow==40&&cHigh==80){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,9.5981e+05);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,9.1751e+05);}
  else if(ptLow==3&&ptHigh==4.5&&cLow==60&&cHigh==100){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,900000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,4000000);}
  else if(ptLow==4.5&&ptHigh==6.5&&cLow==60&&cHigh==100){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,500000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,600000);}
  else if(ptLow==3&&ptHigh==4.5){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,2500000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,12000000);}
  else if(ptLow==3&&ptHigh==6.5){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,3014378);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,13759320);}
  else if(ptLow==4.5&&ptHigh==6.5){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,950000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,2400000);}
  else if(ptLow==15&&ptHigh==30){
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,50000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,20000);}
  else{
  N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,4500000);
  N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",0,12000000);}*/
  //pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *bkg_1order),RooArgList(*N_Jpsi,*N_Bkg));
  //pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*N_JpsiPR,*N_JpsiNP,*N_Bkg));



  RooGaussian x_constraint("x_constraint","x_constraint",*x_A,RooConst(xA_MC_value),RooConst(xA_MC_value_err));
  //RooGaussian alpha_constraint("alpha_constraint","alpha_constraint",alpha_1_A,RooConst(alpha_MC_value),RooConst(alpha_MC_value_err));
  //RooGaussian n_constraint("n_constraint","n_constraint",n_1_A,RooConst(n_MC_value),RooConst(n_MC_value_err));
  //RooGaussian f_constraint("f_constraint","f_constraint",*f,RooConst(f_MC_value),RooConst(f_MC_value_err));
  //RooGaussian sigma1_constraint("sigma1_constraint","sigma1_constraint",sigma_1_A,RooConst(sigma_MC_value),RooConst(sigma_MC_value_err));

  RooProdPdf model("model","model with constraint",RooArgSet(*pdfMASS_Jpsi,RooArgSet(x_constraint)));
  //RooProdPdf model("model","model with constraint",RooArgSet(*pdfMASS_Jpsi,RooArgSet(f_constraint,sigma1_constraint)));
  //RooProdPdf model("model","model with constraint",RooArgSet(*pdfMASS_Jpsi,RooArgSet(f_constraint,sigma1_constraint)));
  RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(model, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
  //RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
  //RooAddPdf* pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
  ws->import(*pdfMASS_Tot);

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
  //RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Constrain(RooArgSet(*f,sigma_1_A)),Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
  //RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Constrain(RooArgSet(*f,sigma_1_A)),Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
  RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
  cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,VisualizeError(*fitMass,1),FillColor(kOrange));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Tot"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Components(RooArgSet(*pdfMASS_bkg, *cb_1_A)), LineColor(44));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A, Components(RooArgSet(*pdfMASS_bkg, *cb_2_A)), LineColor(8));
  //model.plotOn(myPlot2_A,Name("model"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
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
  Yup = YMax*TMath::Power((YMax/YMin), (0.6/(1.0-0.1-0.4)));
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
  leg_B->Draw("same");

  /*drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
    if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
    else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
    drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
    drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
    drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);*/

  if(yLow==0)drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f, Cent. %d - %d%s ",ptLow, ptHigh, yHigh, cLow/2, cHigh/2, "%"),text_x,text_y,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f; Cent. %d - %d%s", ptLow, ptHigh, yLow, yHigh, cLow/2, cHigh/2, "%"), text_x,text_y,text_color,text_size);
  drawText(Form("N_{#psi(2S)} = %.f #pm %.f,  N_{Bkg} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError(),ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError())
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
  //Fixed
  drawText(Form("m_{#psi(2S)} = %.4f #pm %.4f",ws->var("m_{J/#Psi}")->getVal(),ws->var("m_{J/#Psi}")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("#alpha_{#psi(2S)} = %.4f (fixed)",ws->var("alpha_1_A")->getVal()),text_x,text_y-y_diff*3,text_color,text_size);
  drawText(Form("f_{#psi(2S)} = %.4f (fixed)",ws->var("f")->getVal()),text_x,text_y-y_diff*4,text_color,text_size);
  drawText(Form("n_{#psi(2S)} = %.4f (fixed)",ws->var("n_1_A")->getVal()),text_x,text_y-y_diff*5,text_color,text_size);
  //drawText(Form("#sigma1_{#psi(2S)} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{#psi(2S)} = %.3f (fixed)",(ws->var("sigma_1_A")->getVal())*1000,(ws->var("sigma_1_A")->getError())*1000,ws->var("x_A")->getVal()),text_x,text_y-y_diff*6,text_color,text_size);
  //Free
  //drawText(Form("m_{#psi(2S)} = %.4f #pm %.4f",ws->var("m_{J/#Psi}")->getVal(),ws->var("m_{J/#Psi}")->getError()),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("#alpha_{#psi(2S)} = %.4f #pm %.4f",ws->var("alpha_1_A")->getVal(),ws->var("alpha_1_A")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("f_{#psi(2S)} = %.4f #pm %.4f",ws->var("f")->getVal(),ws->var("f")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  //drawText(Form("n_{#psi(2S)} = %.4f #pm %.4f",ws->var("n_1_A")->getVal(),ws->var("n_1_A")->getError()),text_x,text_y-y_diff*5,text_color,text_size);
  drawText(Form("#sigma1_{#psi(2S)} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{#psi(2S)} = %.3f #pm %.3f",(ws->var("sigma_1_A")->getVal())*1000,(ws->var("sigma_1_A")->getError())*1000,ws->var("x_A")->getVal(),ws->var("x_A")->getError()),text_x,text_y-y_diff*6,text_color,text_size);
  //drawText(Form("(#sigma2/#sigma1)_{#psi(2S)} = %.3f (fixed)",ws->var("x_A")->getVal()),text_x,text_y-y_diff*7,text_color,text_size);


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
  outFile = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
  c_A->SaveAs(Form("figs/2DFit_%s/Mass/Mass_Fixed_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  pdfMASS_Tot->Write();
  pdfMASS_bkg->Write();
  datasetMass->Write();
  //ws->Write();
  outh->Write();
  fitMass->Write();
  outFile->Close();
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
