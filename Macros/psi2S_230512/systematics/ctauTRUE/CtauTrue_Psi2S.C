#include "../../../../rootFitHeaders.h"
#include "../../../../commonUtility.h"
#include "../../../../JpsiUtility.h"
#include "../../../../cutsAndBin.h"
#include "../../../../CMS_lumi_v2mass.C"
#include "../../../../tdrstyle.C"
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
#include "TStopwatch.h"

using namespace std;
using namespace RooFit;

void CtauTrue_Psi2S(
    float ptLow=4.5, float ptHigh=6.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    float muPtCut=0.0,
    int PR=2, //0=PR, 1=NP, 2=Inc.
    float ctauCut=0.08
    )
{
  TStopwatch *t = new TStopwatch;
  t->Start();
  
  TString DATE="No_Weight";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/CtauTrue",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauTrue",DATE.Data()), kTRUE);

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
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* f1; TFile* f2; 
  TString kineCut;
  TString kineCutMC;
  TString SigCut;
  TString BkgCut;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);
  cout << kineLabel << endl;

  //f1 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco_3_6p5_201005.root");
  //f2 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco_6p5_50_201005.root");
  //f1 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco.root");
  //f2 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco.root");
  //f1 = new TFile("../../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco.root");
  f1 = new TFile("../../../../skimmedFiles/OniaRooDataSet_Psi2S_miniAOD_GENONLY_NonPrompt_20230517.root","read");
  
  kineCutMC = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>3.3 && mass<4.1" ,ptLow, ptHigh, yLow, yHigh);
  
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  //kineCutMC = accCut+kineCutMC;
  kineCutMC = kineCutMC;

  RooDataSet *datasetMC;
  datasetMC = (RooDataSet*)f1->Get("dataset");
  //if(ptHigh <= 6.5) {datasetMC = (RooDataSet*)f1->Get("dataset");}
  //else {datasetMC = (RooDataSet*)f2->Get("dataset");}

  RooWorkspace *ws = new RooWorkspace("workspace");
  RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
  wsmc->import(*datasetMC);

  //wsmc->import(*datasetMC2);
  RooCategory tp("tp","tp");
  tp.defineType("C");
  //tp.defineType("D");
  //RooDataSet *datasetMCW = new RooDataSet("datasetW","A sample",RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y")), *(wsmc->var("pt1")), *(wsmc->var("pt2")), *(wsmc->var("eta1")), *(wsmc->var("eta2"))), Import(*datasetMC) );
  RooDataSet *datasetMCW = new RooDataSet("datasetW","A sample",RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y")), *(wsmc->var("weight")), *(wsmc->var("pt1")), *(wsmc->var("pt2")), *(wsmc->var("eta1")), *(wsmc->var("eta2"))), Import(*datasetMC), WeightVar(*wsmc->var("weight")));
  wsmc->import(*datasetMCW);
  RooDataSet *reducedDS_MC = (RooDataSet*)datasetMCW->reduce(RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))), kineCutMC.Data() );
  //RooDataSet *reducedDS_MC = (RooDataSet*)datasetMC->reduce(RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))), kineCutMC.Data() );
  //RooDataSet* reducedDS_MC = new RooDataSet("reducedDS_MC","reducedDS_MC",RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))),Index(tp),Import("C",*datasetMC1),Import("D",*datasetMC2));
  reducedDS_MC->SetName("reducedDS_MC");
  reducedDS_MC->Print();
  ws->import(*reducedDS_MC);
  //ws->var("ctau3Dtrue")->setRange(0, 5);

  //Double_t ctauTrueMax; Double_t ctauTrueMin;
  //ws->data("reducedDS_MC")->getRange(*ws->var("ctau3Dtrue"), ctauTrueMin, ctauTrueMax);
  //ctauTrueMin -= 0.00001;  ctauTrueMax += 0.00001;
  //cout << "Range from data: ctauTrueMin: " << ctauTrueMin << "  ctauTrueMax: " << ctauTrueMax << endl;
  //ws->var("ctau3Dtrue")->setRange("CtauTrueWindow", ctauTrueMin, ctauTrueMax);
  //Double_t ctau3DTrueMax; Double_t ctau3DTrueMin;
  //ctau3DTrueMax = (double)(ceil(ctauTrueMax));
  //ctau3DTrueMin = (double)(floor(ctauTrueMin));
  Double_t ctau3DTrueMin=1e-6; Double_t ctau3DTrueMax =3;

  //ws->var("ctau3Dtrue")->setRange("ctauTrueRange", ctauTrueMin, ctauTrueMax);
  //ws->var("ctau3Dtrue")->setRange(ctau3DTrueMin, ctau3DTrueMax);
  ws->var("ctau3Dtrue")->setRange(1e-6, ctau3DTrueMax);
  ws->var("ctau3Dtrue")->setRange("ctauTrueRange", 1e-6, ctau3DTrueMax);
  ws->var("ctau3Dtrue")->Print();

  //***********************************************************************
  //**************************** CTAU TRUE FIT ****************************
  //***********************************************************************
  cout << endl << "************ Start MC Ctau True Fit ***************" << endl << endl;
  //MC NP ctau true
  double entries_True = ws->data("reducedDS_MC")->numEntries();
  if(ptLow==6.5&&ptHigh==9){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.642, 1e-6, 1.]");
  ws->factory("lambdaDSS2[0.65, 1e-6, 1.]");
  ws->factory("lambdaDSS3[0.1, 1e-6, 1.]");
  ws->factory("fDSS[0.3, 1e-6, 1.]");
  ws->factory("fDSS1[0.3, 1e-6.,1.]");}
  else if(ptLow==3.5&&ptHigh==5){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*1.01));
  ws->factory("lambdaDSS[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS2[0.13, 1e-6, 1.]");
  ws->factory("fDSS[0.96, 1e-6, 1.]");}
  else if(ptLow==4&&ptHigh==6.5){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.3, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.4, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.5, 1e-6, 1.0]");
  ws->factory("fDSS[0.5, 1e-6, 1.]");
  ws->factory("fDSS1[0.5, 1e-6.,1.]");}
  else if(ptLow==6.5&&ptHigh==12){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.2, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.4, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.5, 1e-6, 1.0]");
  ws->factory("fDSS[0.43, 1e-6, 1.]");
  ws->factory("fDSS1[0.5, 1e-6.,1.]");}
  else if(ptLow==9&&ptHigh==12){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.8, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.42, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.8, 1e-6, 1.0]");
  ws->factory("fDSS[0.7, 1e-6, 1.]");
  ws->factory("fDSS1[0.6, 1e-6.,1.]");}
  else if(ptLow==12&&ptHigh==15){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.3568, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.2, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.8, 1e-6, 1.0]");
  ws->factory("fDSS[0.7, 1e-6, 1.]");
  ws->factory("fDSS1[0.6, 1e-6.,1.]");}
  else if(ptLow==15&&ptHigh==20){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.3568, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.2, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.8, 1e-6, 1.0]");
  ws->factory("fDSS[0.7, 1e-6, 1.]");
  ws->factory("fDSS1[0.6, 1e-6.,1.]");}
  else if(ptLow==20&&ptHigh==50){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.2953, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.459, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.3204, 1e-6, 1.0]");
  ws->factory("fDSS[0.3546, 1e-6, 1.]");
  ws->factory("fDSS1[0.7351, 1e-6.,1.]");}
  else if(ptLow==3&&ptHigh==50){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.13, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.4749, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.21, 1e-6, 1.0]");
  ws->factory("fDSS[0.5, 1e-6, 1.]");
  ws->factory("fDSS1[0.8, 1e-6.,1.]");}
  else if(ptLow==6.5&&ptHigh==50){
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.13, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.4749, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.21, 1e-6, 1.0]");
  ws->factory("fDSS[0.5, 1e-6, 1.]");
  ws->factory("fDSS1[0.8, 1e-6.,1.]");}
  else if(cLow==20&&cHigh==60){
    ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
    ws->factory("lambdaDSS[0.23, 1e-6, 1.0]");
    ws->factory("lambdaDSS2[0.8, 1e-6, 1.0]");
    ws->factory("lambdaDSS3[.21, 1e-6, 1.0]");
    ws->factory("fDSS[0.1, 1e-6, 1.]");
    ws->factory("fDSS1[0.3, 1e-6.,1.]");}
  else if(cLow==60&&cHigh==100){
    ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
    ws->factory("lambdaDSS[0.23, 1e-6, 1.0]");
    ws->factory("lambdaDSS2[0.8, 1e-6, 1.0]");
    ws->factory("lambdaDSS3[.21, 1e-6, 1.0]");
    ws->factory("fDSS[0.1, 1e-6, 1.]");
    ws->factory("fDSS1[0.3, 1e-6.,1.]");}
  else {
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, 0., entries_True*2));
  ws->factory("lambdaDSS[0.34, 1e-6, 1.0]");
  ws->factory("lambdaDSS2[0.362, 1e-6, 1.0]");
  ws->factory("lambdaDSS3[.241, 1e-6, 1.0]");
  ws->factory("fDSS[0.8, 0., 1.]");
  ws->factory("fDSS1[0.8, 0., 1.]");}
  /*
  ws->factory("lambdaDSS[0.5769, 0.5768, 0.577]");
  ws->factory("lambdaDSS2[0.3680, 0.3675, 0.3681]");
  ws->factory("lambdaDSS3[.582, 0.581, 0.583]");
  ws->factory("fDSS[0.8, 0., 2.]");
  ws->factory("fDSS1[0.8, 0., 2.]"); */
  //else {
  //ws->factory("lambdaDSS[0.8, 0.001, 2.0]");
  //ws->factory("lambdaDSS2[0.9, 0.001, 2.0]");
  //ws->factory("lambdaDSS3[1., 0.001, 2.0]");
  //ws->factory("fDSS[0.8, 0., 1.]");
  //ws->factory("fDSS2[0.8, 0., 1.]");
  //ws->factory("f[0.8, 0., 1.]");
  //}
  ws->factory("sigmaMC[0.001, 0.000000001, 1.0]");
  ws->factory("ctauMC[0.0, 0.0, 0.0]");
  // create the PDF
  ws->factory(Form("TruthModel::%s(%s)", "pdfCTAUTRUERES", "ctau3Dtrue"));

  //ws->factory(Form("GaussModel::%s(%s, %s, %s, One, One)", "pdfCTAUTRUERES", "ctau3Dtrue", "ctauMC","sigmaMC"));

  cout << endl << "************ 2 ***************" << endl << endl;
  //ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "ctauTruePdf", "ctau3Dtrue",  
  //      "lambdaDSS", 
  //      "pdfCTAUTRUERES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUEDSS1",
        "ctau3Dtrue",
        "lambdaDSS",
        "pdfCTAUTRUERES"));

  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUEDSS2",
        "ctau3Dtrue",
        "lambdaDSS2",
        "pdfCTAUTRUERES"));


  cout << endl << "************ 3 ***************" << endl << endl;
  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUTRUE1", "pdfCTAUTRUEDSS1", "pdfCTAUTRUEDSS2", "fDSS"));
  cout << endl << "************ 4 ***************" << endl << endl;
//  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUTRUE", "pdfCTAUTRUE2", "pdfCTAUTRUEDSS3", "fDSS1"));
  // 1 expo
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfCTAUTRUETot","pdfCTAUTRUE1", "N_Jpsi_MC"));
  //3 expo
  //ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfCTAUTRUETot","pdfCTAUTRUE", "N_Jpsi_MC"));
  cout << endl << "************ 5 ***************" << endl << endl;

  //ws->var("ctau3Dtrue")->setRange(0.1, 6);
  RooAbsPdf *ctauTrueModel = ctauTrueModel = new RooAddPdf("TrueModel_Tot", "TrueModel_Tot", *ws->pdf("pdfCTAUTRUETot"), *ws->var("N_Jpsi_MC"));
  ws->import(*ctauTrueModel);
  ws->pdf("TrueModel_Tot")->setNormRange("CtauTrueRange");
  //TH1 *hCTrue = (TH1*)ctauTrueModel->createHistogram("ctau3Dtrue",50,50);

  RooDataSet* dataToFit = (RooDataSet*)reducedDS_MC->reduce(Form("ctau3Dtrue>=%.f&&ctau3Dtrue<%.f",ctau3DTrueMin,ctau3DTrueMax))->Clone("reducedDS_MC");
  ws->import(*dataToFit, Rename("dataToFit"));
  TCanvas* c_D =  new TCanvas("canvas_D","My plots",4,420,550,520);
  c_D->cd();
  TPad *pad_D_1 = new TPad("pad_D_1", "pad_D_1", 0, 0.16, 0.98, 1.0);
  pad_D_1->SetTicks(1,1);
  pad_D_1->Draw(); pad_D_1->cd();

  bool isWeighted = ws->data("dataToFit")->isWeighted();
  pad_D_1->cd();
  gPad->SetLogy();
  //RooFitResult* fitCtauTrue = ws->pdf("TrueModel_Tot")->fitTo(*reducedDS_MC, SumW2Error(true), Range(ctau3DTrueMin, ctau3DTrueMax), Extended(kTRUE), NumCPU(nCPU), PrintLevel(3), Save());
  RooFitResult* fitCtauTrue = ws->pdf("TrueModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(nCPU), Range("ctauTrueRange"), PrintLevel(-1));
  fitCtauTrue->Print("v");
  //int nBins = min(int( round((ctau3DTrueMax - ctau3DTrueMin)/0.1) ), 1000);
  int nBins = 70;
  //int nBins = min(int( round((ws->var("ctau3Dtrue")->getMax() - ws->var("ctau3Dtrue")->getMin())/0.1) ), 100);
  RooPlot* myPlot_D = ws->var("ctau3Dtrue")->frame(Bins(nBins), Range(ctau3DTrueMin, ctau3DTrueMax)); // bins
  //RooPlot* myPlot_D = wsmc->var("ctau3Dtrue")->frame(Bins(nBins), Range(ctau3DTrueMin, ctau3DTrueMax)); // bins
  //RooPlot* myPlot_D = wsmc->var("ctau3Dtrue")->frame(Bins(nBins), Range(-1, ctau3DTrueMax)); // bins
  //ws->data("reducedDS_MC")->plotOn(myPlot_D,Name("mcHist_D"));
  myPlot_D->SetTitle("");
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  ws->data("reducedDS_MC")->plotOn(myPlot2_D,Name("MCHist_Tot"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7));
  ws->pdf("TrueModel_Tot")->plotOn(myPlot2_D,Name("MCpdf_Tot"), 
      Normalization(ws->data("reducedDS_MC")->sumEntries(), RooAbsReal::NumEvent), 
      Precision(1e-4), LineColor(kRed+2), NormRange("ctauTrueRange"), Range("ctauTrueRange"));
  //myPlot2_D->GetYaxis()->SetRangeUser(10e-2, ws->data("reducedDS_MC")->sumEntries()*10);
  TH1* h = ws->data("reducedDS_MC")->createHistogram("hist", *ws->var("ctau3Dtrue"), Binning(myPlot2_D->GetNbinsX(),myPlot2_D->GetXaxis()->GetXmin(),myPlot2_D->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));

  Double_t Yup(0.),Ydown(0.);
    Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/0.6)));
    Yup = YMax*TMath::Power((YMax/YMin), (0.3/0.6));

  myPlot2_D->GetYaxis()->SetRangeUser(Ydown,Yup);
  myPlot2_D->GetXaxis()->SetRangeUser(-1, 7);
  myPlot2_D->GetXaxis()->CenterTitle();
  myPlot2_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  myPlot2_D->SetFillStyle(4000);
  myPlot2_D->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_D->GetXaxis()->SetLabelSize(0);
  myPlot2_D->GetXaxis()->SetTitleSize(0);
  myPlot2_D->Draw();

  TLegend* leg_D = new TLegend(text_x+0.25,text_y+0.04,text_x+0.4,text_y-0.1); leg_D->SetTextSize(text_size);
  leg_D->SetTextFont(43);
  leg_D->SetBorderSize(0);
  leg_D->AddEntry(myPlot2_D->findObject("MCHist_Tot"),"Data","pe");
  //leg_D->AddEntry(myPlot2_D->findObject("MCHist_Tot_NoW"),"Data w/o WF","pe");
  leg_D->AddEntry(myPlot2_D->findObject("MCpdf_Tot"),"Total fit","fl");
  //leg_D->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
  leg_D->Draw("same");

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);

  //drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x+0.5,text_y,text_color,text_size);
  //if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  //else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  drawText(Form("#lambdaDSS = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError() ),text_x+0.5,text_y-y_diff,text_color,text_size);
  drawText(Form("#lambdaDSS2 = %.4f #pm %.4f", ws->var("lambdaDSS2")->getVal(), ws->var("lambdaDSS2")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("fDSS = %.4f #pm %.4f", ws->var("fDSS")->getVal(), ws->var("fDSS")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);

  TPad *pad_D_2 = new TPad("pad_D_2", "pad_D_2", 0, 0.006, 0.98, 0.227);
  c_D->cd();
  pad_D_2->Draw();
  pad_D_2->cd();
  pad_D_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_D_2->SetBottomMargin(0.67);
  pad_D_2->SetBottomMargin(0.4);
  pad_D_2->SetFillStyle(4000);
  pad_D_2->SetFrameFillStyle(4000);
  pad_D_2->SetTicks(1,1);

  RooPlot* frameTMP = (RooPlot*)myPlot2_D->Clone("TMP");
  RooHist* hpull_D = frameTMP->pullHist(0,0, true);
  hpull_D->SetMarkerSize(0.8);
  //RooPlot* pullFrame_D = ws->var("ctau3Dtrue")->frame(Title("Pull Distribution"), Range(-1, ctau3DTrueMax)) ;
  //RooPlot* pullFrame_D = ws->var("ctau3Dtrue")->frame(Title("Pull Distribution"), Range(ctau3DTrueMin, ctau3DTrueMax)) ;
  RooPlot* pullFrame_D = ws->var("ctau3Dtrue")->frame(Title("Pull Distribution"), Range(ctau3DTrueMin, ctau3DTrueMax)) ;
  pullFrame_D->addPlotable(hpull_D,"PX");
  pullFrame_D->SetTitle("");
  pullFrame_D->SetTitleSize(0);
  pullFrame_D->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_D->GetYaxis()->SetTitle("Pull") ;
  pullFrame_D->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_D->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_D->GetYaxis()->SetRangeUser(-7,7);
  pullFrame_D->GetYaxis()->CenterTitle();
  pullFrame_D->GetYaxis()->SetTickSize(0.04);
  pullFrame_D->GetYaxis()->SetNdivisions(404);

  pullFrame_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  //pullFrame_D->GetXaxis()->SetRangeUser(-1, 7);
  pullFrame_D->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_D->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_D->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_D->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_D->GetXaxis()->CenterTitle();
  pullFrame_D->GetXaxis()->SetTickSize(0.03);
  pullFrame_D->Draw() ;

  printChi2(ws, pad_D_2, frameTMP, fitCtauTrue, "ctau3DTrue", "MCHist_Tot", "MCpdf_Tot", nBins, false);

  TLine *lD = new TLine(-1, 0, 6., 0);
  lD->SetLineStyle(1);
  lD->Draw("same");
  pad_D_2->Update();

  c_D->Update();
  c_D->SaveAs(Form("figs/2DFit_%s/CtauTrue/ctauTrue_%s_%s.pdf",DATE.Data(),bCont.Data(),kineLabel.Data()));

  //TH1 *h1 = (TH1*)TrueModel_Tot->createHistogram("ctau3Dtrue",50,50);
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_%s_%s.root",DATE.Data(),bCont.Data(),kineLabel.Data()),"RECREATE");
  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauTrue->floatParsFinal());
  ctauTrueModel->Write();
  fitCtauTrue->Write();
  outFile->Close();
  //TFile *hFile = new TFile(Form("ctauTreuHist_%s.root",kineLabel.Data()),"recreate");
  //hCTrue->Write();
  //hFile->Close();
  t -> Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
  cout << endl << "************ Finished MC Ctau True Fit ***************" << endl << endl;

}
