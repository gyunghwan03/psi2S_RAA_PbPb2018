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

void CtauRes(
    double ptLow=3, double ptHigh=4.5,
    double yLow=1.6, double yHigh=2.4,
    int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false
    )
{
  TString DATE;
  DATE="No_Weight";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/CtauRes",DATE.Data()),kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauRes",DATE.Data()),kTRUE);

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

  TFile* f1; TFile* fMass; TFile* fCErr;
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);

  fMass = new TFile(Form("../../roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Tot = (RooDataSet*)fCErr->Get("dataw_Tot");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Tot);
  cout << "******** New Combined Dataset ***********" << endl;
  ws->var("ctau3DRes")->Print();

  //***********************************************************************
  //**************************** CTAU RES FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start Res Fit *****************" << endl << endl;

  double ctauErrMin;
  double ctauErrMax;
  ctauErrMin = ws->var("ctau3DErr")->getMin();  ctauErrMax = ws->var("ctau3DErr")->getMax();

  RooDataSet *ctauResCutDS =(RooDataSet*)dataw_Tot->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")),*(ws->var("ctau3DErr"))),Form("ctau3DRes<0&&ctau3DErr>%f&&ctau3DErr<%f",ctauErrMin,ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  cout<<"N_Jpsi: "<<ws->var("N_Jpsi")->getVal()<<"+/-"<<ws->var("N_Jpsi")->getError()<<endl;
  if (ptLow==5&&ptHigh==6.5) {
    ws->var("N_Jpsi")->setMax(1000000);
  }  
  // create the variables for this model
  int nGauss = 3;
  if(ptLow==25&&ptHigh==50) {nGauss = 2;}
  ws->factory("One[1.0]");
  ws->factory("ctauRes_mean[0.0]");
  //else {ws->factory("ctau1_CtauRes[0.]");  ws->factory("s1_CtauRes[.6, 0.01, 1.1]");}
  if(ptLow==3&&ptHigh==6.5){
  ws->factory("ctau1_CtauRes[0.]");  
  ws->factory("ctau2_CtauRes[0.]");  //ws->factory("s2_CtauRes[2., 1e-6, 10.]");
  ws->factory("ctau3_CtauRes[0.]");  //ws->factory("s3_CtauRes[3,  1e-6, 10.]");
  ws->factory("ctau4_CtauRes[0.]");  //ws->factory("s4_CtauRes[5.37, 0., 10.]");
  ws->factory("s1_CtauRes[0.6, 1e-3, 1.0]");
  ws->factory("rS21_CtauRes[1.96, 1e-3, 5.0]");
  ws->factory("rS32_CtauRes[2.57, 1e-3, 5.0]");}
  else if(ptLow==30&&ptHigh==50){
  ws->factory("ctau1_CtauRes[0.]");  
  ws->factory("ctau2_CtauRes[0.]");  //ws->factory("s2_CtauRes[2., 1e-6, 10.]");
  ws->factory("ctau3_CtauRes[0.]");  //ws->factory("s3_CtauRes[3,  1e-6, 10.]");
  ws->factory("ctau4_CtauRes[0.]");  //ws->factory("s4_CtauRes[5.37, 0., 10.]");
  ws->factory("s1_CtauRes[0.7, 1e-3, 1.0]");
  ws->factory("rS21_CtauRes[1.96, 1e-3, 5.0]");
  ws->factory("rS32_CtauRes[2.57, 1e-3, 5.0]");
  } else if(ptLow==20&&ptHigh==25){
    ws->factory("ctau1_CtauRes[0.]");  
    ws->factory("ctau2_CtauRes[0.]");  //ws->factory("s2_CtauRes[2., 1e-6, 10.]");
    ws->factory("ctau3_CtauRes[0.]");  //ws->factory("s3_CtauRes[3,  1e-6, 10.]");
    ws->factory("ctau4_CtauRes[0.]");  //ws->factory("s4_CtauRes[5.37, 0., 10.]");
    ws->factory("s1_CtauRes[0.5, 1e-3, 1.0]");
    ws->factory("rS21_CtauRes[0.16, 1e-3, 5.0]");
    ws->factory("rS32_CtauRes[1.7, 1e-3, 4]");
  } else if(ptLow==25&&ptHigh==50){
    ws->factory("ctau1_CtauRes[0.]");  
    ws->factory("ctau2_CtauRes[0.]");  //ws->factory("s2_CtauRes[2., 1e-6, 10.]");
    ws->factory("ctau3_CtauRes[0.]");  //ws->factory("s3_CtauRes[3,  1e-6, 10.]");
    ws->factory("ctau4_CtauRes[0.]");  //ws->factory("s4_CtauRes[5.37, 0., 10.]");
    ws->factory("s1_CtauRes[0.21, 1e-3, 1.0]");
    ws->factory("rS21_CtauRes[1.27, 1e-3, 4.5]");
    ws->factory("rS32_CtauRes[1.7, 1e-3, 4]");
    ws->var("N_Jpsi")->setVal(1200);
  } else{
  ws->factory("ctau1_CtauRes[0.]");  
  ws->factory("ctau2_CtauRes[0.]");  //ws->factory("s2_CtauRes[2., 1e-6, 10.]");
  ws->factory("ctau3_CtauRes[0.]");  //ws->factory("s3_CtauRes[3,  1e-6, 10.]");
  ws->factory("ctau4_CtauRes[0.]");  //ws->factory("s4_CtauRes[5.37, 0., 10.]");
  ws->factory("s1_CtauRes[0.7, 1e-3, 0.95]");
  ws->factory("rS21_CtauRes[1.96, 1e-3, 5.0]");
  ws->factory("rS32_CtauRes[2.57, 1e-3, 5.0]");}

  ws->factory("rS43_CtauRes[1.5, 1.0, 10.0]");
  ws->factory("RooFormulaVar::s2_CtauRes('@0*@1',{rS21_CtauRes,s1_CtauRes})");
  ws->factory("RooFormulaVar::s3_CtauRes('@0*@1',{rS32_CtauRes,s2_CtauRes})");
  ws->factory("RooFormulaVar::s4_CtauRes('@0*@1',{rS43_CtauRes,s3_CtauRes})");
 
  //  ws->factory("f_CtauRes[0.4, 1e-6, 1.]");ws->factory("f2_CtauRes[0.2, 1e-6, 1.]");ws->factory("f3_CtauRes[0.7, 0., 1.]");}
  if(ptLow==3&&ptHigh==50) {
    ws->factory("f_CtauRes[0.2, 1e-6, 1.]");ws->factory("f2_CtauRes[0.5, 1e-6, 1.]");ws->factory("f3_CtauRes[0.8, 1e-6, 1.]");
    }
  else if(ptLow==25&&ptHigh==50) {
    ws->factory("f_CtauRes[0.7, 1e-6, 1.]");ws->factory("f2_CtauRes[0.6, 1e-6, 1.]");ws->factory("f3_CtauRes[0.8, 1e-6, 1.]");
  } 
  else {
    ws->factory("f_CtauRes[0.28, 0., 1.]");ws->factory("f2_CtauRes[0.3, 0., 1.]");ws->factory("f3_CtauRes[0.5, 0., 1.]");
  }
  
  
  // create the three PDFs
  TString varName="ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel1_ctauRes", varName.Data(),
        "ctau1_CtauRes", //"ctau1_CtauRes",
        "s1_CtauRes"
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel2_ctauRes", varName.Data(),
        "ctau2_CtauRes", //"ctau2_CtauRes",
        "s2_CtauRes"
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel3_ctauRes", varName.Data(),
        "ctau3_CtauRes", //"ctau3_CtauRes",
        "s3_CtauRes"
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel4_ctauRes", varName.Data(),
        "ctau4_CtauRes", //"ctau3_CtauRes",
        "s4_CtauRes"
        ));
  // combine the two PDFs
  if(nGauss==4){
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
  }
  else if(nGauss==3){
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if(nGauss==2){
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel2_ctauRes", "f_CtauRes"));}

  RooAddPdf* GaussModel_Tot = new RooAddPdf("GaussModel_Tot");
  RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  ws->import(*ctauResModel);

  TCanvas* c_C =  new TCanvas("canvas_C","My plots",1108,4,550,520);
  c_C->cd();
  TPad *pad_C_1 = new TPad("pad_C_1", "pad_C_1", 0, 0.16, 0.98, 1.0);
  pad_C_1->SetTicks(1,1);
  pad_C_1->Draw(); pad_C_1->cd();
  RooPlot* myPlot_C = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
  myPlot_C->SetTitle("");

  TH1D* hTot = (TH1D*)ws->data("dataw_Tot")->createHistogram(("hTot"), *ws->var("ctau3DRes"),Binning(nCtauResBins,ctauResLow,ctauResHigh));
  double ctauResMin=-10;
  //double ctauResMax = hTot->GetBinCenter(hTot->FindLastBinAbove(1,1));
  double ctauResMax=0;
  cout<<"NBins: "<<hTot->GetNbinsX()<<endl;
  for(int i=0; i<hTot->GetNbinsX()/2; i++){
    //cout<<"Content: "<<hTot->GetBinContent(i)<<endl;
    if(hTot->GetBinContent(i)<=0&&hTot->GetBinContent(i+1)<=0){
      //cout<<"#####"<<i<<": "<<hTot->GetBinLowEdge(i+2)<<endl;
      ctauResMin = hTot->GetBinLowEdge(i+2);
    }
    //if(hTot->GetBinContent(i)>1)ctauResMax = hTot->GetBinCenter(i)+hTot->GetBinWidth(i);
  }
  
  //if (ptLow==12&&ptHigh==50) ctauResMin=-5.6;
  //if (ptLow==15&&ptHigh==20) ctauResMin=-6.42;
  if (ptLow==25&&ptHigh==50) ctauResMin=-4.75;

  ws->var("ctau3DRes")->setRange("ctauResWindow",ctauResMin,0);
  cout<<"Fit Range: "<<ctauResMin<<" - 0"<<endl;
  RooDataSet* dataToFit = (RooDataSet*)ctauResCutDS->reduce(Form("ctau3DRes>=%.f&&ctau3DRes<=0",ctauResMin))->Clone("ctauResCutDS");
  //RooDataSet* dataToFit = (RooDataSet*)ctauResCutDS->reduce(Form("ctau3DRes>=%.f&&ctau3DRes<=0",-5.0))->Clone("ctauResCutDS");
  ws->import(*dataToFit, Rename("dataToFit"));
  pad_C_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_C = (RooPlot*)myPlot_C->Clone();
  bool isWeighted = ws->data("dataToFit")->isWeighted();
  //isWeighted = false;
  //RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(false), Extended(kTRUE), NumCPU(4), PrintLevel(-1));
  RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1), Minimizer("Minuit"), Strategy(2));
  //RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(4), PrintLevel(-1));
  ws->import(*fitCtauRes);
  //setFixedVarsToContantVars(ws);
  //ws->data("dataToFit")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0),
  //    MarkerSize(.7), LineColor(kBlack), MarkerColor(kBlack));
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0),
      MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_ctauRes"), Precision(1e-6), NormRange("ctauResWindow"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm1"), Precision(1e-6), NormRange("ctauResWindow"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm2"), Precision(1e-6), NormRange("ctauResWindow"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm3"), Precision(1e-6), NormRange("ctauResWindow"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue+2));
  if(nGauss==4){ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm4"), Precision(1e-6), NormRange("ctauResWindow"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta+2));}

  //y-axis range
  TH1* h = ws->data("ctauResCutDS")->createHistogram("hist", *ws->var("ctau3DRes"), Binning(myPlot_C->GetNbinsX(),myPlot_C->GetXaxis()->GetXmin(),myPlot_C->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  Yup = YMax*TMath::Power((YMax/YMin), (0.5/(1.0-0.5-0.2)));
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.2/(1.0-0.5-0.2))));
  myPlot2_C->GetYaxis()->SetRangeUser(Ydown,Yup);
  cout<<Ydown<<" to "<<Yup<<endl;
  cout<<"###"<<endl;

  Double_t outTot = ws->data("ctauResCutDS")->numEntries();
  Double_t outRes = ws->data("dataToFit")->numEntries();
  cout<<"Tot evt: ("<<outTot<<")"<<endl;
  cout<<"Res evt: ("<<outRes<<")"<<endl;
  cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<endl;

  if (outRes>0.0) {
    //TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine   *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown*TMath::Power((Yup/Ydown),0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    myPlot2_C->addObject(minline);
    TLine   *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown*TMath::Power((Yup/Ydown),0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    myPlot2_C->addObject(maxline);}
  myPlot2_C->GetXaxis()->CenterTitle();
  myPlot2_C->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  myPlot2_C->SetFillStyle(4000);
  myPlot2_C->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_C->GetXaxis()->SetLabelSize(0);
  myPlot2_C->GetXaxis()->SetTitleSize(0);
  myPlot2_C->Draw();

  //cout<<ws->var("ctau3DRes")->getMin()<<", "<<ws->var("ctau3DRes")->getMax()<<endl;
  cout<<"Min: "<<ctauResMin<<", Max: "<<ctauResMax<<endl;
  //vector<double> rangeErr; rangeErr.push_back(ws->var("ctau3DRes")->getMin()); rangeErr.push_back(ws->var("ctau3DRes")->getMax());

  TLegend* leg_C = new TLegend(text_x+0.29,text_y+0.03,text_x+0.39,text_y-0.17); leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(myPlot2_C->findObject("dataHist_ctauRes"),"Data","pe");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_ctauRes"),"Total PDF","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm1"),"Gauss 1","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm2"),"Gauss 2","l");
  if(nGauss==3)leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm3"),"Gauss 3","l");
  if(nGauss==4)leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm4"),"Gauss 4","l");
  leg_C->Draw("same");
  //cout<<"s2/s1: "<<ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal()<<endl;
  //cout<<"s3/s2: "<<ws->var("s3_CtauRes")->getVal()/ws->var("s2_CtauRes")->getVal()<<endl;
  cout<<"s1: "<<ws->var("s1_CtauRes")->getVal()<<endl;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot-outRes)*100/outTot, "%", outTot-outRes),text_x,text_y-y_diff*3,text_color,text_size);
  //cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<endl;

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError() ),text_x+0.5,text_y,text_color,text_size);
  drawText(Form("s1_{Res} = %.4f #pm %.4f", ws->var("s1_CtauRes")->getVal(), ws->var("s1_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
  drawText(Form("(s2/s1)_{Res} = %.4f #pm %.4f", ws->var("rS21_CtauRes")->getVal(), ws->var("rS21_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  if(nGauss==3){
	  drawText(Form("(s3/s2)_{Res} = %.4f #pm %.4f", ws->var("rS32_CtauRes")->getVal(), ws->var("rS32_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
	  drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
	  drawText(Form("f2_{Res} = %.4f #pm %.4f", ws->var("f2_CtauRes")->getVal(), ws->var("f2_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);}
  else if(nGauss==2) {drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);}
  TPad *pad_C_2 = new TPad("pad_C_2", "pad_C_2", 0, 0.006, 0.98, 0.227);
  c_C->cd();
  pad_C_2->Draw();
  pad_C_2->cd();
  pad_C_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_C_2->SetBottomMargin(0.67);
  pad_C_2->SetBottomMargin(0.4);
  pad_C_2->SetFillStyle(4000);
  pad_C_2->SetFrameFillStyle(4000);
  pad_C_2->SetTicks(1,1);

  RooPlot* frameTMP_C = (RooPlot*)myPlot2_C->Clone("TMP");
  RooHist* hpull_C = frameTMP_C->pullHist("dataHist_ctauRes","modelHist_ctauRes", true);
  hpull_C->SetMarkerSize(0.8);
  //RooPlot* pullFrame_C = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)) ;
  RooPlot* pullFrame_C = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(myPlot_C->GetNbinsX()),Range(myPlot_C->GetXaxis()->GetXmin(),myPlot_C->GetXaxis()->GetXmax())) ;
  pullFrame_C->addPlotable(hpull_C,"PX") ;
  pullFrame_C->SetTitle("");
  pullFrame_C->SetTitleSize(0);
  pullFrame_C->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_C->GetYaxis()->SetTitle("Pull") ;
  pullFrame_C->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_C->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_C->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_C->GetYaxis()->CenterTitle();

  pullFrame_C->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  pullFrame_C->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_C->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_C->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_C->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_C->GetXaxis()->CenterTitle();

  pullFrame_C->GetYaxis()->SetTickSize(0.04);
  pullFrame_C->GetYaxis()->SetNdivisions(404);
  pullFrame_C->GetXaxis()->SetTickSize(0.03);
  pullFrame_C->Draw() ;

  TLine *lC = new TLine(ctauResLow,0, ctauResHigh,0);
  lC->SetLineStyle(1);
  lC->Draw("same");

  printChi2(ws, pad_C_2, frameTMP_C, fitCtauRes, "ctau3DRes", "dataHist_ctauRes", "modelHist_ctauRes", nCtauResBins, false);
  pad_C_2->Update();
  cout << endl << "************* Finished Sig Res Fit ****************" << endl << endl;

  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauRes->floatParsFinal());
  RooDataSet *datasetRes = new RooDataSet("datasetRes","dataset with Resolution Fit result", *fitargs);

  datasetRes->Print("V");
  ctauResModel->Print("V");

  c_C->Update();
  c_C->SaveAs(Form("figs/2DFit_%s/CtauRes/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
  ctauResModel->Write();
  fitCtauRes->Write(); 
  fitCtauRes->Print("V"); 
  //GaussModel_Tot->Write();
  //	ctauResCutDS->Write();
  //datasetRes->Write();
  outFile->Close();
}
