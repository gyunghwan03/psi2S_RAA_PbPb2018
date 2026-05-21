#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"

using namespace std;

valErr getYield(float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0);
double getFrac(float ptLow, float ptHigh, float yLow, float yHigh);
double getFracErr(float ptLow, float ptHigh, float yLow, float yHigh);
void dndpt_y1p6_2p4_Jpsi_pp_noCtau_2exp(int PR=0, int WRITE=1) {

  TString fname;
  if(PR==0) fname = "Prompt";
  else if(PR==1) fname = "NonPrompt";
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);

  TH1::SetDefaultSumw2();

  const int nPtBins=5;
  double ptBin[nPtBins+1]={3.5,5,6.5,9,12,40};
  const int nPtBinsMC=nPtBins;
  double ptBinMC[nPtBinsMC+1];
  for (int i=0; i<=nPtBins; ++i) ptBinMC[i]=ptBin[i];
  const int nYBins=6;
  double yBin[nYBins+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  double frac[nPtBins];
  double fracErr[nPtBins];
  for(int ipt=0;ipt<nPtBins;ipt++) {
    frac[ipt]=getFrac(ptBin[ipt],ptBin[ipt+1],1.6,2.4);
    fracErr[ipt]=getFracErr(ptBin[ipt],ptBin[ipt+1],1.6,2.4);
    cout << ptBin[ipt] << " - " << ptBin[ipt+1]
         << " frac : " << frac[ipt] << " +/- " << fracErr[ipt] << endl;
  }

  float massLow = 2.6;
  float massHigh = 3.5;
  double ptMin = ptBinMC[0];
  double ptMax = ptBinMC[nPtBinsMC];

  TH1D* hptData=new TH1D("hptData",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData1=new TH1D("hptData1",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hfracData=new TH1D("hfracData",";p_{T}(GeV/c);b_fraction",nPtBins,ptBin);
  TH1D* hptMC  = new TH1D("hptMC","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D* hptMC1 = new TH1D("hptMC1","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D* hptMC2 = new TH1D("hptMC2","; p_{T} (GeV/c) ; ", 100, 0, 50 );

  TChain *tree = new TChain("mmepevt");
  TString f1;
  if(PR==0)  f1 ="../skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root";
  else if(PR==1) f1 ="../skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_miniAOD_251103.root";
  tree->Add(f1.Data());

  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu];
  Int_t event;
  Int_t nDimu;
  float vz;
  int recoQQsign[nMaxDimu];

  TBranch *b_event;
  TBranch *b_nDimu;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_recoQQsign;
  TBranch *b_pt;
  TBranch *b_y;

  tree -> SetBranchAddress("event", &event, &b_event);
  tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree -> SetBranchAddress("vz", &vz, &b_vz);
  tree -> SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
  tree -> SetBranchAddress("mass", mass, &b_mass);
  tree -> SetBranchAddress("y", y, &b_y);
  tree -> SetBranchAddress("pt", pt, &b_pt);

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);

    for(int j=0; j<nDimu; j++){
      if (  !( (mass[j] > massLow)
            && (mass[j] < massHigh)
            && ( pt[j] > ptMin)
            && ( pt[j] < ptMax)
            && ( fabs(y[j]) > 1.6)
            && ( fabs(y[j]) < 2.4) )
         )
        continue;
      hptMC->Fill      ( pt[j] );
      hptMC1->Fill     ( pt[j] );
      hptMC2->Fill     ( pt[j] );
    }
  }

  for(int ipt=1;ipt<=nPtBins;ipt++)
  {
    valErr yieldAA = getYield(ptBin[ipt-1],ptBin[ipt],1.6,2.4);

    cout << "yield, pt " << ptBin[ipt-1] << " - " << ptBin[ipt]
         << " : " << yieldAA.val << " +/- " << yieldAA.err
         << ", Frac: " << frac[ipt-1] << " +/- " << fracErr[ipt-1] << endl;

    if(PR==0){
      hptData->SetBinContent(ipt,yieldAA.val*(1-frac[ipt-1]));
      hptData->SetBinError(ipt,yieldAA.err);
      hptData1->SetBinContent(ipt,yieldAA.val*(1-frac[ipt-1]));
      hptData1->SetBinError(ipt,yieldAA.err);
      hfracData->SetBinContent(ipt,frac[ipt-1]);
      hfracData->SetBinError(ipt,fracErr[ipt-1]);
    }
    if(PR==1){
      hptData->SetBinContent(ipt,yieldAA.val*frac[ipt-1]);
      hptData->SetBinError(ipt,yieldAA.err);
      hptData1->SetBinContent(ipt,yieldAA.val*frac[ipt-1]);
      hptData1->SetBinError(ipt,yieldAA.err);
    }
  }

  if (hptMC->Integral() <= 0 || hptData->Integral() <= 0) {
    cout << "ERROR: zero integral histogram, stop." << endl;
    return;
  }

  hptMC->Scale(1./hptMC->Integral());
  hptMC1->Scale(1./hptMC1->Integral());
  hptMC2->Scale(1./hptMC2->Integral());

  hptData->Scale(1./hptData->Integral());
  hptData1->Scale(1./hptData1->Integral());

  TH1ScaleByWidth(hptMC);
  TH1ScaleByWidth(hptData);

  handsomeTH1(hptMC,1);
  handsomeTH1(hptData,1);
  handsomeTH1(hptData1,1);

  TF1* fitRatio1;
  fitRatio1 = new TF1("fitRatio1","[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]",3,40);
  fitRatio1->SetParameters(1.0, 0.30, 0.5, 0.08, 0.0);
  fitRatio1->SetParLimits(1, 0.0, 5.0);
  fitRatio1->SetParLimits(3, 0.0, 5.0);

  TLegend *leg1 = new TLegend(0.65,0.75,0.85,0.85);
  leg1->AddEntry(hptData,"Data","p");
  leg1->AddEntry(hptMC,"MC","l");
  leg1->SetLineColor(kWhite);

  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
  pad_A_1->SetTicks(1,1);
  pad_A_1->Draw(); pad_A_1->cd();
  pad_A_1->SetFillColor(0);
  pad_A_1->SetBorderMode(0);
  pad_A_1->SetBorderSize(2);
  pad_A_1->SetTicks(1,1);
  pad_A_1->SetTopMargin(0.05646528);
  pad_A_1->SetFrameBorderMode(0);
  pad_A_1->SetFrameBorderMode(0);
  hptData->Draw();
  hptMC->Draw("same hist");
  leg1->Draw("same");
  hptData->SetAxisRange(0.,hptMC->GetMaximum()+0.05,"Y");
  hptData->GetXaxis()->SetLabelSize(0);
  hptData->GetYaxis()->SetTitleSize(0.04);
  hptData->GetYaxis()->SetTitleOffset(1.00);
  hptData->GetYaxis()->SetTitle("dN/dp_{T}");

  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
  c_A->cd();
  pad_A_2->Draw();
  pad_A_2->cd();
  pad_A_2->SetFillColor(0);
  pad_A_2->SetBorderMode(0);
  pad_A_2->SetBorderSize(2);
  pad_A_2->SetTicks(1,1);
  pad_A_2->SetBottomMargin(0.4361001);

  hptData1->Divide(hptMC1);
  hptData1->Draw();
  hptData1->GetXaxis()->SetTitleOffset(1.2) ;
  hptData1->GetXaxis()->SetTitleSize(0.15) ;
  hptData1->GetXaxis()->CenterTitle();
  hptData1->GetXaxis()->SetLabelOffset(0.04) ;
  hptData1->GetXaxis()->SetLabelSize(0.15) ;
  hptData1->GetXaxis()->SetTickSize(0.03);
  hptData1->GetYaxis()->SetTickSize(0.04);
  hptData1->GetYaxis()->SetNdivisions(404);
  hptData1->GetYaxis()->SetTitle("Data/MC");
  hptData1->GetYaxis()->SetTitleOffset(0.25) ;
  hptData1->GetYaxis()->SetTitleSize(0.15) ;
  hptData1->GetYaxis()->CenterTitle();
  hptData1->GetYaxis()->SetLabelSize(0.15) ;
  hptData1->GetYaxis()->SetTickSize(0.04);
  hptData1->GetYaxis()->SetNdivisions(404);
  hptData1->SetAxisRange(0,4.,"Y");
  hptData1->Fit(fitRatio1,"IE","",3,12);
  TFitResultPtr r = hptData1->Fit(fitRatio1,"S","",3,40);
  r.Get()->Print("V");
  jumSun(4,1,ptMax,1);

  TCanvas* c2 =  new TCanvas("c2","",604, 0, 800, 600);
  hptMC2->Draw("same hist");

  TCanvas* c_3 =  new TCanvas("c_3","b_fraction",4,4,550,520);
  c_3->cd();
  hfracData->Draw();
  hfracData->GetYaxis()->SetRangeUser(0,1);
  hfracData->SetMarkerStyle(21);
  hfracData->SetLineColor(kRed+2);
  hfracData->SetMarkerColor(kRed+2);
  c_3->SaveAs("./fraction_vs_pt.pdf");

  if(WRITE==1&&PR==0){
    TFile *fJpsipb = new TFile("./ratioDataMC_pp_Jpsi_DATA_noCtau_y1p6_2p4_260514_2exp.root","RECREATE");
    fJpsipb->cd();
    hptData1->SetName("WeightFactor");
    hptData1->Write();
    fitRatio1->SetName("dataMC_Ratio1");
    fitRatio1->Write();
    c_A->Write();
  }
  else if(WRITE==1&&PR==1){
    TFile *fJpsipb = new TFile("./ratioDataMC_pp_BtoJpsi_DATA_noCtau_y1p6_2p4_260514_2exp.root","RECREATE");
    fJpsipb->cd();
    hptData1->SetName("WeightFactor");
    hptData1->Write();
    fitRatio1->SetName("dataMC_Ratio1");
    fitRatio1->Write();
    c_A->Write();
  }

  if(WRITE==1){
    c_A->SaveAs(Form("./dNdpt_plot_%s_Jpsi_pp_y1p6_2p4_noCtau_260514_2exp.pdf",fname.Data()));
    c_A->SaveAs(Form("./dNdpt_plot_%s_Jpsi_pp_y1p6_2p4_noCtau_260514_2exp.png",fname.Data()));
  }
}

valErr getYield(float ptLow, float ptHigh, float yLow, float yHigh) {
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);
  TFile* inf = new TFile(Form("../Macros/pp_Jpsi/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel.Data()));

  valErr ret;
  ret.val = 0;
  ret.err = 0;

  if (!inf || inf->IsZombie()) {
    cout << "cannot open file: " << (inf ? inf->GetName() : "null") << endl;
    return ret;
  }

  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  if (!fitResults) {
    cout << "missing fitResults in: " << inf->GetName() << endl;
    inf->Close();
    return ret;
  }

  ret.val = fitResults->GetBinContent(1);
  ret.err = fitResults->GetBinError(1);
  inf->Close();
  return ret;
}

double getFrac(float ptLow, float ptHigh, float yLow, float yHigh) {
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);
  TFile* inf = new TFile(Form("../Macros/pp_Jpsi/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel.Data()));
  if (!inf || inf->IsZombie()) {
    cout << "cannot open fraction file: " << (inf ? inf->GetName() : "null") << endl;
    return 0;
  }
  TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
  if (!fitResults) {
    cout << "missing 2DfitResults in: " << inf->GetName() << endl;
    inf->Close();
    return 0;
  }
  double frac = fitResults->GetBinContent(1);
  inf->Close();
  return frac;
}

double getFracErr(float ptLow, float ptHigh, float yLow, float yHigh) {
  TString kineLabel = getKineLabelpp(ptLow, ptHigh, yLow, yHigh, 0.0);
  TFile* inf = new TFile(Form("../Macros/pp_Jpsi/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel.Data()));
  if (!inf || inf->IsZombie()) {
    cout << "cannot open fraction file: " << (inf ? inf->GetName() : "null") << endl;
    return 0;
  }
  TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
  if (!fitResults) {
    cout << "missing 2DfitResults in: " << inf->GetName() << endl;
    inf->Close();
    return 0;
  }
  double fracErr = fitResults->GetBinError(1);
  inf->Close();
  return fracErr;
}
