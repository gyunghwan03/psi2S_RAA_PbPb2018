#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <TStyle.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../JpsiUtility.h"

using namespace std;

void draw_ctauCut()
{
    gStyle->SetOptStat(0);
    TFile *f_pt_PR = new TFile("./roots_ctau3D/ctau3D_cut_ptBin_PRMC_cent0-180.root");
    TFile *f_pt_NP = new TFile("./roots_ctau3D/ctau3D_cut_ptBin_NPMC_cent0-180.root");
    
    TH1D *h_pt_PR_for = (TH1D*)f_pt_PR->Get("h_ctau_ptBin_for");
    TH1D *h_pt_NP_for = (TH1D*)f_pt_NP->Get("h_ctau_ptBin_for");
    TH1D *h_pt_PR_mid = (TH1D*)f_pt_PR->Get("h_ctau_ptBin_mid");
    TH1D *h_pt_NP_mid = (TH1D*)f_pt_NP->Get("h_ctau_ptBin_mid");    

    // Forward rapidity 피팅 함수: a + b/pT
    TF1* f1_fwd = new TF1("f1_fwd","[0]+[1]/x", 3.5, 40);
    f1_fwd->SetParameters(0.01, 0.1);  // 초기값 설정
    f1_fwd->SetParNames("a", "b");
    f1_fwd->SetLineColor(kRed);
    f1_fwd->SetLineWidth(2);

    // Mid rapidity 피팅 함수: a + b/pT
    TF1* f1_mid = new TF1("f1_mid","[0]+[1]/x", 6.5, 40);
    f1_mid->SetParameters(0.01, 0.1);  // 초기값 설정
    f1_mid->SetParNames("a", "b");
    f1_mid->SetLineColor(kRed);
    f1_mid->SetLineWidth(2);

    TCanvas* c_PR_fwd = new TCanvas("c_PR_fwd","",550,520);
    c_PR_fwd->cd();
    h_pt_PR_for->SetMarkerColor(kBlue+2);
    h_pt_PR_for->SetMarkerStyle(20);
    h_pt_PR_for->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_pt_PR_for->GetXaxis()->CenterTitle();
    h_pt_PR_for->GetYaxis()->CenterTitle();
    h_pt_PR_for->GetYaxis()->SetTitle("l_{J/#psi} cut (mm)");
    h_pt_PR_for->GetYaxis()->SetRangeUser(0, 0.1);
    h_pt_PR_for->Fit("f1_fwd", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_for->Fit("f1_fwd", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_for->Fit("f1_fwd", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_for->Draw("P");
    f1_fwd->Draw("same");
    
    // 피팅 결과 출력
    cout << "=== Forward rapidity fit result ===" << endl;
    cout << "a = " << f1_fwd->GetParameter(0) << " +/- " << f1_fwd->GetParError(0) << endl;
    cout << "b = " << f1_fwd->GetParameter(1) << " +/- " << f1_fwd->GetParError(1) << endl;
    cout << "Chi2/NDF = " << f1_fwd->GetChisquare() << "/" << f1_fwd->GetNDF() << endl;

    drawText(Form("#chi^{2}/NDF = %.2f/%d", f1_fwd->GetChisquare(), f1_fwd->GetNDF()), 0.15, 0.85, 0.03);
    drawText(Form("a = %.4f #pm %.4f", f1_fwd->GetParameter(0), f1_fwd->GetParError(0)), 0.15, 0.80, 0.03);
    drawText(Form("b = %.4f #pm %.4f", f1_fwd->GetParameter(1), f1_fwd->GetParError(1)), 0.15, 0.75, 0.03);
    
    c_PR_fwd->SaveAs("./figs/ctauCut_ptBin_PR_fwd.pdf");

    TCanvas* c_PR_mid = new TCanvas("c_PR_mid","",550,520);
    c_PR_mid->cd();
    h_pt_PR_mid->SetMarkerColor(kBlue+2);
    h_pt_PR_mid->SetMarkerStyle(20);
    h_pt_PR_mid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_pt_PR_mid->GetXaxis()->CenterTitle();
    h_pt_PR_mid->GetYaxis()->CenterTitle();
    h_pt_PR_mid->GetYaxis()->SetTitle("l_{J/#psi} cut (mm)");
    h_pt_PR_mid->GetYaxis()->SetRangeUser(0, 0.1);
    h_pt_PR_mid->Fit("f1_mid", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_mid->Fit("f1_mid", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_mid->Fit("f1_mid", "RM");  // "R": 함수 범위 내에서 피팅
    h_pt_PR_mid->Draw("P");
    f1_mid->Draw("same");
    
    // 피팅 결과 출력
    cout << "=== Mid rapidity fit result ===" << endl;
    cout << "a = " << f1_mid->GetParameter(0) << " +/- " << f1_mid->GetParError(0) << endl;
    cout << "b = " << f1_mid->GetParameter(1) << " +/- " << f1_mid->GetParError(1) << endl;
    cout << "Chi2/NDF = " << f1_mid->GetChisquare() << "/" << f1_mid->GetNDF() << endl;
    
    c_PR_mid->SaveAs("./figs/ctauCut_ptBin_PR_mid.pdf");
}