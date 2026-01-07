#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <TStyle.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"

using namespace std;
void compare_Eff()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickY(1);

    TFile *f1 = new TFile("./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_251125.root"); // w\o ctau cut
    TFile *f2 = new TFile("./roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_260107_test.root"); // w/ ctau cut

    TH1D *h_eff1 = (TH1D*)f1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");
    TH1D *h_eff2 = (TH1D*)f2->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");

    TH1D *h_eff3 = (TH1D*)f1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6");
    TH1D *h_eff4 = (TH1D*)f2->Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6");

    TCanvas* c_pT_fwd = new TCanvas("c_pT_fwd","",550,520);
    c_pT_fwd->cd();
    TPad *pad_A_1 = new TPad("pad_A_1","",0,0.16,0.98,1);
    pad_A_1->SetBottomMargin(0.056);
    pad_A_1->SetTicks(1,1);
    pad_A_1->Draw();
    pad_A_1->cd();
    h_eff1->SetLineColor(kBlue+2);
    h_eff1->SetMarkerColor(kBlue+2);
    h_eff1->SetMarkerStyle(20);
    h_eff2->SetLineColor(kRed+2);
    h_eff2->SetMarkerColor(kRed+2);
    h_eff2->SetMarkerStyle(21);
    h_eff1->SetTitle("J/#psi Efficiency Comparison with/without ctau cut");
    h_eff1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_eff1->GetYaxis()->SetTitle("Efficiency");
    h_eff1->GetYaxis()->SetRangeUser(0, 1.2);
    h_eff1->GetXaxis()->SetRangeUser(3.5, 40);
    h_eff2->GetXaxis()->SetRangeUser(3.5, 40);
    h_eff1->Draw("PE");
    h_eff2->Draw("PE same"); 

    auto legend = new TLegend(0.55, 0.6, 0.84, 0.84);
    legend->AddEntry(h_eff1, "Without ctau cut", "lep");
    legend->AddEntry(h_eff2, "With ctau cut", "lep");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->Draw("E");
    
    TPad *pad_A_2 = new TPad("pad_A_2","",0,0.06,0.98,0.227);
    c_pT_fwd->cd();
    pad_A_2->Draw();
    pad_A_2->cd();
    pad_A_2->SetTopMargin(0);
    pad_A_2->SetBottomMargin(0.34);
    pad_A_2->SetFillColor(0);
    pad_A_2->SetBorderMode(0);
    pad_A_2->SetBorderSize(2);
    TH1D *h_ratio = (TH1D*)h_eff2->Clone("h_ratio");
    h_ratio->Divide(h_eff1);
    h_ratio->SetTitle("");
    h_ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_ratio->GetYaxis()->SetTitle("Ratio");
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->GetYaxis()->SetRangeUser(0.73, 1.03);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.17);
    h_ratio->GetYaxis()->SetTitleOffset(0.28);
    h_ratio->GetYaxis()->SetLabelSize(0.16);
    h_ratio->GetYaxis()->SetLabelOffset(0.005);
    h_ratio->GetXaxis()->SetTitleSize(0.16);
    h_ratio->GetXaxis()->SetLabelSize(0.16);
    h_ratio->GetXaxis()->SetTitleOffset(1.);
    h_ratio->GetXaxis()->SetTickLength(0.11);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetMarkerStyle(20);
    h_ratio->Draw("PE");   
    jumSun(3.5,0.9,40,0.9);
    c_pT_fwd->SaveAs("./figs/compare_Eff_pT_fwd.png");
    c_pT_fwd->SaveAs("./figs/compare_Eff_pT_fwd.pdf");

    TCanvas* c_pT_mid = new TCanvas("c_pT_mid","",550,520);
    c_pT_mid->cd();
    TPad *pad_B_1 = new TPad("pad_B_1","",0,0.16,0.98,1);
    pad_B_1->SetBottomMargin(0.056);
    pad_B_1->SetTicks(1,1);
    pad_B_1->Draw();
    pad_B_1->cd();
    h_eff3->SetLineColor(kBlue+2);
    h_eff3->SetMarkerColor(kBlue+2);
    h_eff3->SetMarkerStyle(20);
    h_eff4->SetLineColor(kRed+2);
    h_eff4->SetMarkerColor(kRed+2);
    h_eff4->SetMarkerStyle(21);
    h_eff3->SetTitle("J/#psi Efficiency Comparison with/without ctau cut");
    h_eff3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_eff3->GetYaxis()->SetTitle("Efficiency");
    h_eff3->GetYaxis()->SetRangeUser(0, 1.2);
    h_eff3->GetXaxis()->SetRangeUser(6.5, 40);
    h_eff4->GetXaxis()->SetRangeUser(6.5, 40);
    h_eff3->Draw("PE");
    h_eff4->Draw("PE same"); 

    auto legend_B = new TLegend(0.55, 0.6, 0.84, 0.84);
    legend_B->AddEntry(h_eff3, "Without ctau cut", "lep");
    legend_B->AddEntry(h_eff4, "With ctau cut", "lep");
    legend_B->SetBorderSize(0);
    legend_B->SetFillColor(0);
    legend_B->Draw("E");
    
    TPad *pad_B_2 = new TPad("pad_B_2","",0,0.06,0.98,0.227);
    c_pT_mid->cd();
    pad_B_2->Draw();
    pad_B_2->cd();
    pad_B_2->SetTopMargin(0);
    pad_B_2->SetBottomMargin(0.34);
    pad_B_2->SetFillColor(0);
    pad_B_2->SetBorderMode(0);
    pad_B_2->SetBorderSize(2);
    TH1D *h_ratio_mid = (TH1D*)h_eff4->Clone("h_ratio_mid");
    h_ratio_mid->Divide(h_eff3);
    h_ratio_mid->SetTitle("");
    h_ratio_mid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_ratio_mid->GetYaxis()->SetTitle("Ratio");
    h_ratio_mid->GetYaxis()->CenterTitle();
    h_ratio_mid->GetXaxis()->CenterTitle();
    h_ratio_mid->GetYaxis()->SetRangeUser(0.73, 1.03);
    h_ratio_mid->GetYaxis()->SetNdivisions(505);
    h_ratio_mid->GetYaxis()->SetTitleSize(0.17);
    h_ratio_mid->GetYaxis()->SetTitleOffset(0.28);
    h_ratio_mid->GetYaxis()->SetLabelSize(0.16);
    h_ratio_mid->GetYaxis()->SetLabelOffset(0.005);
    h_ratio_mid->GetXaxis()->SetTitleSize(0.16);
    h_ratio_mid->GetXaxis()->SetLabelSize(0.16);
    h_ratio_mid->GetXaxis()->SetTitleOffset(1.);
    h_ratio_mid->GetXaxis()->SetTickLength(0.11);
    h_ratio_mid->SetLineColor(kBlack);
    h_ratio_mid->SetMarkerColor(kBlack);
    h_ratio_mid->SetMarkerStyle(20);
    h_ratio_mid->Draw("PE");   
    jumSun(6.5,0.9,40,0.9);
    c_pT_mid->SaveAs("./figs/compare_Eff_pT_mid.png");
    c_pT_mid->SaveAs("./figs/compare_Eff_pT_mid.pdf");

}