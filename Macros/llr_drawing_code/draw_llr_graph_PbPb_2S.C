#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include "TFile.h"
#include "TH1.h"
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"

void draw_llr_graph_PbPb_2S()
{
    gStyle->SetOptStat(0);
    setTDRStyle();

    // Forward pT
    // 4 bins
    auto canvas_pt = new TCanvas("canvas_pt", "", 700, 700);
    const int n_forward1 = 4;
    double forward_x1[n_forward1] = {5, 7.75, 10.5, 31};
    double forward_ex1[n_forward1] = {1.5, 1.25, 1.5, 19}; // error x points
    double forward_y1[n_forward1] = {2, 1, 1, 1};
    double forward_ey1[n_forward1] = {0.1, 0.1, 0.1, 0.1}; // error y points -> Added for cosmetics.

    TGraph *graph_forward1 = new TGraphErrors(n_forward1, forward_x1, forward_y1, forward_ex1, forward_ey1);
    graph_forward1->SetTitle("LLR Test PbPb #psi(2S);p_{T} (GeV/c);nParBkg");
    graph_forward1->SetMarkerStyle(27);
    graph_forward1->SetMarkerColor(kRed);
    graph_forward1->SetLineColor(kRed);
    graph_forward1->GetXaxis()->SetLimits(0, 50);
    graph_forward1->SetMinimum(0);
    graph_forward1->SetMaximum(4);

    TGaxis::SetMaxDigits(2); // y축 틱을 정수만 표시하기
    //canvas_pt->SetLeftMargin(0.15);

    graph_forward1->Draw("AP");

    // Mid-rapidity pT
    // 6 + 1 bins
    const int n_mid1 = 6;
    double mid_x1[n_mid1] = {7.75, 10.5, 13.5, 17.5, 22.5, 37.5};
    double mid_ex1[n_mid1] = {1.25, 1.5, 1.5, 2.5, 2.5, 12.5};
    double mid_y1[n_mid1] = {1, 1, 1, 1, 1, 1};
    double mid_ey1[n_mid1] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

    TGraph *graph_mid1 = new TGraphErrors(n_mid1, mid_x1, mid_y1, mid_ex1, mid_ey1);
    

    graph_mid1->SetMarkerStyle(28);
    graph_mid1->SetMarkerColor(kBlue);
    graph_mid1->SetLineColor(kBlue);
    graph_mid1->GetXaxis()->SetLimits(0, 50);
    graph_mid1->SetMinimum(0);
    graph_mid1->SetMaximum(4);

    TGaxis::SetMaxDigits(2); // digit limit of yaxis
    //canvas_mid->SetLeftMargin(0.15);

    graph_mid1->Draw("P same");

    TLegend *legend = new TLegend(0.6, 0.85, 0.94, 0.92);
    legend->AddEntry(graph_forward1, "1.6 < |y| < 2.4, cent. 10 - 90 %", "lp");
    legend->AddEntry(graph_mid1, "|y| < 1.6, cent. 10 - 90 %", "lp");
    legend->SetTextSize(0.02); 
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    canvas_pt->SaveAs("llr_result_PbPb_2S_pT.pdf");


    // ===== Centrality ===== //
    // Forward cent
    // 6 bins
    auto canvas_cent = new TCanvas("canvas_cent", "", 700, 700);
    canvas_cent->cd();
    const int n_forward2 = 6; // 3.5 ~ 50
    double forward_x2[n_forward2] = {5, 15, 25, 35, 45, 70};
    double forward_ex2[n_forward2] = {5, 5, 5, 5, 5, 20};
    double forward_y2[n_forward2] = {2, 2, 2, 2, 2, 2};
    double forward_ey2[n_forward2] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    TGraph *graph_forward2 = new TGraphErrors(n_forward2, forward_x2, forward_y2, forward_ex2, forward_ey2);
    TGaxis::SetMaxDigits(3);
    graph_forward2->SetTitle("LLR Test PbPb #psi(2S);centrality (%);nParBkg");
    graph_forward2->SetMarkerStyle(27);
    graph_forward2->SetMarkerColor(kRed);
    graph_forward2->SetLineColor(kRed);
    graph_forward2->GetXaxis()->SetLimits(0, 100);
    graph_forward2->SetMinimum(0);
    graph_forward2->SetMaximum(4);

    graph_forward2->Draw("AP");


    // Mid-rapidity cent
    // 6 bins
    const int n_mid2 = 6; 
    double mid_x2[n_mid2] = {5, 15, 25, 35, 45, 70};
    double mid_ex2[n_mid2] = {5, 5, 5, 5, 5, 20};
    double mid_y2[n_mid2] = {1, 1, 1, 1, 1, 1};
    double mid_ey2[n_mid2] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

    TGraph *graph_mid2 = new TGraphErrors(n_mid2, mid_x2, mid_y2, mid_ex2, mid_ey2);

    graph_mid2->SetMarkerStyle(28);
    graph_mid2->SetMarkerColor(kBlue);
    graph_mid2->SetLineColor(kBlue);

    graph_mid2->Draw("P same");

    TLegend *legend_cent = new TLegend(0.62, 0.85, 0.94, 0.92);
    legend_cent->AddEntry(graph_forward2, "3.5 < p_{T} < 50, 1.6 < |y| < 2.4", "lp");
    legend_cent->AddEntry(graph_mid2, "6.5 < p_{T} < 50, |y| < 1.6", "lp");
    
    legend_cent->SetTextSize(0.02); 
    legend_cent->SetFillColor(0);
    legend_cent->SetFillStyle(0);
    legend_cent->SetBorderSize(0);
    legend_cent->Draw();

    canvas_cent->SaveAs("llr_result_PbPb_2S_cent.pdf");
}