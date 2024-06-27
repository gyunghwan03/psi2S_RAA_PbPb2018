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

void draw_llr_graph_pp_2S()
{
    gStyle->SetOptStat(0);
    setTDRStyle();

    // Forward
    // 4 + 1 bins
    auto canvas_forward = new TCanvas("canvas_forward", "", 700, 700);
    const int n_forward1 = 4;
    const int n_forward2 = 1; // 3.5 ~ 50
    double forward_x1[n_forward1] = {5, 7.75, 10.5, 31};
    double forward_ex1[n_forward1] = {1.5, 1.25, 1.5, 19};
    double forward_y1[n_forward1] = {3, 3, 3, 3};
    double forward_ey1[n_forward1] = {0.1, 0.1, 0.1, 0.1}; // Added for cosmetics.

    double forward_x2[n_forward2] = {26.75};
    double forward_ex2[n_forward2] = {23.25};
    double forward_y2[n_forward2] = {3};
    double forward_ey2[n_forward2] = {0.1};

    TGraph *graph_forward1 = new TGraphErrors(n_forward1, forward_x1, forward_y1, forward_ex1, forward_ey1);
    graph_forward1->SetTitle("LLR Test pp #psi(2S);p_{T} (GeV/c);nParBkg");
    graph_forward1->SetMarkerStyle(27);
    graph_forward1->SetMarkerColor(kRed);
    graph_forward1->SetLineColor(kRed);
    graph_forward1->GetXaxis()->SetLimits(0, 50);
    graph_forward1->SetMinimum(0);
    graph_forward1->SetMaximum(4);

    TGaxis::SetMaxDigits(2); // y축 틱을 정수만 표시하기
    //canvas_forward->SetLeftMargin(0.15);

    graph_forward1->Draw("AP");

    
    TGraph *graph_forward2 = new TGraphErrors(n_forward2, forward_x2, forward_y2, forward_ex2, forward_ey2);

    graph_forward2->SetMarkerStyle(28);
    graph_forward2->SetMarkerColor(kRed);
    graph_forward2->SetLineColor(kRed);
    graph_forward2->Draw("P same");


    // Mid-rapidity
    // 6 + 1 bins
    const int n_mid1 = 6;
    const int n_mid2 = 1; // 3.5 ~ 50
    double mid_x1[n_mid1] = {7.75, 10.5, 13.5, 17.5, 22.5, 37.5};
    double mid_ex1[n_mid1] = {1.25, 1.5, 1.5, 2.5, 2.5, 12.5};
    double mid_y1[n_mid1] = {3, 3, 3, 2, 2, 2};
    double mid_ey1[n_mid1] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
    
    double mid_x2[n_mid2] = {28.25};
    double mid_ex2[n_mid2] = {21.75};
    double mid_y2[n_mid2] = {3};
    double mid_ey2[n_mid2] = {0.03};

    TGraph *graph_mid1 = new TGraphErrors(n_mid1, mid_x1, mid_y1, mid_ex1, mid_ey1);

    graph_mid1->SetMarkerStyle(27);
    graph_mid1->SetMarkerColor(kBlue);
    graph_mid1->SetLineColor(kBlue);
    //graph_mid1->SetLineStyle(2);
    graph_mid1->GetXaxis()->SetLimits(0, 50);
    graph_mid1->SetMinimum(0);
    graph_mid1->SetMaximum(4);

    TGaxis::SetMaxDigits(2); // y축 틱을 정수만 표시하기
    //canvas_mid->SetLeftMargin(0.15);

    graph_mid1->Draw("P same");

    
    TGraph *graph_mid2 = new TGraphErrors(n_mid2, mid_x2, mid_y2, mid_ex2, mid_ey2);

    graph_mid2->SetMarkerStyle(28);
    graph_mid2->SetMarkerColor(kBlue);
    graph_mid2->SetLineColor(kBlue);
    //graph_mid2->SetLineStyle(2);
    graph_mid2->Draw("P same");


    TLegend *legend = new TLegend(0.538, 0.8, 0.94, 0.92);
    SetLegendStyle(legend);
	//drawText(" #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	//drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	//drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    //CMS_lumi_v2mass(c_mid_pt_PR,iPeriod,iPos);	
    legend->AddEntry(graph_forward1, "1.6 < |y| < 2.4", "lp"); // 1 = line
    legend->AddEntry(graph_forward2, "3.5 < p_{T} < 50, 1.6 < |y| < 2.4", "lp");
    legend->AddEntry(graph_mid1, "|y| < 1.6", "lp");
    legend->AddEntry(graph_mid2, "6.5 < p_{T} < 50, |y| < 1.6", "lp");
    //legend->SetTextSize(0.02); 
    //legend->SetFillColor(0);
    //legend->SetFillStyle(0);
    //legend->SetBorderSize(0);
    legend->Draw();

    canvas_forward->SaveAs("llr_result_pp_2S.pdf");
}