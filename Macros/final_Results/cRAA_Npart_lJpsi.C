#include "../../tdrstyle.C"
void cRAA_Npart_lJpsi()
{
    //=========Macro generated from canvas: cRAA/
    //=========  (Tue Jul 16 20:49:19 2024) by ROOT version 6.18/04

    TCanvas *cRAA = new TCanvas("cRAA", "",0,53,700,700);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cRAA->Range(-79.20792,-0.2282927,415.8416,1.527805);
    cRAA->SetFillColor(0);
    cRAA->SetBorderMode(0);
    cRAA->SetBorderSize(2);
    cRAA->SetTickx(1);
    cRAA->SetTicky(1);
    cRAA->SetLeftMargin(0.16);
    cRAA->SetRightMargin(0.04);
    cRAA->SetTopMargin(0.05);
    cRAA->SetBottomMargin(0.13);
    cRAA->SetFrameFillStyle(0);
    cRAA->SetFrameBorderMode(0);
    cRAA->SetFrameFillStyle(0);
    cRAA->SetFrameBorderMode(0);

    // Npart bin
    Double_t x_npart[6] = {27.12,87.19,131,188.2,262.3,356.9};
    // Input Raa for Prompt
    Double_t raa_npart_pr[6] = {
        0.4293478,
        0.4964902,
        0.3728716,
        0.4751304,
        0.6624808,
        0.4181731
    };
    Double_t x_npart_err[6] = {0,0,0,0,0,0};
    // Input Raa Error for Prompt
    Double_t raa_npart_pr_err[6] = {
        0.1164218,
        0.08055286,
        0.08985873,
        0.1165608,
        0.2366457,
        0.1108037
    };

    // Input Raa for Nonprompt
    Double_t raa_npart_npr[6] = {
        0.5175139,
        0.6739952,
        0.3850315,
        0.5122114,
        0.6024247,
        0.2802176
    };
    // Input Raa Error for Nonprompt
    Double_t raa_npart_npr_err[6] = {
        0.1455266,
        0.1194875,
        0.09696515,
        0.1311,
        0.2173889,
        0.07986358
    };

    // Input Raa for Prompt
    Double_t raa_npart_pr_sys[6] = {
        0.4293478,
        0.4964902,
        0.3728716,
        0.4751304,
        0.6624808,
        0.4181731
    };
    Double_t x_npart_err2[6] = {5.0,5.0,5.0,5.0,5.0,5.0};
    //Double_t x_npart_err2[6] = {4.3,4.3,4.3,4.3,4.3,4.3};
    // Input Raa systematic error for Prompt
    Double_t raa_npart_pr_sys_err[6] = {
        0.06700721,
        0.106655,
        0.06215021,
        0.08402691,
        0.1229465,
        0.08527856
    };

    // Input Raa for Nonprompt
    Double_t raa_npart_npr_sys[6] = {
        0.5175139,
        0.6739952,
        0.3850315,
        0.5122114,
        0.6024247,
        0.2802176
    };
    // Input Raa systematic error for Prompt
    Double_t raa_npart_npr_sys_err[6] = {
        0.1546556,
        0.2036263,
        0.1172443,
        0.1614459,
        0.1795797,
        0.09753354};
    TGraphErrors *gRaaNpartPr = new TGraphErrors(6,x_npart,raa_npart_pr,x_npart_err,raa_npart_pr_err);
    gRaaNpartPr->SetFillStyle(1000);

    gRaaNpartPr->SetLineColor(kBlue+1);

    gRaaNpartPr->SetMarkerColor(kBlue+1);
    gRaaNpartPr->SetMarkerStyle(20);
    gRaaNpartPr->SetMarkerSize(1.5);

    TH1F *hRaaNpartPr = new TH1F("hRaaNpartPr","",100,0,400);
    hRaaNpartPr->SetMinimum(0);
    hRaaNpartPr->SetMaximum(1.44);
    hRaaNpartPr->SetDirectory(0);
    hRaaNpartPr->SetStats(0);
    hRaaNpartPr->SetLineStyle(0);
    hRaaNpartPr->SetMarkerStyle(20);
    hRaaNpartPr->GetXaxis()->SetTitle("<N_{Part}>");
    hRaaNpartPr->GetXaxis()->CenterTitle(true);
    hRaaNpartPr->GetXaxis()->SetLabelFont(42);
    hRaaNpartPr->GetXaxis()->SetLabelOffset(0.007);
    hRaaNpartPr->GetXaxis()->SetLabelSize(0.05);
    hRaaNpartPr->GetXaxis()->SetTitleSize(0.06);
    hRaaNpartPr->GetXaxis()->SetTitleOffset(0.9);
    hRaaNpartPr->GetXaxis()->SetTitleFont(42);
    hRaaNpartPr->GetYaxis()->SetTitle("R_{AA}");
    hRaaNpartPr->GetYaxis()->CenterTitle(true);
    hRaaNpartPr->GetYaxis()->SetLabelFont(42);
    hRaaNpartPr->GetYaxis()->SetLabelOffset(0.007);
    hRaaNpartPr->GetYaxis()->SetLabelSize(0.05);
    hRaaNpartPr->GetYaxis()->SetTitleSize(0.06);
    hRaaNpartPr->GetYaxis()->SetTitleOffset(1.25);
    hRaaNpartPr->GetYaxis()->SetTitleFont(42);
    hRaaNpartPr->GetZaxis()->SetLabelFont(42);
    hRaaNpartPr->GetZaxis()->SetLabelOffset(0.007);
    hRaaNpartPr->GetZaxis()->SetLabelSize(0.05);
    hRaaNpartPr->GetZaxis()->SetTitleSize(0.06);
    hRaaNpartPr->GetZaxis()->SetTitleOffset(1);
    hRaaNpartPr->GetZaxis()->SetTitleFont(42);
    gRaaNpartPr->SetHistogram(hRaaNpartPr);

    gRaaNpartPr->Draw("ap");

    TGraphErrors *gRaaNpartNpr = new TGraphErrors(6,x_npart,raa_npart_npr,x_npart_err,raa_npart_npr_err);
    gRaaNpartNpr->SetFillStyle(1000);

    gRaaNpartNpr->SetLineColor(kRed+1);

    gRaaNpartNpr->SetMarkerColor(kRed+1);
    gRaaNpartNpr->SetMarkerStyle(21);
    gRaaNpartNpr->SetMarkerSize(1.5);

    TH1F *hRaaNpartNpr = new TH1F("hRaaNpartNpr","Graph",100,0,389.878);
    hRaaNpartNpr->SetMinimum(0.138408);
    hRaaNpartNpr->SetMaximum(0.8817596);
    hRaaNpartNpr->SetDirectory(0);
    hRaaNpartNpr->SetStats(0);
    hRaaNpartNpr->SetLineStyle(0);
    hRaaNpartNpr->SetMarkerStyle(20);
    hRaaNpartNpr->GetXaxis()->SetLabelFont(42);
    hRaaNpartNpr->GetXaxis()->SetLabelOffset(0.007);
    hRaaNpartNpr->GetXaxis()->SetLabelSize(0.05);
    hRaaNpartNpr->GetXaxis()->SetTitleSize(0.06);
    hRaaNpartNpr->GetXaxis()->SetTitleOffset(0.9);
    hRaaNpartNpr->GetXaxis()->SetTitleFont(42);
    hRaaNpartNpr->GetYaxis()->SetLabelFont(42);
    hRaaNpartNpr->GetYaxis()->SetLabelOffset(0.007);
    hRaaNpartNpr->GetYaxis()->SetLabelSize(0.05);
    hRaaNpartNpr->GetYaxis()->SetTitleSize(0.06);
    hRaaNpartNpr->GetYaxis()->SetTitleOffset(1.25);
    hRaaNpartNpr->GetYaxis()->SetTitleFont(42);
    hRaaNpartNpr->GetZaxis()->SetLabelFont(42);
    hRaaNpartNpr->GetZaxis()->SetLabelOffset(0.007);
    hRaaNpartNpr->GetZaxis()->SetLabelSize(0.05);
    hRaaNpartNpr->GetZaxis()->SetTitleSize(0.06);
    hRaaNpartNpr->GetZaxis()->SetTitleOffset(1);
    hRaaNpartNpr->GetZaxis()->SetTitleFont(42);
    gRaaNpartNpr->SetHistogram(hRaaNpartNpr);

    gRaaNpartNpr->Draw("p");

    TGraphErrors *gRaaNpartPrSys = new TGraphErrors(6,x_npart,raa_npart_pr_sys,x_npart_err2,raa_npart_pr_sys_err);

    gRaaNpartPrSys->SetFillColor(kBlue+1);
    gRaaNpartPrSys->SetFillStyle(1000);

    gRaaNpartPrSys->SetLineColor(kBlue+1);
    gRaaNpartPrSys->SetMarkerStyle(20);

    TH1F *hRaaNpartPrSys = new TH1F("hRaaNpartPrSys","Graph",100,0,395.038);
    hRaaNpartPrSys->SetMinimum(0);
    hRaaNpartPrSys->SetMaximum(1.44);
    hRaaNpartPrSys->SetDirectory(0);
    hRaaNpartPrSys->SetStats(0);
    hRaaNpartPrSys->SetLineStyle(0);
    hRaaNpartPrSys->SetMarkerStyle(20);
    hRaaNpartPrSys->GetXaxis()->SetLabelFont(42);
    hRaaNpartPrSys->GetXaxis()->SetLabelOffset(0.007);
    hRaaNpartPrSys->GetXaxis()->SetLabelSize(0.05);
    hRaaNpartPrSys->GetXaxis()->SetTitleSize(0.06);
    hRaaNpartPrSys->GetXaxis()->SetTitleOffset(0.9);
    hRaaNpartPrSys->GetXaxis()->SetTitleFont(42);
    hRaaNpartPrSys->GetYaxis()->SetLabelFont(42);
    hRaaNpartPrSys->GetYaxis()->SetLabelOffset(0.007);
    hRaaNpartPrSys->GetYaxis()->SetLabelSize(0.05);
    hRaaNpartPrSys->GetYaxis()->SetTitleSize(0.06);
    hRaaNpartPrSys->GetYaxis()->SetTitleOffset(1.25);
    hRaaNpartPrSys->GetYaxis()->SetTitleFont(42);
    hRaaNpartPrSys->GetZaxis()->SetLabelFont(42);
    hRaaNpartPrSys->GetZaxis()->SetLabelOffset(0.007);
    hRaaNpartPrSys->GetZaxis()->SetLabelSize(0.05);
    hRaaNpartPrSys->GetZaxis()->SetTitleSize(0.06);
    hRaaNpartPrSys->GetZaxis()->SetTitleOffset(1);
    hRaaNpartPrSys->GetZaxis()->SetTitleFont(42);
    gRaaNpartPrSys->SetHistogram(hRaaNpartPrSys);

    gRaaNpartPrSys->Draw("5");

    TGraphErrors *gRaaNpartNprSys = new TGraphErrors(6,x_npart,raa_npart_npr_sys,x_npart_err2,raa_npart_npr_sys_err);

    gRaaNpartNprSys->SetFillColor(kRed+1);
    gRaaNpartNprSys->SetFillStyle(1000);

    gRaaNpartNprSys->SetLineColor(kRed+1);
    gRaaNpartNprSys->SetMarkerStyle(20);

    TH1F *hRaaNpartNprSys = new TH1F("hRaaNpartNprSys","Graph",100,0,395.038);
    hRaaNpartNprSys->SetMinimum(0.1131903);
    hRaaNpartNprSys->SetMaximum(0.9471153);
    hRaaNpartNprSys->SetDirectory(0);
    hRaaNpartNprSys->SetStats(0);
    hRaaNpartNprSys->SetLineStyle(0);
    hRaaNpartNprSys->SetMarkerStyle(20);
    hRaaNpartNprSys->GetXaxis()->SetLabelFont(42);
    hRaaNpartNprSys->GetXaxis()->SetLabelOffset(0.007);
    hRaaNpartNprSys->GetXaxis()->SetLabelSize(0.05);
    hRaaNpartNprSys->GetXaxis()->SetTitleSize(0.06);
    hRaaNpartNprSys->GetXaxis()->SetTitleOffset(0.9);
    hRaaNpartNprSys->GetXaxis()->SetTitleFont(42);
    hRaaNpartNprSys->GetYaxis()->SetLabelFont(42);
    hRaaNpartNprSys->GetYaxis()->SetLabelOffset(0.007);
    hRaaNpartNprSys->GetYaxis()->SetLabelSize(0.05);
    hRaaNpartNprSys->GetYaxis()->SetTitleSize(0.06);
    hRaaNpartNprSys->GetYaxis()->SetTitleOffset(1.25);
    hRaaNpartNprSys->GetYaxis()->SetTitleFont(42);
    hRaaNpartNprSys->GetZaxis()->SetLabelFont(42);
    hRaaNpartNprSys->GetZaxis()->SetLabelOffset(0.007);
    hRaaNpartNprSys->GetZaxis()->SetLabelSize(0.05);
    hRaaNpartNprSys->GetZaxis()->SetTitleSize(0.06);
    hRaaNpartNprSys->GetZaxis()->SetTitleOffset(1);
    hRaaNpartNprSys->GetZaxis()->SetTitleFont(42);
    gRaaNpartNprSys->SetHistogram(hRaaNpartNprSys);

    gRaaNpartNprSys->Draw("5");

    TLegend *leg = new TLegend(0.68,0.72,0.8,0.82,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.029);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(4000);
    TLegendEntry *entry=leg->AddEntry("gRaaNpartPr","Prompt #psi(2S)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);

    entry->SetMarkerColor(kBlue+1);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1.5);
    entry->SetTextFont(42);
    entry=leg->AddEntry("gRaaNpartNpr","Non Prompt #psi(2S)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);

    entry->SetMarkerColor(kRed+1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1.5);
    entry->SetTextFont(42);
    leg->Draw();
    TLine *line = new TLine(0,1,400,1);
    line->SetLineStyle(7);
    line->Draw();
    TBox *box = new TBox(0,0.9771962,15,1.022804);

    box->SetFillColor(kGray);
    box->Draw();
    TLatex *   tex = new TLatex(0.2,0.85,"3.5 < p_{T} < 40 GeV/c");
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(25);
    tex->Draw();
    tex = new TLatex(0.2,0.799,"1.6 < |y| < 2.4");
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(25);
    tex->Draw();
    tex = new TLatex(0.968,0.9528,"PbPb 1.61 nb^{-1}, pp 300.0 pb^{-1} (5.02 TeV)");
    tex->SetNDC();
    tex->SetTextAlign(31);
    tex->SetTextFont(42);
    tex->SetTextSize(0.03735);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.92164,0.9095026,"CMS");
    tex->SetNDC();
    tex->SetTextAlign(33);
    tex->SetTextFont(61);
    tex->SetTextSize(0.0466875);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.92164,0.8534775,"Preliminary");
    tex->SetNDC();
    tex->SetTextAlign(33);
    tex->SetTextFont(52);
    tex->SetTextSize(0.0354825);
    tex->SetLineWidth(2);
    tex->Draw();
    cRAA->Modified();
    cRAA->cd();
    cRAA->SetSelected(cRAA);

    gRaaNpartPrSys->SetLineColor(kBlue-4);
    gRaaNpartPrSys->SetFillColorAlpha(kBlue-10,0.40);
    gRaaNpartNprSys->SetLineColor(kRed-4);
    gRaaNpartNprSys->SetFillColorAlpha(kRed-10,0.40);


    gRaaNpartPr->SetMarkerColor(kBlue+1);
    gRaaNpartNpr->SetMarkerColor(kRed+1);
    gRaaNpartNprSys->Draw("5");
    gRaaNpartPrSys->Draw("5");
    gRaaNpartNpr->Draw("pz same");
    gRaaNpartPr->Draw("pz same");

    cRAA->Modified();
    cRAA->cd();
    cRAA->SetSelected(cRAA);

    cRAA->SaveAs("plot_raa_npart.png");
    cRAA->SaveAs("plot_raa_npart.pdf");
    cRAA->Modified();

}
