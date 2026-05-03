#include "tdrstyle.C"
void cRAA_pT()
{
    //=========Macro generated from canvas: cRAA/
    //=========  (Tue Jul 16 20:43:30 2024) by ROOT version 6.18/04
    TCanvas *cRAA = new TCanvas("cRAA", "",2883,362,700,700);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cRAA->Range(-7.920792,-0.2282927,41.58416,1.527805);
    cRAA->SetFillColor(0);
    cRAA->SetBorderMode(0);
    cRAA->SetBorderSize(2);
    cRAA->SetTickx(1);
    cRAA->SetTicky(1);
    cRAA->SetLeftMargin(0.16);
    cRAA->SetRightMargin(0.032);
    cRAA->SetTopMargin(0.05);
    cRAA->SetBottomMargin(0.13);
    cRAA->SetFrameFillStyle(0);
    cRAA->SetFrameBorderMode(0);
    cRAA->SetFrameFillStyle(0);
    cRAA->SetFrameBorderMode(0);

    // pt bin
    Double_t x_pt[6] = {7.75,10.5,13.5,17.5,22.5,32.5};
    // Input Raa for Prompt
    Double_t raa_pt_pr[6] = {
        0.2434115,
        0.1785764,
        0.1934771,
        0.2114064,
        0.2413684,
        0.2736468
    };
    Double_t x_pt_err[6] = {0,0,0,0,0,0};
    // Input Raa Error for Prompt
    Double_t raa_pt_pr_err[6] = {
        0.03577576,
        0.01294611,
        0.01636791,
        0.02284889,
        0.04192479,
        0.05337822
    };

    // Input Raa for Nprompt 
    Double_t raa_pt_npr[6] = {
        0.2754807,
        0.217654,
        0.2577349,
        0.2550198,
        0.3157402,
        0.3558492
    };
    // Input Raa Error for Nprompt 
    Double_t raa_pt_npr_err[6] = {
        0.04177038,
        0.01639407,
        0.02114623,
        0.02562665,
        0.04474357,
        0.05656255
    };


    // Input Raa for Prompt
    Double_t raa_pt_pr_sys[6] = {
        0.2434115,
        0.1785764,
        0.1934771,
        0.2114064,
        0.2413684,
        0.2736468
    };
    // Input Raa sys x bin error for Prompt
    Double_t x_pt_err2[6] = {
        1.25,
        1.5,
        1.5,
        2.5,
        2.5,
        7.5
    };
    // Input Raa sys error for Prompt
    Double_t raa_pt_pr_sys_err[6] = {
        0.02599687,
        0.01977726,
        0.03341218,
        0.02190423,
        0.02362113,
        0.02887368
    };

    // Input Raa for NonPrompt
    Double_t raa_pt_npr_sys[6] = {
        0.2754807,
        0.217654,
        0.2577349,
        0.2550198,
        0.3157402,
        0.3558492
    };
    // Input Raa sys error for NonPrompt
    Double_t raa_pt_npr_sys_err[6] = {
        0.04852139,
        0.03220003,
        0.04826296,
        0.02937523,
        0.0272676,
        0.03492241
    };
    TGraphErrors *gRaaPtPr = new TGraphErrors(6,x_pt,raa_pt_pr,x_pt_err,raa_pt_pr_err);
    gRaaPtPr->SetName("gRaaPtPr");
    gRaaPtPr->SetTitle("");
    gRaaPtPr->SetFillStyle(1000);

    gRaaPtPr->SetLineColor(kBlue+2);

    gRaaPtPr->SetMarkerColor(kBlue+2);
    gRaaPtPr->SetMarkerStyle(20);
    gRaaPtPr->SetMarkerSize(1.4);

    TH1F *hRaaPtPr = new TH1F("hRaaPtPr","",100,0,40);
    hRaaPtPr->SetMinimum(0);
    hRaaPtPr->SetMaximum(1.44);
    hRaaPtPr->SetDirectory(0);
    hRaaPtPr->SetStats(0);
    hRaaPtPr->SetLineStyle(0);
    hRaaPtPr->SetMarkerStyle(20);
    hRaaPtPr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRaaPtPr->GetXaxis()->CenterTitle(true);
    hRaaPtPr->GetXaxis()->SetLabelFont(42);
    hRaaPtPr->GetXaxis()->SetLabelOffset(0.007);
    hRaaPtPr->GetXaxis()->SetLabelSize(0.05);
    hRaaPtPr->GetXaxis()->SetTitleSize(0.06);
    hRaaPtPr->GetXaxis()->SetTitleOffset(0.9);
    hRaaPtPr->GetXaxis()->SetTitleFont(42);
    hRaaPtPr->GetYaxis()->SetTitle("R_{AA}");
    hRaaPtPr->GetYaxis()->CenterTitle(true);
    hRaaPtPr->GetYaxis()->SetLabelFont(42);
    hRaaPtPr->GetYaxis()->SetLabelOffset(0.007);
    hRaaPtPr->GetYaxis()->SetLabelSize(0.05);
    hRaaPtPr->GetYaxis()->SetTitleSize(0.06);
    hRaaPtPr->GetYaxis()->SetTitleOffset(1.25);
    hRaaPtPr->GetYaxis()->SetTitleFont(42);
    hRaaPtPr->GetZaxis()->SetLabelFont(42);
    hRaaPtPr->GetZaxis()->SetLabelOffset(0.007);
    hRaaPtPr->GetZaxis()->SetLabelSize(0.05);
    hRaaPtPr->GetZaxis()->SetTitleSize(0.06);
    hRaaPtPr->GetZaxis()->SetTitleOffset(1);
    hRaaPtPr->GetZaxis()->SetTitleFont(42);
    gRaaPtPr->SetHistogram(hRaaPtPr);

    gRaaPtPr->Draw("ap");

    TGraphErrors *gRaaPtNpr = new TGraphErrors(6,x_pt,raa_pt_npr,x_pt_err,raa_pt_npr_err);
    gRaaPtNpr->SetFillStyle(1000);

    gRaaPtNpr->SetLineColor(kRed+2);

    gRaaPtNpr->SetMarkerStyle(21);
    gRaaPtNpr->SetMarkerSize(1.4);

    TH1F *hRaaPtNpr = new TH1F("hRaaPtNpr","Graph",100,5.275,34.975);
    hRaaPtNpr->SetMinimum(0.1801448);
    hRaaPtNpr->SetMaximum(0.4335269);
    hRaaPtNpr->SetDirectory(0);
    hRaaPtNpr->SetStats(0);
    hRaaPtNpr->SetLineStyle(0);
    hRaaPtNpr->SetMarkerStyle(20);
    hRaaPtNpr->GetXaxis()->SetLabelFont(42);
    hRaaPtNpr->GetXaxis()->SetLabelOffset(0.007);
    hRaaPtNpr->GetXaxis()->SetLabelSize(0.05);
    hRaaPtNpr->GetXaxis()->SetTitleSize(0.06);
    hRaaPtNpr->GetXaxis()->SetTitleOffset(0.9);
    hRaaPtNpr->GetXaxis()->SetTitleFont(42);
    hRaaPtNpr->GetYaxis()->SetLabelFont(42);
    hRaaPtNpr->GetYaxis()->SetLabelOffset(0.007);
    hRaaPtNpr->GetYaxis()->SetLabelSize(0.05);
    hRaaPtNpr->GetYaxis()->SetTitleSize(0.06);
    hRaaPtNpr->GetYaxis()->SetTitleOffset(1.25);
    hRaaPtNpr->GetYaxis()->SetTitleFont(42);
    hRaaPtNpr->GetZaxis()->SetLabelFont(42);
    hRaaPtNpr->GetZaxis()->SetLabelOffset(0.007);
    hRaaPtNpr->GetZaxis()->SetLabelSize(0.05);
    hRaaPtNpr->GetZaxis()->SetTitleSize(0.06);
    hRaaPtNpr->GetZaxis()->SetTitleOffset(1);
    hRaaPtNpr->GetZaxis()->SetTitleFont(42);

    gRaaPtNpr->SetHistogram(hRaaPtNpr);

    TGraphErrors *gRaaPtPrSys = new TGraphErrors(6,x_pt,raa_pt_pr_sys,x_pt_err2,raa_pt_pr_sys_err);

    gRaaPtPrSys->SetFillStyle(1000);
    gRaaPtPrSys->SetMarkerStyle(20);

    TH1F *hRaaPtPrSys = new TH1F("hRaaPtPrSys","Graph",100,0,40);
    hRaaPtPrSys->SetMinimum(0);
    hRaaPtPrSys->SetMaximum(1.44);
    hRaaPtPrSys->SetDirectory(0);
    hRaaPtPrSys->SetStats(0);
    hRaaPtPrSys->SetLineStyle(0);
    hRaaPtPrSys->SetMarkerStyle(20);
    hRaaPtPrSys->GetXaxis()->SetLabelFont(42);
    hRaaPtPrSys->GetXaxis()->SetLabelOffset(0.007);
    hRaaPtPrSys->GetXaxis()->SetLabelSize(0.05);
    hRaaPtPrSys->GetXaxis()->SetTitleSize(0.06);
    hRaaPtPrSys->GetXaxis()->SetTitleOffset(0.9);
    hRaaPtPrSys->GetXaxis()->SetTitleFont(42);
    hRaaPtPrSys->GetYaxis()->SetLabelFont(42);
    hRaaPtPrSys->GetYaxis()->SetLabelOffset(0.007);
    hRaaPtPrSys->GetYaxis()->SetLabelSize(0.05);
    hRaaPtPrSys->GetYaxis()->SetTitleSize(0.06);
    hRaaPtPrSys->GetYaxis()->SetTitleOffset(1.25);
    hRaaPtPrSys->GetYaxis()->SetTitleFont(42);
    hRaaPtPrSys->GetZaxis()->SetLabelFont(42);
    hRaaPtPrSys->GetZaxis()->SetLabelOffset(0.007);
    hRaaPtPrSys->GetZaxis()->SetLabelSize(0.05);
    hRaaPtPrSys->GetZaxis()->SetTitleSize(0.06);
    hRaaPtPrSys->GetZaxis()->SetTitleOffset(1);
    hRaaPtPrSys->GetZaxis()->SetTitleFont(42);


    gRaaPtPrSys->SetHistogram(hRaaPtPrSys);
    gRaaPtPrSys->SetLineColor(kBlue-4);
    gRaaPtPrSys->SetFillColorAlpha(kBlue-10,0.40);
    //gRaaPtPrSys->SetFillColorAlpha(kBlue-9,0.40);


    TGraphErrors *gRaaPtNprSys = new TGraphErrors(6,x_pt,raa_pt_npr_sys,x_pt_err2,raa_pt_npr_sys_err);

    gRaaPtNprSys->SetFillStyle(1000);
    gRaaPtNprSys->SetMarkerStyle(20);

    TH1F *hRaaPtNprSys = new TH1F("hRaaPtNprSys","Graph",100,0,40);
    hRaaPtNprSys->SetMinimum(0);
    hRaaPtNprSys->SetMaximum(1.44);
    hRaaPtNprSys->SetDirectory(0);
    hRaaPtNprSys->SetStats(0);
    hRaaPtNprSys->SetLineStyle(0);
    hRaaPtNprSys->SetMarkerStyle(20);
    hRaaPtNprSys->GetXaxis()->SetLabelFont(42);
    hRaaPtNprSys->GetXaxis()->SetLabelOffset(0.007);
    hRaaPtNprSys->GetXaxis()->SetLabelSize(0.05);
    hRaaPtNprSys->GetXaxis()->SetTitleSize(0.06);
    hRaaPtNprSys->GetXaxis()->SetTitleOffset(0.9);
    hRaaPtNprSys->GetXaxis()->SetTitleFont(42);
    hRaaPtNprSys->GetYaxis()->SetLabelFont(42);
    hRaaPtNprSys->GetYaxis()->SetLabelOffset(0.007);
    hRaaPtNprSys->GetYaxis()->SetLabelSize(0.05);
    hRaaPtNprSys->GetYaxis()->SetTitleSize(0.06);
    hRaaPtNprSys->GetYaxis()->SetTitleOffset(1.25);
    hRaaPtNprSys->GetYaxis()->SetTitleFont(42);
    hRaaPtNprSys->GetZaxis()->SetLabelFont(42);
    hRaaPtNprSys->GetZaxis()->SetLabelOffset(0.007);
    hRaaPtNprSys->GetZaxis()->SetLabelSize(0.05);
    hRaaPtNprSys->GetZaxis()->SetTitleSize(0.06);
    hRaaPtNprSys->GetZaxis()->SetTitleOffset(1);
    hRaaPtNprSys->GetZaxis()->SetTitleFont(42);
    hRaaPtNprSys->SetLineColor(kRed-4);
    hRaaPtNprSys->SetFillColorAlpha(kRed-9,0.40);
    gRaaPtNprSys->SetHistogram(hRaaPtNprSys);


    gRaaPtNprSys->SetHistogram(hRaaPtNprSys);
    gRaaPtNprSys->SetLineColor(kRed-4);
    gRaaPtNprSys->SetFillColorAlpha(kRed-10,0.40);
    //gRaaPtNprSys->SetFillColorAlpha(kRed-9,0.40);


    TLegend *leg = new TLegend(0.68,0.72,0.8,0.82,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.029);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(4000);
    TLegendEntry *entry=leg->AddEntry("gRaaPtPr","Prompt #psi(2S)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);

    entry->SetMarkerColor(kBlue+2);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1.4);
    entry->SetTextFont(42);
    entry=leg->AddEntry("gRaaPtNpr","Nonprompt #psi(2S)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);

    entry->SetMarkerColor(kRed+2);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1.4);
    entry->SetTextFont(42);
    leg->Draw();
    TLine *line = new TLine(0,1,40,1);
    line->SetLineStyle(7);
    line->Draw();
    TBox *box = new TBox(0,0.9683139,2,1.031686);

    box->SetFillColor(kGray+2);
    box->Draw();
    TLatex *   tex = new TLatex(0.25,0.85,"Cent. 0-90%");
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(25);
    tex->Draw();
    tex = new TLatex(0.25,0.779,"|y| < 1.6");
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

    gRaaPtPr->SetMarkerColor(kBlue+1);
    gRaaPtNpr->SetMarkerColor(kRed+1);
    gRaaPtNprSys->Draw("5");
    gRaaPtPrSys->Draw("5");
    gRaaPtNpr->Draw("pz same");
    gRaaPtPr->Draw("pz same");

    cRAA->Modified();
    cRAA->cd();
    cRAA->SetSelected(cRAA);

    cRAA->SaveAs("plot_raa_pt.png");
    cRAA->SaveAs("plot_raa_pt.pdf");
    cRAA->Modified();
}
