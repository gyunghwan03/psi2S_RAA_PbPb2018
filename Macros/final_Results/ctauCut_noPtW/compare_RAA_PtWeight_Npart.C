#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "../../../Style.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"

namespace {
TGraphErrors *makeGraphFromHist(TH1D *h, const double *x, int n, Color_t color, int markerStyle)
{
    double y[6] = {0., 0., 0., 0., 0., 0.};
    double ey[6] = {0., 0., 0., 0., 0., 0.};

    for (int i = 0; i < n; i++) {
        y[i] = h->GetBinContent(i + 1);
        ey[i] = h->GetBinError(i + 1);
    }

    TGraphErrors *g = new TGraphErrors(n, x, y, 0, ey);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(1.4);
    g->SetLineWidth(2);
    return g;
}

void drawCommonText(const char *stateLabel)
{
    drawText(stateLabel, 0.20, 0.84, 1, 26);
    drawText("|y| < 1.6, 6.5 < p_{T} < 40 GeV/c", 0.20, 0.79, 1, 20);
}
}

void compare_RAA_PtWeight_Npart()
{
    gStyle->SetOptStat(0);
    setTDRStyle();

    const int iPeriod = 101;
    const int iPos = 33;
    const int nCentBins = 6;
    const double nPart[nCentBins] = {27.12, 87.19, 131.0, 188.2, 262.3, 356.9};

    TFile *fNoPtW_Jpsi = TFile::Open("roots/RAA_JPsi_midRap_Npart.root");
    TFile *fNoPtW_psi2S = TFile::Open("roots/RAA_psi2S_midRap_Npart_pt6p5_40.root");
    TFile *fPtW_Jpsi = TFile::Open("../ctauCut_260311/roots/RAA_JPsi_midRap_Npart.root");
    TFile *fPtW_psi2S = TFile::Open("../ctauCut_260311/roots/RAA_psi2S_midRap_Npart_pt6p5_40.root");

    if (!fNoPtW_Jpsi || !fNoPtW_psi2S || !fPtW_Jpsi || !fPtW_psi2S) {
        cout << "[ERROR] Failed to open one or more input ROOT files." << endl;
        return;
    }

    TH1D *hJpsiPR_noPtW = (TH1D *)fNoPtW_Jpsi->Get("hRAA_PR");
    TH1D *hJpsiNP_noPtW = (TH1D *)fNoPtW_Jpsi->Get("hRAA_NP");
    TH1D *hPsi2SPR_noPtW = (TH1D *)fNoPtW_psi2S->Get("hRAA_PR");
    TH1D *hPsi2SNP_noPtW = (TH1D *)fNoPtW_psi2S->Get("hRAA_NP");

    TH1D *hJpsiPR_PtW = (TH1D *)fPtW_Jpsi->Get("hRAA_PR");
    TH1D *hJpsiNP_PtW = (TH1D *)fPtW_Jpsi->Get("hRAA_NP");
    TH1D *hPsi2SPR_PtW = (TH1D *)fPtW_psi2S->Get("hRAA_PR");
    TH1D *hPsi2SNP_PtW = (TH1D *)fPtW_psi2S->Get("hRAA_NP");

    if (!hJpsiPR_noPtW || !hJpsiNP_noPtW || !hPsi2SPR_noPtW || !hPsi2SNP_noPtW ||
        !hJpsiPR_PtW || !hJpsiNP_PtW || !hPsi2SPR_PtW || !hPsi2SNP_PtW) {
        cout << "[ERROR] Failed to get hRAA_PR/hRAA_NP from one or more files." << endl;
        return;
    }

    TGraphErrors *gJpsiPR_noPtW = makeGraphFromHist(hJpsiPR_noPtW, nPart, nCentBins, kBlue + 1, 24);
    TGraphErrors *gJpsiPR_PtW = makeGraphFromHist(hJpsiPR_PtW, nPart, nCentBins, kBlue + 1, 20);
    TGraphErrors *gPsi2SPR_noPtW = makeGraphFromHist(hPsi2SPR_noPtW, nPart, nCentBins, kRed + 1, 25);
    TGraphErrors *gPsi2SPR_PtW = makeGraphFromHist(hPsi2SPR_PtW, nPart, nCentBins, kRed + 1, 21);

    TCanvas *cPrompt = new TCanvas("cPrompt_ptwComp", "", 800, 700);
    cPrompt->cd();

    gJpsiPR_PtW->GetXaxis()->SetTitle("<N_{Part}>");
    gJpsiPR_PtW->GetYaxis()->SetTitle("R_{AA}");
    gJpsiPR_PtW->GetXaxis()->SetLimits(0., 400.);
    gJpsiPR_PtW->SetMinimum(0.);
    gJpsiPR_PtW->SetMaximum(1.6);
    gJpsiPR_PtW->SetTitle("");
    gJpsiPR_PtW->Draw("AP");
    gJpsiPR_noPtW->Draw("P");
    gPsi2SPR_PtW->Draw("P");
    gPsi2SPR_noPtW->Draw("P");

    TLegend *legPR = new TLegend(0.54, 0.67, 0.88, 0.88);
    SetLegendStyle(legPR);
    legPR->SetTextSize(0.030);
    legPR->AddEntry(gJpsiPR_PtW, "Prompt J/#psi (PtW)", "pe");
    legPR->AddEntry(gJpsiPR_noPtW, "Prompt J/#psi (no PtW)", "pe");
    legPR->AddEntry(gPsi2SPR_PtW, "Prompt #psi(2S) (PtW)", "pe");
    legPR->AddEntry(gPsi2SPR_noPtW, "Prompt #psi(2S) (no PtW)", "pe");
    legPR->Draw("same");

    TLine *linePR = new TLine(0., 1., 400., 1.);
    linePR->SetLineStyle(2);
    linePR->SetLineColor(kGray + 2);
    linePR->Draw("same");

    drawCommonText("Prompt");
    CMS_lumi_v2mass(cPrompt, iPeriod, iPos);
    cPrompt->SaveAs("figs/compare_RAA_PtWeight_prompt_Npart.pdf");

    TGraphErrors *gJpsiNP_noPtW = makeGraphFromHist(hJpsiNP_noPtW, nPart, nCentBins, kBlue + 1, 24);
    TGraphErrors *gJpsiNP_PtW = makeGraphFromHist(hJpsiNP_PtW, nPart, nCentBins, kBlue + 1, 20);
    TGraphErrors *gPsi2SNP_noPtW = makeGraphFromHist(hPsi2SNP_noPtW, nPart, nCentBins, kRed + 1, 25);
    TGraphErrors *gPsi2SNP_PtW = makeGraphFromHist(hPsi2SNP_PtW, nPart, nCentBins, kRed + 1, 21);

    TCanvas *cNonPrompt = new TCanvas("cNonPrompt_ptwComp", "", 800, 700);
    cNonPrompt->cd();

    gJpsiNP_PtW->GetXaxis()->SetTitle("<N_{Part}>");
    gJpsiNP_PtW->GetYaxis()->SetTitle("R_{AA}");
    gJpsiNP_PtW->GetXaxis()->SetLimits(0., 400.);
    gJpsiNP_PtW->SetMinimum(0.);
    gJpsiNP_PtW->SetMaximum(1.6);
    gJpsiNP_PtW->SetTitle("");
    gJpsiNP_PtW->Draw("AP");
    gJpsiNP_noPtW->Draw("P");
    gPsi2SNP_PtW->Draw("P");
    gPsi2SNP_noPtW->Draw("P");

    TLegend *legNP = new TLegend(0.54, 0.67, 0.88, 0.88);
    SetLegendStyle(legNP);
    legNP->SetTextSize(0.030);
    legNP->AddEntry(gJpsiNP_PtW, "Nonprompt J/#psi (PtW)", "pe");
    legNP->AddEntry(gJpsiNP_noPtW, "Nonprompt J/#psi (no PtW)", "pe");
    legNP->AddEntry(gPsi2SNP_PtW, "Nonprompt #psi(2S) (PtW)", "pe");
    legNP->AddEntry(gPsi2SNP_noPtW, "Nonprompt #psi(2S) (no PtW)", "pe");
    legNP->Draw("same");

    TLine *lineNP = new TLine(0., 1., 400., 1.);
    lineNP->SetLineStyle(2);
    lineNP->SetLineColor(kGray + 2);
    lineNP->Draw("same");

    drawCommonText("Nonprompt");
    CMS_lumi_v2mass(cNonPrompt, iPeriod, iPos);
    cNonPrompt->SaveAs("figs/compare_RAA_PtWeight_nonprompt_Npart.pdf");

    cout << "[INFO] Saved:" << endl;
    cout << "  - figs/compare_RAA_PtWeight_prompt_Npart.pdf" << endl;
    cout << "  - figs/compare_RAA_PtWeight_nonprompt_Npart.pdf" << endl;
}
