#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TChain.h>
#include <TStyle.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TROOT.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../JpsiUtility.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "../Style.h"

using namespace std;

namespace {
    TGraphErrors* makeGraphFromHist(const TH1D* h, double minPt)
    {
        if (!h) return nullptr;
        auto gr = new TGraphErrors();
        int idx = 0;
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double x = h->GetBinCenter(i);
            double y = h->GetBinContent(i);
            if (x < minPt) continue;
            if (y <= 0) continue;
            double ex = 0.0;
            double ey = h->GetBinError(i);
            if (ey <= 0) {
                ey = std::max(1e-4, 0.05 * y);
            }
            gr->SetPoint(idx, x, y);
            gr->SetPointError(idx, ex, ey);
            idx++;
        }
        return gr;
    }

    void drawCtAuPtPlot(const TString &tag,
                        const TString &stateLabel,
                        const TString &collisionLabel,
                        const TString &rootPath,
                        const TString &yTitle,
                        double fwdMinPt,
                        double midMinPt,
                        const TString &outDir)
    {
        TFile *f = TFile::Open(rootPath, "READ");
        if (!f || f->IsZombie()) {
            cout << "[draw_ctauCut] Cannot open file: " << rootPath << endl;
            return;
        }

        TH1D *h_pt_for = (TH1D*)f->Get("hpt_ctau_1");
        TH1D *h_pt_mid = (TH1D*)f->Get("hpt_ctau_2");
        if (!h_pt_for || !h_pt_mid) {
            cout << "[draw_ctauCut] Missing histograms in: " << rootPath << endl;
            f->Close();
            return;
        }

        TGraphErrors *gr_for = makeGraphFromHist(h_pt_for, fwdMinPt);
        TGraphErrors *gr_mid = makeGraphFromHist(h_pt_mid, midMinPt);
        if (!gr_for || !gr_mid || gr_for->GetN() < 2 || gr_mid->GetN() < 2) {
            cout << "[draw_ctauCut] Not enough points for fit in: " << rootPath << endl;
            f->Close();
            return;
        }

        float pos_x = 0.70;
        float pos_y = 0.85;
        float pos_y_diff = 0.051;
        int text_color = 1;
        float text_size = 12;

        int iPos = 11;
        int iPeriod = 33;        
        if(collisionLabel=="PbPb") { iPeriod = 2;}
        else if(collisionLabel=="pp") { iPeriod = 1;}

        TString fwdName = Form("f1_fwd_%s", tag.Data());
        TString midName = Form("f1_mid_%s", tag.Data());

        TF1* f1_fwd = new TF1(fwdName, "[0]+[1]/x", fwdMinPt, 40);
        f1_fwd->SetParameters(0.01, 0.3);
        f1_fwd->SetParNames("a", "b");
        f1_fwd->SetLineColor(kRed);
        f1_fwd->SetLineWidth(2);

        TF1* f1_mid = new TF1(midName, "[0]+[1]/x", midMinPt, 40);
        f1_mid->SetParameters(0.01, 0.2);
        f1_mid->SetParNames("a", "b");
        f1_mid->SetLineColor(kRed);
        f1_mid->SetLineWidth(2);

        TCanvas* c_fwd = new TCanvas(Form("c_%s_fwd", tag.Data()), "", 550, 520);
        c_fwd->cd();
        gPad->SetLeftMargin(0.15);
        gr_for->SetMarkerColor(kBlue+2);
        gr_for->SetMarkerStyle(20);
        gr_for->SetLineColor(kBlue+2);
        gr_for->SetTitle("");
        gr_for->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        gr_for->GetXaxis()->SetLimits(0,40);
        gr_for->GetXaxis()->CenterTitle();
        gr_for->GetYaxis()->CenterTitle();
        gr_for->GetYaxis()->SetTitle(yTitle);
        gr_for->GetYaxis()->SetRangeUser(0, 0.06);
        gr_for->Draw("AP");
        gr_for->Fit(fwdName, "RM");
        f1_fwd->Draw("same");

        cout << Form("=== %s Forward rapidity fit result ===", stateLabel.Data()) << endl;
        cout << "a = " << f1_fwd->GetParameter(0) << " +/- " << f1_fwd->GetParError(0) << endl;
        cout << "b = " << f1_fwd->GetParameter(1) << " +/- " << f1_fwd->GetParError(1) << endl;
        cout << "Chi2/NDF = " << f1_fwd->GetChisquare() << "/" << f1_fwd->GetNDF() << endl;

        drawText(Form("#chi^{2}/NDF = %.2f/%d", f1_fwd->GetChisquare(), f1_fwd->GetNDF()), pos_x, pos_y, text_color, text_size);
        drawText(Form("a = %.4f #pm %.4f", f1_fwd->GetParameter(0), f1_fwd->GetParError(0)), pos_x, pos_y-pos_y_diff, text_color, text_size);
        drawText(Form("b = %.4f #pm %.4f", f1_fwd->GetParameter(1), f1_fwd->GetParError(1)), pos_x, pos_y-pos_y_diff*2, text_color, text_size);

        CMS_lumi_v2mass(c_fwd, iPeriod, iPos);
        c_fwd->Update();

        c_fwd->SaveAs(Form("%s/ctauCut_ptBin_%s_fwd.pdf", outDir.Data(), tag.Data()));

        TCanvas* c_mid = new TCanvas(Form("c_%s_mid", tag.Data()), "", 550, 520);
        c_mid->cd();
        gPad->SetLeftMargin(0.15);
        gr_mid->SetMarkerColor(kBlue+2);
        gr_mid->SetMarkerStyle(20);
        gr_mid->SetLineColor(kBlue+2);
        gr_mid->SetTitle("");
        gr_mid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        gr_mid->GetXaxis()->SetLimits(0,40);
        gr_mid->GetXaxis()->CenterTitle();
        gr_mid->GetYaxis()->CenterTitle();
        gr_mid->GetYaxis()->SetTitle(yTitle);
        gr_mid->GetYaxis()->SetRangeUser(0, 0.06);
        gr_mid->Draw("AP");
        gr_mid->Fit(midName, "RM");
        f1_mid->Draw("same");

        cout << Form("=== %s Mid rapidity fit result ===", stateLabel.Data()) << endl;
        cout << "a = " << f1_mid->GetParameter(0) << " +/- " << f1_mid->GetParError(0) << endl;
        cout << "b = " << f1_mid->GetParameter(1) << " +/- " << f1_mid->GetParError(1) << endl;
        cout << "Chi2/NDF = " << f1_mid->GetChisquare() << "/" << f1_mid->GetNDF() << endl;

        drawText(Form("#chi^{2}/NDF = %.2f/%d", f1_mid->GetChisquare(), f1_mid->GetNDF()), pos_x, pos_y, text_color, text_size);
        drawText(Form("a = %.4f #pm %.4f", f1_mid->GetParameter(0), f1_mid->GetParError(0)), pos_x, pos_y-pos_y_diff, text_color, text_size);
        drawText(Form("b = %.4f #pm %.4f", f1_mid->GetParameter(1), f1_mid->GetParError(1)), pos_x, pos_y-pos_y_diff*2, text_color, text_size);

        CMS_lumi_v2mass(c_mid, iPeriod, iPos);
        c_mid->Update();

        c_mid->SaveAs(Form("%s/ctauCut_ptBin_%s_mid.pdf", outDir.Data(), tag.Data()));

        f->Close();
    }
}

void draw_ctauCut()
{
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();
    setTDRStyle();
    gStyle->SetOptFit(0);
    writeExtraText = true;

    TString outDir = "./figs_ctauCut";
    gSystem->mkdir(outDir, kTRUE);

    drawCtAuPtPlot("pp_Jpsi", "pp J/#psi", "pp",
                   "./roots_1S_pp/ctau3D_cut_ptBin_PRMC_y0.0-2.4.root",
                   "l_{J/#psi} cut (mm)", 3.5, 6.5, outDir);

    drawCtAuPtPlot("pp_psi2S", "pp #psi(2S)", "pp",
                   "./roots_2S_pp/ctau3D_cut_ptBin_PRMC_y0.0-2.4.root",
                   "l_{#psi(2S)} cut (mm)", 3.5, 6.5, outDir);

    drawCtAuPtPlot("PbPb_Jpsi", "PbPb J/#psi", "PbPb",
                   "./roots_1S_PbPb/ctau3D_cut_ptBin_PRMC_y0.0-2.4.root",
                   "l_{J/#psi} cut (mm)", 3.5, 6.5, outDir);

    drawCtAuPtPlot("PbPb_psi2S", "PbPb #psi(2S)", "PbPb",
                   "./roots_2S_PbPb/ctau3D_cut_ptBin_PRMC_y0.0-2.4.root",
                   "l_{#psi(2S)} cut (mm)", 3.5, 6.5, outDir);
}