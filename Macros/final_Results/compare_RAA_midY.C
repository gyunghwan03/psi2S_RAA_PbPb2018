#include "TH1.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "../Style.h"

void compare_RAA_midY()
{
    gStyle->SetOptStat(0);
    setTDRStyle();
    int iPeriod = 101;
    int iPos = 33;

	float pos_x = 0.25;
    float pos_y = 0.85;
    float pos_y_diff = 0.071;
    int text_color = 1;
    float text_size = 25;

    const int nPtBins=6;
	const int nCentBins=6;
	double binWidth[nPtBins]; double x[nPtBins]; double binWidth_oldPt[5]; double x_oldPt[5];
	double binWidth_old[nCentBins]; double x_old[nCentBins];
	double binWidth_new[nCentBins]; double x_new[nCentBins];

	double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,50};
	double ptBin_old[6] = {6.5,9,12,15,20,30};
	double NpartBin_new[nCentBins+1] = {27.12,87.19,131.0,188.2,262.3,356,9};
	double NpartBin_old[nCentBins] = {21.9,86.9,131.4,189.2,264.2,358.8};

	double RAA_pT_PR_old[5] = {0.109,0.128,0.172,0.142,0.123};
	double RAA_pT_PRerr_old[5] = {0.057,0.034,0.048,0.061,0.034};

	double RAA_pT_PR_new[nPtBins+1]; double RAA_pT_PRerr_new[nPtBins+1];

	double RAA_cent_PR_old[nCentBins] = {0.303,0.191,0.245,0.214,0.116};
	double RAA_cent_PRerr_old[nCentBins] = {0.092,0.092,0.074,0.072,0.053};

	double RAA_cent_PR_new[nCentBins+1]; double RAA_cent_PRerr_new[nCentBins+1];

	TFile *inputFile_pt = new TFile("./roots/RAA_psi2S_midRap_pT_ver2.root");
	TH1D *hRAA_pT_PR_new = (TH1D*) inputFile_pt -> Get("hRAA_PR");
	TFile *inputFile_cent = new TFile("./roots/RAA_psi2S_midRap_Npart_ver2.root");
	TH1D *hRAA_Npart_PR_new = (TH1D*) inputFile_cent -> Get("hRAA_PR");

	for(int i=0; i<nPtBins; i++) {
		x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		RAA_pT_PR_new[i]=hRAA_pT_PR_new->GetBinContent(i+1);
		RAA_pT_PRerr_new[i]=hRAA_pT_PR_new->GetBinError(i+1);
		cout << "RAA = " << RAA_pT_PR_new[i] << endl;
	}
	for(int i=0; i<5; i++){
		x_oldPt[i] = (ptBin_old[i+1]+ptBin_old[i])/2;
		binWidth_oldPt[i] = (ptBin_old[i+1]-ptBin_old[i])/2;
	}


	for(int i=0; i<nCentBins; i++){
		RAA_cent_PR_new[i] = hRAA_Npart_PR_new->GetBinContent(i+1);
		RAA_cent_PRerr_new[i] = hRAA_Npart_PR_new->GetBinError(i+1);
	}

	TGraphErrors *gRaa_pT_old = new TGraphErrors(5,x_oldPt,RAA_pT_PR_old,binWidth_oldPt,RAA_pT_PRerr_old);
    TGraphErrors *gRaa_pT_new = new TGraphErrors(nPtBins,x,RAA_pT_PR_new,binWidth,RAA_pT_PRerr_new);

	gRaa_pT_old->SetPoint(4,(30+20)/2,0.123);
	gRaa_pT_old->SetPointError(4,(30-20)/2,0.109);


	TCanvas *c1 = new TCanvas("c1","",800,700);
	c1->cd();

	gRaa_pT_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gRaa_pT_old->GetXaxis()->CenterTitle();
    gRaa_pT_old->GetYaxis()->SetTitle("R_{AA}");
    gRaa_pT_old->SetTitle("");
    gRaa_pT_old->GetYaxis()->CenterTitle();
    gRaa_pT_old->GetXaxis()->SetLimits(0.,50);
    gRaa_pT_old->SetMinimum(0);
    gRaa_pT_old->SetMaximum(1.44);

	gRaa_pT_old->SetMarkerColor(kGray+2);
    gRaa_pT_old->SetLineColor(kGray+2);
    gRaa_pT_old->SetMarkerStyle(25);
    gRaa_pT_old->SetMarkerSize(1.4);
    gRaa_pT_new->SetMarkerColor(kBlue+2);
    gRaa_pT_new->SetLineColor(kBlue+2);
    gRaa_pT_new->SetMarkerStyle(8);
    gRaa_pT_new->SetMarkerSize(1.4);

	gRaa_pT_old->Draw("AP");
    gRaa_pT_new->Draw("P");
    TLegend *leg = new TLegend(0.66,0.72,0.78,0.82);
    SetLegendStyle(leg);
    leg->AddEntry(gRaa_pT_old,"Prompt #psi(2S) HIN-16-025", "p");
    leg->AddEntry(gRaa_pT_new,"Prompt #psi(2S) NEW", "p");
    leg->Draw("SAME");
    jumSun(0,1,50,1);

	drawText("Cent. 0-90%", pos_x, pos_y, text_color, text_size);
    drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(c1,iPeriod,iPos);

	c1->SaveAs("./figs/RAA_Comparison_pT_midY.pdf");
	c1->SaveAs("./figs/RAA_Comparison_pT_midY.png");

	TGraphErrors *gRaa_cent_old = new TGraphErrors(nCentBins-1,NpartBin_old,RAA_cent_PR_old,0,RAA_cent_PRerr_old);
    TGraphErrors *gRaa_cent_new = new TGraphErrors(nCentBins,NpartBin_new,RAA_cent_PR_new,0,RAA_cent_PRerr_new);
	TGraphErrors *gRaa_cent_lim = new TGraphErrors();
	TGraphErrors *gRaa_cent_lim2 = new TGraphErrors();

	gRaa_cent_lim->SetPoint(0,358.8,0);
	gRaa_cent_lim->SetPointError(0,0,0.138);
	gRaa_cent_lim2->SetPoint(0,358.8,0.138);
	gRaa_cent_lim2->SetPointError(0,5,0);



	TCanvas *c2 = new TCanvas("c2","",800,700);
	c2->cd();

	TArrow *arrow;
	arrow = new TArrow(358.8,0.,358.8,0.01,0.02,"<");

	gRaa_cent_old->GetXaxis()->SetTitle("N_{part}");
	gRaa_cent_old->GetXaxis()->CenterTitle();
	gRaa_cent_old->GetYaxis()->SetTitle("R_{AA}");
	gRaa_cent_old->SetTitle("");
	gRaa_cent_old->GetYaxis()->CenterTitle();
	gRaa_cent_old->GetXaxis()->SetLimits(0.,400);
	gRaa_cent_old->SetMinimum(0);
	gRaa_cent_old->SetMaximum(1.44);

	gRaa_cent_old->SetMarkerColor(kGray+2);
	gRaa_cent_old->SetLineColor(kGray+2);
	gRaa_cent_old->SetMarkerStyle(25);
	gRaa_cent_old->SetMarkerSize(1.4);
	gRaa_cent_lim->SetMarkerColor(kGray+2);
	gRaa_cent_lim->SetLineColor(kGray+2);
	gRaa_cent_lim->SetMarkerStyle(25);
	gRaa_cent_lim->SetMarkerSize(1.4);
	gRaa_cent_new->SetMarkerColor(kBlue+2);
	gRaa_cent_new->SetLineColor(kBlue+2);
	gRaa_cent_new->SetMarkerStyle(8);
	gRaa_cent_new->SetMarkerSize(1.4);

	gRaa_cent_lim2->SetMarkerSize(0);
	gRaa_cent_lim2->SetLineColor(kGray+2);

	arrow->SetAngle(40);
	arrow->SetLineWidth(1);
	arrow->SetLineColor(kGray+2);

	gRaa_cent_old->Draw("AP");
	arrow->Draw("");
	gRaa_cent_lim2->Draw("L");
	gRaa_cent_lim->Draw("L");
	gRaa_cent_new->Draw("P");
    TLegend *leg1 = new TLegend(0.66,0.72,0.78,0.82);
    SetLegendStyle(leg1);
    leg1->AddEntry(gRaa_cent_old,"Prompt #psi(2S) HIN-16-025", "p");
    leg1->AddEntry(gRaa_cent_new,"Prompt #psi(2S) NEW", "p");
    leg1->Draw("SAME");
    jumSun(0,1,400,1);

	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y, text_color, text_size);
    drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(c2,iPeriod,iPos);

	c2->SaveAs("./figs/RAA_Comparison_Npart_midY.pdf");
	c2->SaveAs("./figs/RAA_Comparison_Npart_midY.png");

}
