#include "TH1.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "../Style.h"

void compare_RAA_forY()
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

    const int nPtBins=4;
	const int nCentBins=5;
	double binWidth[nPtBins]; double x[nPtBins];
	double binWidth_old[nCentBins]; double x_old[nCentBins];
	double binWidth_new[nCentBins]; double x_new[nCentBins];

	double ptBin[nPtBins+1] = {3.5,5,6.5,12,50};
	double NpartBin_new[nCentBins] = {27.12,87.19,131.0,188.2,309.6};
	double NpartBin_old[nCentBins] = {32.7,160.3,311.5};

	double RAA_pT_PR_old[nPtBins+1] = {0,0.135,0.247};
	double RAA_pT_PRerr_old[nPtBins+1] = {0,0.079,0.098};

	double RAA_pT_PR_new[nPtBins+1]; double RAA_pT_PRerr_new[nPtBins+1];

	double RAA_cent_PR_old[nCentBins] = {0,0,0};
	double RAA_cent_PRerr_old[nCentBins] = {0.092,0.092,0.074};

	double RAA_cent_PR_new[nCentBins+1]; double RAA_cent_PRerr_new[nCentBins+1];

	TFile *inputFile_pt = new TFile("./roots/RAA_psi2S_forRap_pT_ver2.root");
	TH1D *hRAA_pT_PR_new = (TH1D*) inputFile_pt -> Get("hRAA_PR");
	TFile *inputFile_cent = new TFile("./roots/RAA_psi2S_forRap_Npart_ver2.root");
	TH1D *hRAA_Npart_PR_new = (TH1D*) inputFile_cent -> Get("hRAA_PR");

	for(int i=0; i<nPtBins; i++) {
		x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		RAA_pT_PR_new[i]=hRAA_pT_PR_new->GetBinContent(i+1);
		RAA_pT_PRerr_new[i]=hRAA_pT_PR_new->GetBinError(i+1);
	}

	for(int i=0; i<nCentBins; i++){
		RAA_cent_PR_new[i] = hRAA_Npart_PR_new->GetBinContent(i+1);
		RAA_cent_PRerr_new[i] = hRAA_Npart_PR_new->GetBinError(i+1);
	}

	TGraphErrors *gRaa_pT_old1 = new TGraphErrors();
	TGraphErrors *gRaa_pT_old2 = new TGraphErrors();
    TGraphErrors *gRaa_pT_new = new TGraphErrors(nPtBins,x,RAA_pT_PR_new,binWidth,RAA_pT_PRerr_new);
	TGraphErrors *gRaa_pT_lim = new TGraphErrors();

	gRaa_pT_lim->SetPoint(0,(6.5+3)/2,0.395);
	gRaa_pT_lim->SetPointError(0,1,0);

	gRaa_pT_old1->SetPoint(0,(6.5+3)/2,0);
	gRaa_pT_old2->SetPoint(0,(12+6.5)/2,0.135);
	gRaa_pT_old2->SetPoint(1,(30+12)/2,0.247);
	gRaa_pT_old1->SetPointError(0,0,0.395);
	gRaa_pT_old2->SetPointError(0,(12-6.5)/2,0.079);
	gRaa_pT_old2->SetPointError(1,(30-12)/2,0.098);

	TCanvas *c1 = new TCanvas("c1","",800,700);
	c1->cd();

	TArrow *arrow;
	arrow = new TArrow((6.5+3)/2,0,(6.5+3)/2,0.01,0.03,"<");

	gRaa_pT_old2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gRaa_pT_old2->GetXaxis()->CenterTitle();
    gRaa_pT_old2->GetYaxis()->SetTitle("R_{AA}");
    gRaa_pT_old2->SetTitle("");
    gRaa_pT_old2->GetYaxis()->CenterTitle();
    gRaa_pT_old2->GetXaxis()->SetLimits(0.,50);
    gRaa_pT_old2->SetMinimum(0);
    gRaa_pT_old2->SetMaximum(1.44);

	gRaa_pT_old1->SetMarkerColor(kGray+2);
    gRaa_pT_old1->SetLineColor(kGray+2);
    gRaa_pT_old1->SetMarkerStyle(25);
    gRaa_pT_old1->SetMarkerSize(1.4);
	gRaa_pT_old2->SetMarkerColor(kGray+2);
    gRaa_pT_old2->SetLineColor(kGray+2);
    gRaa_pT_old2->SetMarkerStyle(25);
    gRaa_pT_old2->SetMarkerSize(1.4);
    gRaa_pT_new->SetMarkerColor(kBlue+2);
    gRaa_pT_new->SetLineColor(kBlue+2);
    gRaa_pT_new->SetMarkerStyle(8);
    gRaa_pT_new->SetMarkerSize(1.4);

	gRaa_pT_lim->SetMarkerSize(0);
	gRaa_pT_lim->SetLineColor(kGray+2);


	gRaa_pT_old2->Draw("AP");
	gRaa_pT_old1->Draw("L");
    gRaa_pT_new->Draw("P");
	gRaa_pT_lim->Draw("L");
	arrow->Draw();
    TLegend *leg = new TLegend(0.66,0.72,0.78,0.82);
    SetLegendStyle(leg);
    leg->AddEntry(gRaa_pT_old1,"Prompt #psi(2S) HIN-16-025", "p");
    leg->AddEntry(gRaa_pT_new,"Prompt #psi(2S) NEW", "p");
    leg->Draw("SAME");
    jumSun(0,1,50,1);

	drawText("Cent. 0-90%", pos_x, pos_y, text_color, text_size);
    drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(c1,iPeriod,iPos);

	TGraphErrors *gRaa_cent_old2 = new TGraphErrors();
    TGraphErrors *gRaa_cent_new = new TGraphErrors(nCentBins,NpartBin_new,RAA_cent_PR_new,0,RAA_cent_PRerr_new);
	TGraphErrors *gRaa_cent_lim1 = new TGraphErrors();
	TGraphErrors *gRaa_cent_lim2 = new TGraphErrors();


	gRaa_cent_old2->SetPoint(0,311.5,0.193);
	gRaa_cent_old2->SetPointError(0,0,0.145);

	gRaa_cent_lim1->SetPoint(0,32.7,0.518);
	gRaa_cent_lim2->SetPoint(0,160.3,0.290);
	gRaa_cent_lim1->SetPointError(0,5,0);
	gRaa_cent_lim2->SetPointError(0,5,0);
	
	c1->SaveAs("./figs/RAA_Comparison_pT_forward.pdf");


	TCanvas *c2 = new TCanvas("c2","",800,700);
	c2->cd();

	TArrow *arrow1;
	arrow1 = new TArrow(160.3,0.290,160.3,0.01,0.03,">");
	TArrow *arrow2;
	arrow2 = new TArrow(32.7,0.518,32.7,0.01,0.03,">");

	gRaa_cent_old2->GetXaxis()->SetTitle("N_{part}");
	gRaa_cent_old2->GetXaxis()->CenterTitle();
	gRaa_cent_old2->GetYaxis()->SetTitle("R_{AA}");
	gRaa_cent_old2->SetTitle("");
	gRaa_cent_old2->GetYaxis()->CenterTitle();
	gRaa_cent_old2->GetXaxis()->SetLimits(0.,400);
	gRaa_cent_old2->SetMinimum(0);
	gRaa_cent_old2->SetMaximum(1.44);

	gRaa_cent_old2->SetMarkerColor(kGray+2);
	gRaa_cent_old2->SetLineColor(kGray+2);
	gRaa_cent_old2->SetMarkerStyle(25);
	gRaa_cent_old2->SetMarkerSize(1.4);
	gRaa_cent_new->SetMarkerColor(kBlue+2);
	gRaa_cent_new->SetLineColor(kBlue+2);
	gRaa_cent_new->SetMarkerStyle(8);
	gRaa_cent_new->SetMarkerSize(1.4);

	gRaa_cent_lim1->SetMarkerSize(0);
	gRaa_cent_lim1->SetLineColor(kGray+2);
	gRaa_cent_lim2->SetMarkerSize(0);
	gRaa_cent_lim2->SetLineColor(kGray+2);

	arrow1->SetAngle(40);
	arrow1->SetLineWidth(1);
	arrow1->SetLineColor(kGray+2);
	arrow2->SetAngle(40);
	arrow2->SetLineWidth(1);
	arrow2->SetLineColor(kGray+2);

	gRaa_cent_old2->Draw("AP");
	arrow1->Draw("");
	arrow2->Draw("");
	gRaa_cent_new->Draw("P");
	gRaa_cent_lim1->Draw("L");
	gRaa_cent_lim2->Draw("L");
    TLegend *leg1 = new TLegend(0.66,0.72,0.78,0.82);
    SetLegendStyle(leg1);
    leg1->AddEntry(gRaa_cent_old2,"Prompt #psi(2S) HIN-16-025", "p");
    leg1->AddEntry(gRaa_cent_new,"Prompt #psi(2S) NEW", "p");
    leg1->Draw("SAME");
    jumSun(0,1,400,1);

	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(c2,iPeriod,iPos);

	c2->SaveAs("./figs/RAA_Comparison_Npart_forward.pdf");

}
