#include "TFile.h"
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"

using namespace std;

void draw_syst_ver2()
{
	gStyle->SetOptStat(0);
    setTDRStyle();

	int iPeriod = 101;
    int iPos = 33;

	const int nFiles = 10;

	//TString SysName[nFiles] = {"Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "Acc", "Eff"};
	//TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "Eff"};
	TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "TNP", "Eff", "Acc"}; //, "b Fraction", "Eff"};
	//TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction"}; //, "b Fraction", "Eff"};

	TFile *in_pt[nFiles];
	TFile *in_cent[nFiles];
	TFile *in_Tot = TFile::Open("syst_roots/total_syst.root");

	in_pt[1] = TFile::Open("syst_roots/syst_pt_sigPDF.root");
	in_pt[2] = TFile::Open("syst_roots/syst_pt_sigPAR.root");
	in_pt[3] = TFile::Open("syst_roots/syst_pt_bkgPDF.root");
	in_pt[4] = TFile::Open("syst_roots/syst_pt_HF.root");
	in_pt[5] = TFile::Open("syst_roots/syst_pt_bFrac.root");
	in_pt[6] = TFile::Open("syst_roots/syst_pt_TNP.root"); // TnP PbPb only
	in_pt[7] = TFile::Open("syst_roots/syst_pt_eff.root");
	in_pt[8] = TFile::Open("syst_roots/syst_pt_acc.root");
	in_pt[9] = TFile::Open("syst_roots/syst_pt_TNP_pp.root"); // TnP pp only

	in_cent[1] = TFile::Open("syst_roots/syst_cent_sigPDF.root");
	in_cent[2] = TFile::Open("syst_roots/syst_cent_sigPAR.root");
	in_cent[3] = TFile::Open("syst_roots/syst_cent_bkgPDF.root");
	in_cent[4] = TFile::Open("syst_roots/syst_cent_HF.root");
	in_cent[5] = TFile::Open("syst_roots/syst_cent_bFrac.root");
	in_cent[6] = TFile::Open("syst_roots/syst_cent_TNP.root"); // TnP PbPb only
	in_cent[7] = TFile::Open("syst_roots/syst_cent_eff.root");	
	in_cent[8] = TFile::Open("syst_roots/syst_cent_acc.root");
	in_cent[9] = TFile::Open("syst_roots/syst_cent_TNP_pp.root"); // TnP pp only

	TH1D *h_mid_pt_PR[nFiles];
	TH1D *h_mid_pt_PR_pp[nFiles];
	TH1D *h_mid_pt_PR_pb[nFiles];
	TH1D *h_mid_pt_NP[nFiles];
	TH1D *h_mid_pt_NP_pp[nFiles];
	TH1D *h_mid_pt_NP_pb[nFiles];
	TH1D *h_fwd_pt_PR[nFiles];
	TH1D *h_fwd_pt_PR_pp[nFiles];
	TH1D *h_fwd_pt_PR_pb[nFiles];
	TH1D *h_fwd_pt_NP[nFiles];
	TH1D *h_fwd_pt_NP_pp[nFiles];
	TH1D *h_fwd_pt_NP_pb[nFiles];
	TH1D *h_mid_cent_PR[nFiles];
	TH1D *h_mid_cent_PR_pp[nFiles];
	TH1D *h_mid_cent_PR_pb[nFiles];
	TH1D *h_mid_cent_NP[nFiles];
	TH1D *h_mid_cent_NP_pp[nFiles];
	TH1D *h_mid_cent_NP_pb[nFiles];
	TH1D *h_fwd_cent_PR[nFiles];
	TH1D *h_fwd_cent_PR_pp[nFiles];
	TH1D *h_fwd_cent_PR_pb[nFiles];
	TH1D *h_fwd_cent_NP[nFiles];
	TH1D *h_fwd_cent_NP_pp[nFiles];
	TH1D *h_fwd_cent_NP_pb[nFiles];

	// Storing total systematic uncertainty
	h_mid_pt_PR[0] = (TH1D*) in_Tot->Get("mid_pt_PR");
	h_mid_pt_PR_pp[0] = (TH1D*) in_Tot->Get("mid_pt_PR_pp");
	h_mid_pt_PR_pb[0] = (TH1D*) in_Tot->Get("mid_pt_PR_pb");
	h_mid_pt_NP[0] = (TH1D*) in_Tot->Get("mid_pt_NP");
	h_mid_pt_NP_pp[0] = (TH1D*) in_Tot->Get("mid_pt_NP_pp");
	h_mid_pt_NP_pb[0] = (TH1D*) in_Tot->Get("mid_pt_NP_pb");
	h_fwd_pt_PR[0] = (TH1D*) in_Tot->Get("fwd_pt_PR");
	h_fwd_pt_PR_pp[0] = (TH1D*) in_Tot->Get("fwd_pt_PR_pp");
	h_fwd_pt_PR_pb[0] = (TH1D*) in_Tot->Get("fwd_pt_PR_pb");
	h_fwd_pt_NP[0] = (TH1D*) in_Tot->Get("fwd_pt_NP");
	h_fwd_pt_NP_pp[0] = (TH1D*) in_Tot->Get("fwd_pt_NP_pp");
	h_fwd_pt_NP_pb[0] = (TH1D*) in_Tot->Get("fwd_pt_NP_pb");
	h_mid_cent_PR[0] = (TH1D*) in_Tot->Get("mid_cent_PR");
	h_mid_cent_PR_pp[0] = (TH1D*) in_Tot->Get("mid_cent_PR_pp");
	h_mid_cent_PR_pb[0] = (TH1D*) in_Tot->Get("mid_cent_PR_pb");
	h_mid_cent_NP[0] = (TH1D*) in_Tot->Get("mid_cent_NP");
	h_mid_cent_NP_pp[0] = (TH1D*) in_Tot->Get("mid_cent_NP_pp");
	h_mid_cent_NP_pb[0] = (TH1D*) in_Tot->Get("mid_cent_NP_pb");
	h_fwd_cent_PR[0] = (TH1D*) in_Tot->Get("fwd_cent_PR");
	h_fwd_cent_PR_pp[0] = (TH1D*) in_Tot->Get("fwd_cent_PR_pp");
	h_fwd_cent_PR_pb[0] = (TH1D*) in_Tot->Get("fwd_cent_PR_pb");
	h_fwd_cent_NP[0] = (TH1D*) in_Tot->Get("fwd_cent_NP");
	h_fwd_cent_NP_pp[0] = (TH1D*) in_Tot->Get("fwd_cent_NP_pp");
	h_fwd_cent_NP_pb[0] = (TH1D*) in_Tot->Get("fwd_cent_NP_pb");

	double mid_pt_PR_max = h_mid_pt_PR[0]->GetBinContent(h_mid_pt_PR[0]->GetMaximumBin()) + 0.1;
	double mid_pt_NP_max = h_mid_pt_NP[0]->GetBinContent(h_mid_pt_NP[0]->GetMaximumBin()) + 0.1;
	double fwd_pt_PR_max = h_fwd_pt_PR[0]->GetBinContent(h_fwd_pt_PR[0]->GetMaximumBin()) + 0.1;
	double fwd_pt_NP_max = h_fwd_pt_NP[0]->GetBinContent(h_fwd_pt_NP[0]->GetMaximumBin()) + 0.1;
	double mid_cent_PR_max = h_mid_cent_PR[0]->GetBinContent(h_mid_cent_PR[0]->GetMaximumBin()) + 0.1;
	double mid_cent_PR_pp_max = h_mid_cent_PR_pp[0]->GetBinContent(h_mid_cent_PR_pp[0]->GetMaximumBin()) + 0.1;
	double mid_cent_NP_max = h_mid_cent_NP[0]->GetBinContent(h_mid_cent_NP[0]->GetMaximumBin()) + 0.1;
	double fwd_cent_PR_max = h_fwd_cent_PR[0]->GetBinContent(h_fwd_cent_PR[0]->GetMaximumBin()) + 0.1;
	double fwd_cent_NP_max = h_fwd_cent_NP[0]->GetBinContent(h_fwd_cent_NP[0]->GetMaximumBin()) + 0.1;

	

	for (int i=1; i<nFiles-1; i++) { // If Acceptance is included, nFiles-1
		if (i==6 || i == 8) {
			continue;
		}
		h_mid_pt_PR[i] = (TH1D*) in_pt[i]->Get("mid_PR");
		h_mid_pt_PR_pp[i] = (TH1D*) in_pt[i]->Get("mid_PR_pp");
		h_mid_pt_PR_pb[i] = (TH1D*) in_pt[i]->Get("mid_PR_pb");
		h_mid_pt_NP[i] = (TH1D*) in_pt[i]->Get("mid_NP");
		h_mid_pt_NP_pp[i] = (TH1D*) in_pt[i]->Get("mid_NP_pp");
		h_mid_pt_NP_pb[i] = (TH1D*) in_pt[i]->Get("mid_NP_pb");
		h_fwd_pt_PR[i] = (TH1D*) in_pt[i]->Get("fwd_PR");
		h_fwd_pt_PR_pp[i] = (TH1D*) in_pt[i]->Get("fwd_PR_pp");
		h_fwd_pt_PR_pb[i] = (TH1D*) in_pt[i]->Get("fwd_PR_pb");
		h_fwd_pt_NP[i] = (TH1D*) in_pt[i]->Get("fwd_NP");
		h_fwd_pt_NP_pp[i] = (TH1D*) in_pt[i]->Get("fwd_NP_pp");
		h_fwd_pt_NP_pb[i] = (TH1D*) in_pt[i]->Get("fwd_NP_pb");
		
		h_mid_cent_PR[i] = (TH1D*) in_cent[i]->Get("mid_PR");
		h_mid_cent_PR_pp[i] = (TH1D*) in_cent[i]->Get("mid_PR_pp");
		h_mid_cent_PR_pb[i] = (TH1D*) in_cent[i]->Get("mid_PR_pb");
		h_mid_cent_NP[i] = (TH1D*) in_cent[i]->Get("mid_NP");
		h_mid_cent_NP_pp[i] = (TH1D*) in_cent[i]->Get("mid_NP_pp");
		h_mid_cent_NP_pb[i] = (TH1D*) in_cent[i]->Get("mid_NP_pb");
		h_fwd_cent_PR[i] = (TH1D*) in_cent[i]->Get("fwd_PR");
		h_fwd_cent_PR_pp[i] = (TH1D*) in_cent[i]->Get("fwd_PR_pp");
		h_fwd_cent_PR_pb[i] = (TH1D*) in_cent[i]->Get("fwd_PR_pb");
		h_fwd_cent_NP[i] = (TH1D*) in_cent[i]->Get("fwd_NP");
		h_fwd_cent_NP_pp[i] = (TH1D*) in_cent[i]->Get("fwd_NP_pp");
		h_fwd_cent_NP_pb[i] = (TH1D*) in_cent[i]->Get("fwd_NP_pb");
		

		//cout << h_mid_pt_PR[i]->GetBinContent(1) << endl;
		//cout << h_mid_pt_NP[i]->GetBinContent(1) << endl;
	}

	// Fill 6th and 8th components
	h_mid_pt_PR[6] = (TH1D*) in_pt[6]->Get("mid_PR");
	h_mid_pt_PR_pp[6] = (TH1D*) in_pt[9]->Get("mid_PR");
	h_mid_pt_PR_pb[6] = (TH1D*) in_pt[6]->Get("mid_PR");
	h_mid_pt_NP[6] = (TH1D*) in_pt[6]->Get("mid_NP");
	h_mid_pt_NP_pp[6] = (TH1D*) in_pt[9]->Get("mid_NP");
	h_mid_pt_NP_pb[6] = (TH1D*) in_pt[6]->Get("mid_NP");
	h_fwd_pt_PR[6] = (TH1D*) in_pt[6]->Get("fwd_PR");
	h_fwd_pt_PR_pp[6] = (TH1D*) in_pt[9]->Get("fwd_PR");
	h_fwd_pt_PR_pb[6] = (TH1D*) in_pt[6]->Get("fwd_NP");
	h_fwd_pt_NP[6] = (TH1D*) in_pt[6]->Get("fwd_NP");
	h_fwd_pt_NP_pp[6] = (TH1D*) in_pt[9]->Get("fwd_NP");
	h_fwd_pt_NP_pb[6] = (TH1D*) in_pt[6]->Get("fwd_NP");
	h_mid_cent_PR[6] = (TH1D*) in_cent[6]->Get("mid_PR");
	h_mid_cent_PR_pp[6] = (TH1D*) in_cent[9]->Get("mid_PR");
	h_mid_cent_PR_pb[6] = (TH1D*) in_cent[6]->Get("mid_PR");
	h_mid_cent_NP[6] = (TH1D*) in_cent[6]->Get("mid_NP");
	h_mid_cent_NP_pp[6] = (TH1D*) in_cent[9]->Get("mid_NP");
	h_mid_cent_NP_pb[6] = (TH1D*) in_cent[6]->Get("mid_NP");
	h_fwd_cent_PR[6] = (TH1D*) in_cent[6]->Get("fwd_PR");
	h_fwd_cent_PR_pp[6] = (TH1D*) in_cent[9]->Get("fwd_PR");
	h_fwd_cent_PR_pb[6] = (TH1D*) in_cent[6]->Get("fwd_PR");
	h_fwd_cent_NP[6] = (TH1D*) in_cent[6]->Get("fwd_NP");
	h_fwd_cent_NP_pp[6] = (TH1D*) in_cent[9]->Get("fwd_NP");
	h_fwd_cent_NP_pb[6] = (TH1D*) in_cent[6]->Get("fwd_NP");

	h_mid_pt_PR[8] = (TH1D*) in_pt[8]->Get("mid"); //Acc
	h_mid_pt_PR_pp[8] = (TH1D*) in_pt[8]->Get("mid_pp"); //Acc
	h_mid_pt_PR_pb[8] = (TH1D*) in_pt[8]->Get("mid_pb"); //Acc
	h_mid_pt_NP[8] = (TH1D*) in_pt[8]->Get("mid"); //Acc
	h_mid_pt_NP_pp[8] = (TH1D*) in_pt[8]->Get("mid_ppp"); //Acc
	h_mid_pt_NP_pb[8] = (TH1D*) in_pt[8]->Get("mid_pb"); //Acc
	h_fwd_pt_PR[8] = (TH1D*) in_pt[8]->Get("fwd"); //Acc
	h_fwd_pt_PR_pp[8] = (TH1D*) in_pt[8]->Get("fwd_pp"); //Acc
	h_fwd_pt_PR_pb[8] = (TH1D*) in_pt[8]->Get("fwd_pb"); //Acc
	h_fwd_pt_NP[8] = (TH1D*) in_pt[8]->Get("fwd"); //Acc
	h_fwd_pt_NP_pp[8] = (TH1D*) in_pt[8]->Get("fwd_pp"); //Acc
	h_fwd_pt_NP_pb[8] = (TH1D*) in_pt[8]->Get("fwd_pb"); //Acc
	h_mid_cent_PR[8] = (TH1D*) in_cent[8]->Get("mid");
	h_mid_cent_PR_pp[8] = (TH1D*) in_cent[8]->Get("mid_pp");
	h_mid_cent_PR_pb[8] = (TH1D*) in_cent[8]->Get("mid_pb");
	h_mid_cent_NP[8] = (TH1D*) in_cent[8]->Get("mid");
	h_mid_cent_NP_pp[8] = (TH1D*) in_cent[8]->Get("mid_pp");
	h_mid_cent_NP_pb[8] = (TH1D*) in_cent[8]->Get("mid_pb");
	h_fwd_cent_PR[8] = (TH1D*) in_cent[8]->Get("fwd");
	h_fwd_cent_PR_pp[8] = (TH1D*) in_cent[8]->Get("fwd_pp");
	h_fwd_cent_PR_pb[8] = (TH1D*) in_cent[8]->Get("fwd_pb");
	h_fwd_cent_NP[8] = (TH1D*) in_cent[8]->Get("fwd");
	h_fwd_cent_NP_pp[8] = (TH1D*) in_cent[8]->Get("fwd_pp");
	h_fwd_cent_NP_pb[8] = (TH1D*) in_cent[8]->Get("fwd_pb");

	
	const int NcentBins = 6;
	double centBin[NcentBins+1] = {0,10,20,30,40,50,90};
	h_mid_cent_PR[6] = new TH1D("mid_PR","mid_PR",NcentBins,centBin);
	h_mid_cent_NP[6] = new TH1D("mid_NP","mid_NP",NcentBins,centBin);
	h_fwd_cent_PR[6] = new TH1D("fwd_PR","fwd_PR",NcentBins,centBin);
	h_fwd_cent_NP[6] = new TH1D("fwd_NP","fwd_NP",NcentBins,centBin);

	h_mid_pt_PR[6] = new TH1D("mid_PR","mid_PR",NcentBins,centBin);
	h_mid_pt_NP[6] = new TH1D("mid_NP","mid_NP",NcentBins,centBin);
	h_fwd_pt_PR[6] = new TH1D("fwd_PR","fwd_PR",NcentBins,centBin);
	h_fwd_pt_NP[6] = new TH1D("fwd_NP","fwd_NP",NcentBins,centBin);

	for(int i=1; i<NcentBins+1; i++) {
		// pt
		double TnP_pb = h_mid_pt_PR_pb[6]->GetBinContent(i);
		double TnP_pp = h_mid_pt_PR_pp[6]->GetBinContent(1); // pp input has one TnP bin
		double TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_mid_pt_PR[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_mid_pt_NP_pb[6]->GetBinContent(i);
		TnP_pp = h_mid_pt_NP_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_mid_pt_NP[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_fwd_pt_PR_pb[6]->GetBinContent(i);
		TnP_pp = h_fwd_pt_PR_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_fwd_pt_PR[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_fwd_pt_NP_pb[6]->GetBinContent(i);
		TnP_pp = h_fwd_pt_NP_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_fwd_pt_NP[6]->SetBinContent(i, TnP_tot);

		// cent
		TnP_pb = h_mid_cent_PR_pb[6]->GetBinContent(i);
		TnP_pp = h_mid_cent_PR_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_mid_cent_PR[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_mid_cent_NP_pb[6]->GetBinContent(i);
		TnP_pp = h_mid_cent_NP_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_mid_cent_NP[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_fwd_cent_PR_pb[6]->GetBinContent(i);
		TnP_pp = h_fwd_cent_PR_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_fwd_cent_PR[6]->SetBinContent(i, TnP_tot);
		
		TnP_pb = h_fwd_cent_NP_pb[6]->GetBinContent(i);
		TnP_pp = h_fwd_cent_NP_pp[6]->GetBinContent(1); // pp input has one TnP bin
		TnP_tot = TMath::Sqrt(TMath::Power(TnP_pb,2) + TMath::Power(TnP_pp,2) ); 
		h_fwd_cent_NP[6]->SetBinContent(i, TnP_tot);


	}
	

	TCanvas *c_mid_pt_PR = new TCanvas("c_mid_pt_PR","",700,700);
	TCanvas *c_mid_pt_PR_pp = new TCanvas("c_mid_pt_PR_pp","",700,700);
	TCanvas *c_mid_pt_PR_pb = new TCanvas("c_mid_pt_PR_pb","",700,700);
	TCanvas *c_mid_pt_NP = new TCanvas("c_mid_pt_NP","",700,700);
	TCanvas *c_mid_pt_NP_pp = new TCanvas("c_mid_pt_NP_pp","",700,700);
	TCanvas *c_mid_pt_NP_pb = new TCanvas("c_mid_pt_NP_pb","",700,700);
	TCanvas *c_fwd_pt_PR = new TCanvas("c_fwd_pt_PR","",700,700);
	TCanvas *c_fwd_pt_PR_pp = new TCanvas("c_fwd_pt_PR_pp","",700,700);
	TCanvas *c_fwd_pt_PR_pb = new TCanvas("c_fwd_pt_PR_pb","",700,700);
	TCanvas *c_fwd_pt_NP = new TCanvas("c_fwd_pt_NP","",700,700);
	TCanvas *c_fwd_pt_NP_pp = new TCanvas("c_fwd_pt_NP_pp","",700,700);
	TCanvas *c_fwd_pt_NP_pb = new TCanvas("c_fwd_pt_NP_pb","",700,700);
	TCanvas *c_mid_cent_PR = new TCanvas("c_mid_cent_PR","",700,700);
	TCanvas *c_mid_cent_PR_pp = new TCanvas("c_mid_cent_PR_pp","",700,700);
	TCanvas *c_mid_cent_PR_pb = new TCanvas("c_mid_cent_PR_pb","",700,700);
	TCanvas *c_mid_cent_NP = new TCanvas("c_mid_cent_NP","",700,700);
	TCanvas *c_mid_cent_NP_pp = new TCanvas("c_mid_cent_NP_pp","",700,700);
	TCanvas *c_mid_cent_NP_pb = new TCanvas("c_mid_cent_NP_pb","",700,700);
	TCanvas *c_fwd_cent_PR = new TCanvas("c_fwd_cent_PR","",700,700);
	TCanvas *c_fwd_cent_PR_pp = new TCanvas("c_fwd_cent_PR_pp","",700,700);
	TCanvas *c_fwd_cent_PR_pb = new TCanvas("c_fwd_cent_PR_pb","",700,700);
	TCanvas *c_fwd_cent_NP = new TCanvas("c_fwd_cent_NP","",700,700);
	TCanvas *c_fwd_cent_NP_pp = new TCanvas("c_fwd_cent_NP_pp","",700,700);
	TCanvas *c_fwd_cent_NP_pb = new TCanvas("c_fwd_cent_NP_pb","",700,700);

	TLegend *leg_midptPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptPR);
	TLegend *leg_midptPR_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptPR_pp);
	TLegend *leg_midptPR_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptPR_pb);
	TLegend *leg_midptNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptNP);
	TLegend *leg_midptNP_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptNP_pp);
	TLegend *leg_midptNP_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptNP_pb);
	TLegend *leg_fwdptPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptPR);
	TLegend *leg_fwdptPR_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptPR_pp);
	TLegend *leg_fwdptPR_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptPR_pb);
	TLegend *leg_fwdptNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptNP);
	TLegend *leg_fwdptNP_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptNP_pp);
	TLegend *leg_fwdptNP_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptNP_pb);
	TLegend *leg_midcentPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentPR);
	TLegend *leg_midcentPR_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentPR_pp);
	TLegend *leg_midcentPR_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentPR_pb);
	TLegend *leg_midcentNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentNP);
	TLegend *leg_midcentNP_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentNP_pp);
	TLegend *leg_midcentNP_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentNP_pb);
	TLegend *leg_fwdcentPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentPR);
	TLegend *leg_fwdcentPR_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentPR_pp);
	TLegend *leg_fwdcentPR_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentPR_pb);
	TLegend *leg_fwdcentNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentNP);
	TLegend *leg_fwdcentNP_pp = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentNP_pp);
	TLegend *leg_fwdcentNP_pb = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentNP_pb);

	float pos_x = 0.20;
	float pos_y = 0.87;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 23;

	// mid pt PR
	c_mid_pt_PR->cd();
	for (int i=0; i<nFiles-1; i++) {	
		h_mid_pt_PR[i]->GetYaxis()->SetRangeUser(0,mid_pt_PR_max);
		h_mid_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_PR[i]->SetLineWidth(2);
		h_mid_pt_PR[i]->SetLineColor(i+1);
		h_mid_pt_PR[i]->Draw("SAME");

		leg_midptPR->AddEntry(h_mid_pt_PR[i],SysName[i]);
		leg_midptPR->Draw("SAME");
		
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_PR,iPeriod,iPos);	
	c_mid_pt_PR->SaveAs("./figs/mid_pt_PR.pdf");

	c_mid_pt_PR_pp->cd();
	for (int i=0; i<nFiles-1; i++) {
		if (i==4) {
			// 4: HF
			continue;
		}
		h_mid_pt_PR_pp[i]->GetYaxis()->SetRangeUser(0,mid_pt_PR_max);
		h_mid_pt_PR_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_PR_pp[i]->SetLineWidth(2);
		h_mid_pt_PR_pp[i]->SetLineColor(i+1);
		
		h_mid_pt_PR_pp[i]->Draw("SAME");

		leg_midptPR_pp->AddEntry(h_mid_pt_PR_pp[i],SysName[i]);
		leg_midptPR_pp->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_PR_pp,iPeriod,iPos);	
	c_mid_pt_PR_pp->SaveAs("./figs/mid_pt_PR_pp.pdf");

	c_mid_pt_PR_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_pt_PR_pb[i]->GetYaxis()->SetRangeUser(0,mid_pt_PR_max);
		h_mid_pt_PR_pb[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_PR_pb[i]->SetLineWidth(2);
		h_mid_pt_PR_pb[i]->SetLineColor(i+1);
		h_mid_pt_PR_pb[i]->Draw("SAME");

		leg_midptPR_pb->AddEntry(h_mid_pt_PR_pb[i],SysName[i]);
		leg_midptPR_pb->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_PR_pb,iPeriod,iPos);	
	c_mid_pt_PR_pb->SaveAs("./figs/mid_pt_PR_pb.pdf");

	// mid pt NP
	c_mid_pt_NP->cd();
	for (int i=0; i<nFiles-1; i++) {	
		h_mid_pt_NP[i]->GetYaxis()->SetRangeUser(0,mid_pt_NP_max);
		h_mid_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_NP[i]->SetLineWidth(2);
		h_mid_pt_NP[i]->SetLineColor(i+1);
		h_mid_pt_NP[i]->Draw("SAME");

		leg_midptPR->AddEntry(h_mid_pt_NP[i],SysName[i]);
		leg_midptPR->Draw("SAME");
		
	}
	drawText("Nonprompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_NP,iPeriod,iPos);	
	c_mid_pt_NP->SaveAs("./figs/mid_pt_NP.pdf");

	c_mid_pt_NP_pp->cd();
	for (int i=0; i<nFiles-1; i++) {
		if (i==4) {
			// 4: HF
			continue;
		}
		h_mid_pt_NP_pp[i]->GetYaxis()->SetRangeUser(0,mid_pt_NP_max);
		h_mid_pt_NP_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_NP_pp[i]->SetLineWidth(2);
		h_mid_pt_NP_pp[i]->SetLineColor(i+1);
		h_mid_pt_NP_pp[i]->Draw("SAME");
		
		leg_midptNP_pp->AddEntry(h_mid_pt_NP_pp[i],SysName[i]);
		leg_midptNP_pp->Draw("SAME");
		
	}
	c_mid_pt_NP_pp->SaveAs("./figs/mid_pt_NP_pp.pdf");
	exit(1);
	drawText("NonPompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_NP_pp,iPeriod,iPos);	
	c_mid_pt_NP_pp->SaveAs("./figs/mid_pt_NP_pp.pdf");
	

	c_mid_pt_NP_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_pt_NP_pb[i]->GetYaxis()->SetRangeUser(0,mid_pt_NP_max);
		h_mid_pt_NP_pb[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_NP_pb[i]->SetLineWidth(2);
		h_mid_pt_NP_pb[i]->SetLineColor(i+1);
		h_mid_pt_NP_pb[i]->Draw("SAME");

		leg_midptNP_pb->AddEntry(h_mid_pt_NP_pb[i],SysName[i]);
		leg_midptNP_pb->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_NP_pb,iPeriod,iPos);	
	c_mid_pt_NP_pb->SaveAs("./figs/mid_pt_NP_pb.pdf");


	exit(1);


	// fwd pt PR
	c_fwd_pt_PR->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_PR[i]->GetYaxis()->SetRangeUser(0,fwd_pt_PR_max);
		h_fwd_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_PR[i]->SetLineWidth(2);
		h_fwd_pt_PR[i]->SetLineColor(i+1);
		h_fwd_pt_PR[i]->Draw("SAME");

		leg_fwdptPR->AddEntry(h_fwd_pt_PR[i],SysName[i]);
		leg_fwdptPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_PR,iPeriod,iPos);	
	c_fwd_pt_PR->SaveAs("./figs/fwd_pt_PR.pdf");

	c_fwd_pt_PR_pp->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_PR_pp[i]->GetYaxis()->SetRangeUser(0,fwd_pt_PR_max);
		h_fwd_pt_PR_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_PR_pp[i]->SetLineWidth(2);
		h_fwd_pt_PR_pp[i]->SetLineColor(i+1);
		h_fwd_pt_PR_pp[i]->Draw("SAME");

		leg_fwdptPR->AddEntry(h_fwd_pt_PR_pp[i],SysName[i]);
		leg_fwdptPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_PR_pp,iPeriod,iPos);	
	c_fwd_pt_PR_pp->SaveAs("./figs/fwd_pt_PR_pp.pdf");

	c_fwd_pt_PR_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_PR_pb[i]->GetYaxis()->SetRangeUser(0,fwd_pt_PR_max);
		h_fwd_pt_PR_pb[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_PR_pb[i]->SetLineWidth(2);
		h_fwd_pt_PR_pb[i]->SetLineColor(i+1);
		h_fwd_pt_PR_pb[i]->Draw("SAME");

		leg_fwdptPR->AddEntry(h_fwd_pt_PR_pb[i],SysName[i]);
		leg_fwdptPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_PR_pb,iPeriod,iPos);	
	c_fwd_pt_PR_pb->SaveAs("./figs/fwd_pt_PR_pb.pdf");

	// fwd pt NP
	c_fwd_pt_NP->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_NP[i]->GetYaxis()->SetRangeUser(0,fwd_pt_NP_max);
		h_fwd_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_NP[i]->SetLineWidth(2);
		h_fwd_pt_NP[i]->SetLineColor(i+1);
		h_fwd_pt_NP[i]->Draw("SAME");

		leg_fwdptNP->AddEntry(h_fwd_pt_NP[i],SysName[i]);
		leg_fwdptNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_NP,iPeriod,iPos);	
	c_fwd_pt_NP->SaveAs("./figs/fwd_pt_NP.pdf");

	c_fwd_pt_NP_pp->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_NP_pp[i]->GetYaxis()->SetRangeUser(0,fwd_pt_NP_max);
		h_fwd_pt_NP_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_NP_pp[i]->SetLineWidth(2);
		h_fwd_pt_NP_pp[i]->SetLineColor(i+1);
		h_fwd_pt_NP_pp[i]->Draw("SAME");

		leg_fwdptNP->AddEntry(h_fwd_pt_NP_pp[i],SysName[i]);
		leg_fwdptNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_NP_pp,iPeriod,iPos);	
	c_fwd_pt_NP_pp->SaveAs("./figs/fwd_pt_NP_pp.pdf");

	c_fwd_pt_NP_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_pt_NP_pb[i]->GetYaxis()->SetRangeUser(0,fwd_pt_NP_max);
		h_fwd_pt_NP_pb[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_NP_pb[i]->SetLineWidth(2);
		h_fwd_pt_NP_pb[i]->SetLineColor(i+1);
		h_fwd_pt_NP_pb[i]->Draw("SAME");

		leg_fwdptNP->AddEntry(h_fwd_pt_NP_pb[i],SysName[i]);
		leg_fwdptNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_NP_pb,iPeriod,iPos);	
	c_fwd_pt_NP_pb->SaveAs("./figs/fwd_pt_NP_pb.pdf");

	// mid cent PR
	c_mid_cent_PR->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_cent_PR[i]->GetYaxis()->SetRangeUser(0,mid_cent_PR_max);
		h_mid_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_PR[i]->SetLineWidth(2);
		h_mid_cent_PR[i]->SetLineColor(i+1);
		h_mid_cent_PR[i]->Draw("SAME");

		leg_midcentPR->AddEntry(h_mid_cent_PR[i],SysName[i]);
		leg_midcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_PR,iPeriod,iPos);	
	c_mid_cent_PR->SaveAs("./figs/mid_cent_PR.pdf");

	c_mid_cent_PR_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_cent_PR_pb[i]->GetYaxis()->SetRangeUser(0,mid_cent_PR_max);
		h_mid_cent_PR_pb[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_PR_pb[i]->SetLineWidth(2);
		h_mid_cent_PR_pb[i]->SetLineColor(i+1);
		h_mid_cent_PR_pb[i]->Draw("SAME");

		leg_midcentPR->AddEntry(h_mid_cent_PR_pb[i],SysName[i]);
		leg_midcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_PR_pb,iPeriod,iPos);	
	c_mid_cent_PR_pb->SaveAs("./figs/mid_cent_PR_pb.pdf");

	// mid cent NP
	c_mid_cent_NP->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_cent_NP[i]->GetYaxis()->SetRangeUser(0,mid_cent_NP_max);
		h_mid_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_NP[i]->SetLineWidth(2);
		h_mid_cent_NP[i]->SetLineColor(i+1);
		h_mid_cent_NP[i]->Draw("SAME");

		leg_midcentNP->AddEntry(h_mid_cent_NP[i],SysName[i]);
		leg_midcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_NP,iPeriod,iPos);	
	c_mid_cent_NP->SaveAs("./figs/mid_cent_NP.pdf");

	c_mid_cent_NP_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_mid_cent_NP_pb[i]->GetYaxis()->SetRangeUser(0,mid_cent_NP_max);
		h_mid_cent_NP_pb[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_NP_pb[i]->SetLineWidth(2);
		h_mid_cent_NP_pb[i]->SetLineColor(i+1);
		h_mid_cent_NP_pb[i]->Draw("SAME");

		leg_midcentNP->AddEntry(h_mid_cent_NP_pb[i],SysName[i]);
		leg_midcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_NP_pb,iPeriod,iPos);	
	c_mid_cent_NP_pb->SaveAs("./figs/mid_cent_NP_pb.pdf");

	// fwd cent PR
	c_fwd_cent_PR->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_cent_PR[i]->GetYaxis()->SetRangeUser(0,fwd_cent_PR_max);
		h_fwd_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_PR[i]->SetLineWidth(2);
		h_fwd_cent_PR[i]->SetLineColor(i+1);
		h_fwd_cent_PR[i]->Draw("SAME");

		leg_fwdcentPR->AddEntry(h_fwd_cent_PR[i],SysName[i]);
		leg_fwdcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_PR,iPeriod,iPos);	
	c_fwd_cent_PR->SaveAs("./figs/fwd_cent_PR.pdf");

	c_fwd_cent_PR_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_cent_PR_pb[i]->GetYaxis()->SetRangeUser(0,fwd_cent_PR_max);
		h_fwd_cent_PR_pb[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_PR_pb[i]->SetLineWidth(2);
		h_fwd_cent_PR_pb[i]->SetLineColor(i+1);
		h_fwd_cent_PR_pb[i]->Draw("SAME");

		leg_fwdcentPR->AddEntry(h_fwd_cent_PR_pb[i],SysName[i]);
		leg_fwdcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_PR_pb,iPeriod,iPos);	
	c_fwd_cent_PR_pb->SaveAs("./figs/fwd_cent_PR_pb.pdf");

	// fwd cent NP
	c_fwd_cent_NP->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_cent_NP[i]->GetYaxis()->SetRangeUser(0,fwd_cent_NP_max);
		h_fwd_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_NP[i]->SetLineWidth(2);
		h_fwd_cent_NP[i]->SetLineColor(i+1);
		h_fwd_cent_NP[i]->Draw("SAME");

		leg_fwdcentNP->AddEntry(h_fwd_cent_NP[i],SysName[i]);
		leg_fwdcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_NP,iPeriod,iPos);	
	c_fwd_cent_NP->SaveAs("./figs/fwd_cent_NP.pdf");

	c_fwd_cent_NP_pb->cd();
	for (int i=0; i<nFiles-1; i++) {
		h_fwd_cent_NP_pb[i]->GetYaxis()->SetRangeUser(0,fwd_cent_NP_max);
		h_fwd_cent_NP_pb[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_NP_pb[i]->SetLineWidth(2);
		h_fwd_cent_NP_pb[i]->SetLineColor(i+1);
		h_fwd_cent_NP_pb[i]->Draw("SAME");

		leg_fwdcentNP->AddEntry(h_fwd_cent_NP_pb[i],SysName[i]);
		leg_fwdcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_NP_pb,iPeriod,iPos);	
	c_fwd_cent_NP_pb->SaveAs("./figs/fwd_cent_NP_pb.pdf");
}