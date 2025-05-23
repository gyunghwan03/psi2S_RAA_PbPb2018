// ========== Start Coment ========== //
// You need inputs:
    // sigPDF
    // sigPAR
    // bkgPDF
    
    // ctauERR
    // ctauRES
    // ctauBKG
    // ctauTRUE

    // Corrections

// This code consits of :
// Part 1: pT bins
    // 1 - 1 mid_PR, mid_NP
    // 1 - 2 fwd_PR, fwd_NP
// Part 2: cent. bins
    // 2 - 1 mid_PR, mid_NP
    // 2 - 2 fwd_PR, fwd_NP
    
// You will get
    // One root file
// ========== End of Coment ========== //


#include "TMath.h"
#include "TFile.h"
#include "TH1.h"

using TMath::Power; using TMath::Sqrt;

void compute_DoubleRatio_syst()
{
    
    // Start part 1-1: mid pT bins
    // Open files - pt
    auto in_sigPDF_pt = TFile::Open("syst_roots/syst_pt_sigPDF.root");
    auto in_sigPAR_pt = TFile::Open("syst_roots/syst_pt_sigPAR.root");
    auto in_bkgPDF_pt = TFile::Open("syst_roots/syst_pt_bkgPDF.root");
    auto in_bFrac_pt = TFile::Open("syst_roots/syst_pt_bFrac.root");


    // Open files - cent
    auto in_sigPDF_cent = TFile::Open("syst_roots/syst_cent_sigPDF.root");
    auto in_sigPAR_cent = TFile::Open("syst_roots/syst_cent_sigPAR.root");
    auto in_bkgPDF_cent = TFile::Open("syst_roots/syst_cent_bkgPDF.root");
    auto in_bFrac_cent = TFile::Open("syst_roots/syst_cent_bFrac.root");
    
    
    TFile out_file("syst_roots/DoubleRatio_syst.root", "recreate");
    
    // Assign hists
    auto h_sigPDF_mid_pt_PR = (TH1D*) in_sigPDF_pt->Get("mid_PR");
    auto h_sigPDF_mid_pt_PR_pp = (TH1D*) in_sigPDF_pt->Get("mid_PR_pp");
    auto h_sigPDF_mid_pt_PR_pb = (TH1D*) in_sigPDF_pt->Get("mid_PR_pb");
    auto h_sigPAR_mid_pt_PR = (TH1D*) in_sigPAR_pt->Get("mid_PR");
    auto h_sigPAR_mid_pt_PR_pp = (TH1D*) in_sigPAR_pt->Get("mid_PR_pp");
    auto h_sigPAR_mid_pt_PR_pb = (TH1D*) in_sigPAR_pt->Get("mid_PR_pb");
    auto h_bkgPDF_mid_pt_PR = (TH1D*) in_bkgPDF_pt->Get("mid_PR");
    auto h_bkgPDF_mid_pt_PR_pp = (TH1D*) in_bkgPDF_pt->Get("mid_PR_pp");
    auto h_bkgPDF_mid_pt_PR_pb = (TH1D*) in_bkgPDF_pt->Get("mid_PR_pb");
    auto h_bFrac_mid_pt_PR = (TH1D*) in_bFrac_pt->Get("mid_PR");
    auto h_bFrac_mid_pt_PR_pp = (TH1D*) in_bFrac_pt->Get("mid_PR_pp");
    auto h_bFrac_mid_pt_PR_pb = (TH1D*) in_bFrac_pt->Get("mid_PR_pb");

    // NP hist
    auto h_sigPDF_mid_pt_NP = (TH1D*) in_sigPDF_pt->Get("mid_NP");
    auto h_sigPDF_mid_pt_NP_pp = (TH1D*) in_sigPDF_pt->Get("mid_NP_pp");
    auto h_sigPDF_mid_pt_NP_pb = (TH1D*) in_sigPDF_pt->Get("mid_NP_pb");
    auto h_sigPAR_mid_pt_NP = (TH1D*) in_sigPAR_pt->Get("mid_NP");
    auto h_sigPAR_mid_pt_NP_pp = (TH1D*) in_sigPAR_pt->Get("mid_NP_pp");
    auto h_sigPAR_mid_pt_NP_pb = (TH1D*) in_sigPAR_pt->Get("mid_NP_pb");
    auto h_bkgPDF_mid_pt_NP = (TH1D*) in_bkgPDF_pt->Get("mid_NP");
    auto h_bkgPDF_mid_pt_NP_pp = (TH1D*) in_bkgPDF_pt->Get("mid_NP_pp");
    auto h_bkgPDF_mid_pt_NP_pb = (TH1D*) in_bkgPDF_pt->Get("mid_NP_pb");
    auto h_bFrac_mid_pt_NP = (TH1D*) in_bFrac_pt->Get("mid_NP");
    auto h_bFrac_mid_pt_NP_pp = (TH1D*) in_bFrac_pt->Get("mid_NP_pp");
    auto h_bFrac_mid_pt_NP_pb = (TH1D*) in_bFrac_pt->Get("mid_NP_pb");

    const int NBINS_mid_pt = 6;
    //double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 30, 50};
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    auto h_total_mid_pt_PR = new TH1D("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_PR_pp = new TH1D("mid_pt_PR_pp", "mid_PR_pp", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_PR_pb = new TH1D("mid_pt_PR_pb", "mid_PR_pb", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP = new TH1D("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP_pp = new TH1D("mid_pt_NP_pp", "mid_NP_pp", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP_pb = new TH1D("mid_pt_NP_pb", "mid_NP_pb", NBINS_mid_pt, edges_mid_pt);

    
    for (int hist_idx = 1; hist_idx < NBINS_mid_pt+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_mid_pt_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_mid_pt_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_mid_pt_PR_pp->GetBinContent(hist_idx);
        double bFrac_PR_pp = h_bFrac_mid_pt_PR_pp->GetBinContent(hist_idx);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2) + Power(bFrac_PR_pp,2));
        h_total_mid_pt_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);
		cout << "################ pp PR ################" << endl;
		cout << "mid pt PR : " << tot_syst_PR_pp << endl;

        double tot_syst_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_mid_pt_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_mid_pt_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_mid_pt_PR_pb->GetBinContent(hist_idx);
        double bFrac_PR_pb = h_bFrac_mid_pt_PR_pb->GetBinContent(hist_idx);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2) + Power(bFrac_PR_pb,2));
        h_total_mid_pt_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
        double sigPDF_PR = h_sigPDF_mid_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_pt_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_mid_pt_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2) + Power(bFrac_PR,2));
        h_total_mid_pt_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_mid_pt_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_mid_pt_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_mid_pt_NP_pp->GetBinContent(hist_idx);
        double bFrac_NP_pp = h_bFrac_mid_pt_NP_pp->GetBinContent(hist_idx);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2) + Power(bFrac_NP_pp,2));
        h_total_mid_pt_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_mid_pt_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_mid_pt_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_mid_pt_NP_pb->GetBinContent(hist_idx);
        double bFrac_NP_pb = h_bFrac_mid_pt_NP_pb->GetBinContent(hist_idx);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                + Power(bFrac_NP_pb,2));
        h_total_mid_pt_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        //cout << h_total_mid_pt_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_pt_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_mid_pt_NP->GetBinContent(hist_idx);
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                + Power(bFrac_NP,2));
        h_total_mid_pt_NP->SetBinContent(hist_idx, tot_syst_NP);
    }


    // Start part 1-2: fwd pT bins
    // Assign hists
    auto h_sigPDF_fwd_pt_PR = (TH1D*) in_sigPDF_pt->Get("fwd_PR");
    auto h_sigPDF_fwd_pt_PR_pp = (TH1D*) in_sigPDF_pt->Get("fwd_PR_pp");
    auto h_sigPDF_fwd_pt_PR_pb = (TH1D*) in_sigPDF_pt->Get("fwd_PR_pb");
    auto h_sigPAR_fwd_pt_PR = (TH1D*) in_sigPAR_pt->Get("fwd_PR");
    auto h_sigPAR_fwd_pt_PR_pp = (TH1D*) in_sigPAR_pt->Get("fwd_PR_pp");
    auto h_sigPAR_fwd_pt_PR_pb = (TH1D*) in_sigPAR_pt->Get("fwd_PR_pb");
    auto h_bkgPDF_fwd_pt_PR = (TH1D*) in_bkgPDF_pt->Get("fwd_PR");
    auto h_bkgPDF_fwd_pt_PR_pp = (TH1D*) in_bkgPDF_pt->Get("fwd_PR_pp");
    auto h_bkgPDF_fwd_pt_PR_pb = (TH1D*) in_bkgPDF_pt->Get("fwd_PR_pb");
    auto h_bFrac_fwd_pt_PR = (TH1D*) in_bFrac_pt->Get("fwd_PR");
    auto h_bFrac_fwd_pt_PR_pp = (TH1D*) in_bFrac_pt->Get("fwd_PR_pp");
    auto h_bFrac_fwd_pt_PR_pb = (TH1D*) in_bFrac_pt->Get("fwd_PR_pb");

    // NP hist
    auto h_sigPDF_fwd_pt_NP = (TH1D*) in_sigPDF_pt->Get("fwd_NP");
    auto h_sigPDF_fwd_pt_NP_pp = (TH1D*) in_sigPDF_pt->Get("fwd_NP_pp");
    auto h_sigPDF_fwd_pt_NP_pb = (TH1D*) in_sigPDF_pt->Get("fwd_NP_pb");
    auto h_sigPAR_fwd_pt_NP = (TH1D*) in_sigPAR_pt->Get("fwd_NP");
    auto h_sigPAR_fwd_pt_NP_pp = (TH1D*) in_sigPAR_pt->Get("fwd_NP_pp");
    auto h_sigPAR_fwd_pt_NP_pb = (TH1D*) in_sigPAR_pt->Get("fwd_NP_pb");
    auto h_bkgPDF_fwd_pt_NP = (TH1D*) in_bkgPDF_pt->Get("fwd_NP");
    auto h_bkgPDF_fwd_pt_NP_pp = (TH1D*) in_bkgPDF_pt->Get("fwd_NP_pp");
    auto h_bkgPDF_fwd_pt_NP_pb = (TH1D*) in_bkgPDF_pt->Get("fwd_NP_pb");
    auto h_bFrac_fwd_pt_NP = (TH1D*) in_bFrac_pt->Get("fwd_NP");
    auto h_bFrac_fwd_pt_NP_pp = (TH1D*) in_bFrac_pt->Get("fwd_NP_pp");
    auto h_bFrac_fwd_pt_NP_pb = (TH1D*) in_bFrac_pt->Get("fwd_NP_pb");

    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 40};
    auto h_total_fwd_pt_PR = new TH1D("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_PR_pp = new TH1D("fwd_pt_PR_pp", "fwd_PR_pp", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_PR_pb = new TH1D("fwd_pt_PR_pb", "fwd_PR_pb", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP = new TH1D("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP_pp = new TH1D("fwd_pt_NP_pp", "fwd_NP_pp", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP_pb = new TH1D("fwd_pt_NP_pb", "fwd_NP_pb", NBINS_fwd_pt, edges_fwd_pt);

    
	cout << "Fwd pT" << endl;
    for (int hist_idx = 1; hist_idx < NBINS_fwd_pt+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double bFrac_PR_pp = h_bFrac_fwd_pt_PR_pp->GetBinContent(hist_idx);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                + Power(bFrac_PR_pp,2));
        h_total_fwd_pt_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double bFrac_PR_pb = h_bFrac_fwd_pt_PR_pb->GetBinContent(hist_idx);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                + Power(bFrac_PR_pb,2));
        h_total_fwd_pt_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
        double sigPDF_PR = h_sigPDF_fwd_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_pt_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_fwd_pt_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2));
        h_total_fwd_pt_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double bFrac_NP_pp = h_bFrac_fwd_pt_NP_pp->GetBinContent(hist_idx);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                + Power(bFrac_NP_pp,2));
        h_total_fwd_pt_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double bFrac_NP_pb = h_bFrac_fwd_pt_NP_pb->GetBinContent(hist_idx);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                + Power(bFrac_NP_pb,2));
        h_total_fwd_pt_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
		//
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_pt_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_fwd_pt_NP->GetBinContent(hist_idx);
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                + Power(bFrac_NP,2));
        h_total_fwd_pt_NP->SetBinContent(hist_idx, tot_syst_NP);
    }


    // Start part 2-1: mid cent bins
    // Assign hists
    auto h_sigPDF_mid_cent_PR = (TH1D*) in_sigPDF_cent->Get("mid_PR");
    auto h_sigPDF_mid_cent_PR_pp = (TH1D*) in_sigPDF_cent->Get("mid_PR_pp");
    auto h_sigPDF_mid_cent_PR_pb = (TH1D*) in_sigPDF_cent->Get("mid_PR_pb");
    auto h_sigPAR_mid_cent_PR = (TH1D*) in_sigPAR_cent->Get("mid_PR");
    auto h_sigPAR_mid_cent_PR_pp = (TH1D*) in_sigPAR_cent->Get("mid_PR_pp");
    auto h_sigPAR_mid_cent_PR_pb = (TH1D*) in_sigPAR_cent->Get("mid_PR_pb");
    auto h_bkgPDF_mid_cent_PR = (TH1D*) in_bkgPDF_cent->Get("mid_PR");
    auto h_bkgPDF_mid_cent_PR_pp = (TH1D*) in_bkgPDF_cent->Get("mid_PR_pp");
    auto h_bkgPDF_mid_cent_PR_pb = (TH1D*) in_bkgPDF_cent->Get("mid_PR_pb");
    auto h_bFrac_mid_cent_PR = (TH1D*) in_bFrac_cent->Get("mid_PR");
    auto h_bFrac_mid_cent_PR_pp = (TH1D*) in_bFrac_cent->Get("mid_PR_pp");
    auto h_bFrac_mid_cent_PR_pb = (TH1D*) in_bFrac_cent->Get("mid_PR_pb");

    // NP hist
    auto h_sigPDF_mid_cent_NP = (TH1D*) in_sigPDF_cent->Get("mid_NP");
    auto h_sigPDF_mid_cent_NP_pp = (TH1D*) in_sigPDF_cent->Get("mid_NP_pp");
    auto h_sigPDF_mid_cent_NP_pb = (TH1D*) in_sigPDF_cent->Get("mid_NP_pb");
    auto h_sigPAR_mid_cent_NP = (TH1D*) in_sigPAR_cent->Get("mid_NP");
    auto h_sigPAR_mid_cent_NP_pp = (TH1D*) in_sigPAR_cent->Get("mid_NP_pp");
    auto h_sigPAR_mid_cent_NP_pb = (TH1D*) in_sigPAR_cent->Get("mid_NP_pb");
    auto h_bkgPDF_mid_cent_NP = (TH1D*) in_bkgPDF_cent->Get("mid_NP");
    auto h_bkgPDF_mid_cent_NP_pp = (TH1D*) in_bkgPDF_cent->Get("mid_NP_pp");
    auto h_bkgPDF_mid_cent_NP_pb = (TH1D*) in_bkgPDF_cent->Get("mid_NP_pb");
    auto h_bFrac_mid_cent_NP = (TH1D*) in_bFrac_cent->Get("mid_NP");
    auto h_bFrac_mid_cent_NP_pp = (TH1D*) in_bFrac_cent->Get("mid_NP_pp");
    auto h_bFrac_mid_cent_NP_pb = (TH1D*) in_bFrac_cent->Get("mid_NP_pb");
    

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0,10,20,30,40,50,90};
    auto h_total_mid_cent_PR = new TH1D("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_PR_pp = new TH1D("mid_cent_PR_pp", "mid_PR_pp", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_PR_pb = new TH1D("mid_cent_PR_pb", "mid_PR_pb", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP = new TH1D("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP_pp = new TH1D("mid_cent_NP_pp", "mid_NP_pp", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP_pb = new TH1D("mid_cent_NP_pb", "mid_NP_pb", NBINS_mid_cent, edges_mid_cent);

    for (int hist_idx = 1; hist_idx < NBINS_mid_cent+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        
        double sigPDF_PR_pp = h_sigPDF_mid_cent_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_mid_cent_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_mid_cent_PR_pp->GetBinContent(hist_idx);
        
        double bFrac_PR_pp = h_bFrac_mid_cent_PR_pp->GetBinContent(hist_idx);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                + Power(bFrac_PR_pp,2));
        h_total_mid_cent_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_mid_cent_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_mid_cent_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_mid_cent_PR_pb->GetBinContent(hist_idx);
        double bFrac_PR_pb = h_bFrac_mid_cent_PR_pb->GetBinContent(hist_idx);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                + Power(bFrac_PR_pb,2));
        h_total_mid_cent_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
        double sigPDF_PR = h_sigPDF_mid_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_cent_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_mid_cent_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2));
        h_total_mid_cent_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_mid_cent_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_mid_cent_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_mid_cent_NP_pp->GetBinContent(hist_idx);
        double bFrac_NP_pp = h_bFrac_mid_cent_NP_pp->GetBinContent(hist_idx);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                + Power(bFrac_NP_pp,2));
        h_total_mid_cent_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_mid_cent_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_mid_cent_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_mid_cent_NP_pb->GetBinContent(hist_idx);
        double bFrac_NP_pb = h_bFrac_mid_cent_NP_pb->GetBinContent(hist_idx);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                + Power(bFrac_NP_pb,2));
        h_total_mid_cent_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        //cout << h_total_mid_cent_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_cent_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_mid_cent_NP->GetBinContent(hist_idx);
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                + Power(bFrac_NP,2));
        h_total_mid_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
    }


    // Start part 2-2: fwd cent bins
    // Assign hists
    auto h_sigPDF_fwd_cent_PR = (TH1D*) in_sigPDF_cent->Get("fwd_PR");
    auto h_sigPDF_fwd_cent_PR_pp = (TH1D*) in_sigPDF_cent->Get("fwd_PR_pp");
    auto h_sigPDF_fwd_cent_PR_pb = (TH1D*) in_sigPDF_cent->Get("fwd_PR_pb");
    auto h_sigPAR_fwd_cent_PR = (TH1D*) in_sigPAR_cent->Get("fwd_PR");
    auto h_sigPAR_fwd_cent_PR_pp = (TH1D*) in_sigPAR_cent->Get("fwd_PR_pp");
    auto h_sigPAR_fwd_cent_PR_pb = (TH1D*) in_sigPAR_cent->Get("fwd_PR_pb");
    auto h_bkgPDF_fwd_cent_PR = (TH1D*) in_bkgPDF_cent->Get("fwd_PR");
    auto h_bkgPDF_fwd_cent_PR_pp = (TH1D*) in_bkgPDF_cent->Get("fwd_PR_pp");
    auto h_bkgPDF_fwd_cent_PR_pb = (TH1D*) in_bkgPDF_cent->Get("fwd_PR_pb");
    auto h_bFrac_fwd_cent_PR = (TH1D*) in_bFrac_cent->Get("fwd_PR");
    auto h_bFrac_fwd_cent_PR_pp = (TH1D*) in_bFrac_cent->Get("fwd_PR_pp");
    auto h_bFrac_fwd_cent_PR_pb = (TH1D*) in_bFrac_cent->Get("fwd_PR_pb");

    // NP hist
    auto h_sigPDF_fwd_cent_NP = (TH1D*) in_sigPDF_cent->Get("fwd_NP");
    auto h_sigPDF_fwd_cent_NP_pp = (TH1D*) in_sigPDF_cent->Get("fwd_NP_pp");
    auto h_sigPDF_fwd_cent_NP_pb = (TH1D*) in_sigPDF_cent->Get("fwd_NP_pb");
    auto h_sigPAR_fwd_cent_NP = (TH1D*) in_sigPAR_cent->Get("fwd_NP");
    auto h_sigPAR_fwd_cent_NP_pp = (TH1D*) in_sigPAR_cent->Get("fwd_NP_pp");
    auto h_sigPAR_fwd_cent_NP_pb = (TH1D*) in_sigPAR_cent->Get("fwd_NP_pb");
    auto h_bkgPDF_fwd_cent_NP = (TH1D*) in_bkgPDF_cent->Get("fwd_NP");
    auto h_bkgPDF_fwd_cent_NP_pp = (TH1D*) in_bkgPDF_cent->Get("fwd_NP_pp");
    auto h_bkgPDF_fwd_cent_NP_pb = (TH1D*) in_bkgPDF_cent->Get("fwd_NP_pb");
    auto h_bFrac_fwd_cent_NP = (TH1D*) in_bFrac_cent->Get("fwd_NP");
    auto h_bFrac_fwd_cent_NP_pp = (TH1D*) in_bFrac_cent->Get("fwd_NP_pp");
    auto h_bFrac_fwd_cent_NP_pb = (TH1D*) in_bFrac_cent->Get("fwd_NP_pb");

    const int NBINS_fwd_cent = 4;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 30, 50, 90};
    auto h_total_fwd_cent_PR = new TH1D("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_PR_pp = new TH1D("fwd_cent_PR_pp", "fwd_PR_pp", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_PR_pb = new TH1D("fwd_cent_PR_pb", "fwd_PR_pb", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP = new TH1D("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP_pp = new TH1D("fwd_cent_NP_pp", "fwd_NP_pp", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP_pb = new TH1D("fwd_cent_NP_pb", "fwd_NP_pb", NBINS_fwd_cent, edges_fwd_cent);
    
	cout << "FWD CENT PR" << endl;
    for (int hist_idx = 1; hist_idx < NBINS_fwd_cent+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double bFrac_PR_pp = h_bFrac_fwd_cent_PR_pp->GetBinContent(hist_idx);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                + Power(bFrac_PR_pp,2));
        h_total_fwd_cent_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double bFrac_PR_pb = h_bFrac_fwd_cent_PR_pb->GetBinContent(hist_idx);

        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                + Power(bFrac_PR_pb,2));
        h_total_fwd_cent_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

		//printf("Cent %.f - %.f\n" ,edges_fwd_cent[hist_idx], edges_fwd_cent[hist_idx+1]);
		//cout << "sigPDF pb : " << sigPDF_PR_pb << "\tsigPR pb : " << sigPAR_PR_pb << "\tTNP_PR_pb : " << TNP_PR_pb << "\ttotal : " << tot_syst_PR_pb << endl; 

        double tot_syst_PR = 0;
        double sigPDF_PR = h_sigPDF_fwd_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_cent_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_fwd_cent_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2));
        h_total_fwd_cent_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double bFrac_NP_pp = h_bFrac_fwd_cent_NP_pp->GetBinContent(hist_idx);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                + Power(bFrac_NP_pp,2));
        h_total_fwd_cent_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double bFrac_NP_pb = h_bFrac_fwd_cent_NP_pb->GetBinContent(hist_idx);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                + Power(bFrac_NP_pb,2));
        h_total_fwd_cent_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
		printf("Cent %.f - %.f\n" ,edges_fwd_cent[hist_idx-1], edges_fwd_cent[hist_idx]);
        //cout << h_total_fwd_cent_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_cent_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_fwd_cent_NP->GetBinContent(hist_idx);
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                + Power(bFrac_NP,2));
        h_total_fwd_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
    }
    out_file.cd();
    out_file.Write();
    out_file.Close();
}
