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
using TMath::Power; using TMath::Sqrt;

void compute_total_syst()
{
    
    // Start part 1-1: mid pT bins
    // Open files - pt
    auto in_sigPDF_pt = TFile::Open("syst_roots/syst_pt_sigPDF.root");
    auto in_sigPAR_pt = TFile::Open("syst_roots/syst_pt_sigPAR.root");
    auto in_bkgPDF_pt = TFile::Open("syst_roots/syst_pt_bkgPDF.root");
    auto in_bFrac_pt = TFile::Open("syst_roots/syst_pt_bFrac.root");
    //auto in_acc_pt = TFile::Open("syst_roots/syst_pt_acc.root");
    //auto in_eff_pt = TFile::Open("syst_roots/syst_pt_eff.root");
    auto in_HF_pt = TFile::Open("syst_roots/syst_pt_HF.root");


    // Open files - cent
    auto in_sigPDF_cent = TFile::Open("syst_roots/syst_cent_sigPDF.root");
    auto in_sigPAR_cent = TFile::Open("syst_roots/syst_cent_sigPAR.root");
    auto in_bkgPDF_cent = TFile::Open("syst_roots/syst_cent_bkgPDF.root");
    auto in_bFrac_cent = TFile::Open("syst_roots/syst_cent_bFrac.root");
    //auto in_acc_cent = TFile::Open("syst_roots/syst_cent_acc.root");
    //auto in_eff_cent = TFile::Open("syst_roots/syst_cent_eff.root");
    auto in_HF_cent = TFile::Open("syst_roots/syst_cent_HF.root");
    
    
    TFile out_file("syst_roots/total_syst.root", "recreate");
    
    // Assign hists
    auto h_sigPDF_mid_pt_PR = (TH1D*) in_sigPDF_pt->Get("mid_PR");
    auto h_sigPAR_mid_pt_PR = (TH1D*) in_sigPAR_pt->Get("mid_PR");
    auto h_bkgPDF_mid_pt_PR = (TH1D*) in_bkgPDF_pt->Get("mid_PR");
    auto h_bFrac_mid_pt_PR = (TH1D*) in_bFrac_pt->Get("mid_PR");
    //auto h_acc_mid_pt_PR = (TH1D*) in_acc_pt->Get("mid_PR");
    //auto h_eff_mid_pt_PR = (TH1D*) in_eff_pt->Get("mid_PR");
    auto h_HF_mid_pt_PR = (TH1D*) in_HF_pt->Get("mid_PR");

    // NP hist
    auto h_sigPDF_mid_pt_NP = (TH1D*) in_sigPDF_pt->Get("mid_NP");
    auto h_sigPAR_mid_pt_NP = (TH1D*) in_sigPAR_pt->Get("mid_NP");
    auto h_bkgPDF_mid_pt_NP = (TH1D*) in_bkgPDF_pt->Get("mid_NP");
    auto h_bFrac_mid_pt_NP = (TH1D*) in_bFrac_pt->Get("mid_NP");
    //auto h_acc_mid_pt_NP = (TH1D*) in_acc_pt->Get("mid_NP");
    //auto h_eff_mid_pt_NP = (TH1D*) in_eff_pt->Get("mid_NP");
    auto h_HF_mid_pt_NP = (TH1D*) in_HF_pt->Get("mid_NP");

    const int NBINS_mid_pt = 6;
    //double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 30, 50};
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 50};
    auto h_total_mid_pt_PR = new TH1D("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP = new TH1D("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);

    
    for (int hist_idx = 1; hist_idx < NBINS_mid_pt+1; hist_idx++) {
        double tot_syst_PR = 0;
		double acc_PR = 0; double eff_PR = 0;
		double acc_NP = 0; double eff_NP = 0;
		//double bFrac_PR = 0; double bFrac_NP = 0;
        double sigPDF_PR = h_sigPDF_mid_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_pt_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_mid_pt_PR->GetBinContent(hist_idx);
        //double acc_PR = h_acc_mid_pt_PR->GetBinContent(hist_idx);
        //double eff_PR = h_eff_mid_pt_PR->GetBinContent(hist_idx);
        double HF_PR = h_HF_mid_pt_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2) + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) );
        h_total_mid_pt_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_pt_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_mid_pt_NP->GetBinContent(hist_idx);
        //double acc_NP = h_acc_mid_pt_NP->GetBinContent(hist_idx);
        //double eff_NP = h_eff_mid_pt_NP->GetBinContent(hist_idx);
        double HF_NP = h_HF_mid_pt_NP->GetBinContent(hist_idx);

        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                        + Power(bFrac_NP,2) + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) );
        h_total_mid_pt_NP->SetBinContent(hist_idx, tot_syst_NP);
        //cout << h_total_mid_pt_NP->GetBinContent(hist_idx) << endl;
    }


    // Start part 1-2: fwd pT bins
    // Assign hists
    auto h_sigPDF_fwd_pt_PR = (TH1D*) in_sigPDF_pt->Get("fwd_PR");
    auto h_sigPAR_fwd_pt_PR = (TH1D*) in_sigPAR_pt->Get("fwd_PR");
    auto h_bkgPDF_fwd_pt_PR = (TH1D*) in_bkgPDF_pt->Get("fwd_PR");
    auto h_bFrac_fwd_pt_PR = (TH1D*) in_bFrac_pt->Get("fwd_PR");
    //auto h_acc_fwd_pt_PR = (TH1D*) in_acc_pt->Get("fwd_PR");
    //auto h_eff_fwd_pt_PR = (TH1D*) in_eff_pt->Get("fwd_PR");
    auto h_HF_fwd_pt_PR = (TH1D*) in_HF_pt->Get("fwd_PR");

    // NP hist
    auto h_sigPDF_fwd_pt_NP = (TH1D*) in_sigPDF_pt->Get("fwd_NP");
    auto h_sigPAR_fwd_pt_NP = (TH1D*) in_sigPAR_pt->Get("fwd_NP");
    auto h_bkgPDF_fwd_pt_NP = (TH1D*) in_bkgPDF_pt->Get("fwd_NP");
    auto h_bFrac_fwd_pt_NP = (TH1D*) in_bFrac_pt->Get("fwd_NP");
    //auto h_acc_fwd_pt_NP = (TH1D*) in_acc_pt->Get("fwd_NP");
    //auto h_eff_fwd_pt_NP = (TH1D*) in_eff_pt->Get("fwd_NP");
    auto h_HF_fwd_pt_NP = (TH1D*) in_HF_pt->Get("fwd_NP");

    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 5, 6.5, 12, 50};
    auto h_total_fwd_pt_PR = new TH1D("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP = new TH1D("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);

    
    for (int hist_idx = 1; hist_idx < NBINS_fwd_pt+1; hist_idx++) {
        double tot_syst_PR = 0;
		double acc_PR = 0; double eff_PR = 0;
		double acc_NP = 0; double eff_NP = 0;
		//double bFrac_PR = 0; double bFrac_NP = 0;
        double sigPDF_PR = h_sigPDF_fwd_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_pt_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_fwd_pt_PR->GetBinContent(hist_idx);
        //double acc_PR = h_acc_fwd_pt_PR->GetBinContent(hist_idx);
        //double eff_PR = h_eff_fwd_pt_PR->GetBinContent(hist_idx);
        double HF_PR = h_HF_fwd_pt_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2) + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) );
        h_total_fwd_pt_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_pt_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_fwd_pt_NP->GetBinContent(hist_idx);
        //double acc_NP = h_acc_fwd_pt_NP->GetBinContent(hist_idx);
        //double eff_NP = h_eff_fwd_pt_NP->GetBinContent(hist_idx);
        double HF_NP = h_HF_fwd_pt_NP->GetBinContent(hist_idx);

        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                        + Power(bFrac_NP,2) + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) );
        h_total_fwd_pt_NP->SetBinContent(hist_idx, tot_syst_NP);
        //cout << h_total_mid_pt_NP->GetBinContent(hist_idx) << endl;
    }


    // Start part 2-1: mid cent bins
    // Assign hists
    // Assign hists
    auto h_sigPDF_mid_cent_PR = (TH1D*) in_sigPDF_cent->Get("mid_PR");
    auto h_sigPAR_mid_cent_PR = (TH1D*) in_sigPAR_cent->Get("mid_PR");
    auto h_bkgPDF_mid_cent_PR = (TH1D*) in_bkgPDF_cent->Get("mid_PR");
    auto h_bFrac_mid_cent_PR = (TH1D*) in_bFrac_cent->Get("mid_PR");
    //auto h_acc_mid_cent_PR = (TH1D*) in_acc_cent->Get("mid_PR");
    //auto h_eff_mid_cent_PR = (TH1D*) in_eff_cent->Get("mid_PR");
    auto h_HF_mid_cent_PR = (TH1D*) in_HF_cent->Get("mid_PR");

    // NP hist
    auto h_sigPDF_mid_cent_NP = (TH1D*) in_sigPDF_cent->Get("mid_NP");
    auto h_sigPAR_mid_cent_NP = (TH1D*) in_sigPAR_cent->Get("mid_NP");
    auto h_bkgPDF_mid_cent_NP = (TH1D*) in_bkgPDF_cent->Get("mid_NP");
    auto h_bFrac_mid_cent_NP = (TH1D*) in_bFrac_cent->Get("mid_NP");
    //auto h_acc_mid_cent_NP = (TH1D*) in_acc_cent->Get("mid_NP");
    //auto h_eff_mid_cent_NP = (TH1D*) in_eff_cent->Get("mid_NP");
    auto h_HF_mid_cent_NP = (TH1D*) in_HF_cent->Get("mid_NP");

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    auto h_total_mid_cent_PR = new TH1D("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP = new TH1D("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);

    for (int hist_idx = 1; hist_idx < NBINS_mid_cent+1; hist_idx++) {
        double tot_syst_PR = 0;
		double acc_PR = 0; double eff_PR = 0;
		double acc_NP = 0; double eff_NP = 0;
		//double bFrac_PR = 0; double bFrac_NP = 0;
        double sigPDF_PR = h_sigPDF_mid_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_cent_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_mid_cent_PR->GetBinContent(hist_idx);
        //double acc_PR = h_acc_mid_cent_PR->GetBinContent(hist_idx);
        //double eff_PR = h_eff_mid_cent_PR->GetBinContent(hist_idx);
        double HF_PR = h_HF_mid_cent_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2) + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) );
        h_total_mid_cent_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_cent_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_mid_cent_NP->GetBinContent(hist_idx);
        //double acc_NP = h_acc_mid_cent_NP->GetBinContent(hist_idx);
        //double eff_NP = h_eff_mid_cent_NP->GetBinContent(hist_idx);
        double HF_NP = h_HF_mid_cent_NP->GetBinContent(hist_idx);

        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                        + Power(bFrac_NP,2) + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) );
        h_total_mid_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
        //cout << h_total_mid_cent_NP->GetBinContent(hist_idx) << endl;
    }


    // Start part 2-2: fwd cent bins
    // Assign hists
    auto h_sigPDF_fwd_cent_PR = (TH1D*) in_sigPDF_cent->Get("fwd_PR");
    auto h_sigPAR_fwd_cent_PR = (TH1D*) in_sigPAR_cent->Get("fwd_PR");
    auto h_bkgPDF_fwd_cent_PR = (TH1D*) in_bkgPDF_cent->Get("fwd_PR");
    auto h_bFrac_fwd_cent_PR = (TH1D*) in_bFrac_cent->Get("fwd_PR");
    //auto h_acc_fwd_cent_PR = (TH1D*) in_acc_cent->Get("fwd_PR");
    //auto h_eff_fwd_cent_PR = (TH1D*) in_eff_cent->Get("fwd_PR");
    auto h_HF_fwd_cent_PR = (TH1D*) in_HF_cent->Get("fwd_PR");

    // NP hist
    auto h_sigPDF_fwd_cent_NP = (TH1D*) in_sigPDF_cent->Get("fwd_NP");
    auto h_sigPAR_fwd_cent_NP = (TH1D*) in_sigPAR_cent->Get("fwd_NP");
    auto h_bkgPDF_fwd_cent_NP = (TH1D*) in_bkgPDF_cent->Get("fwd_NP");
    auto h_bFrac_fwd_cent_NP = (TH1D*) in_bFrac_cent->Get("fwd_NP");
    //auto h_acc_fwd_cent_NP = (TH1D*) in_acc_cent->Get("fwd_NP");
    //auto h_eff_fwd_cent_NP = (TH1D*) in_eff_cent->Get("fwd_NP");
    auto h_HF_fwd_cent_NP = (TH1D*) in_HF_cent->Get("fwd_NP");

    const int NBINS_fwd_cent = 6;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    auto h_total_fwd_cent_PR = new TH1D("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP = new TH1D("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    
    for (int hist_idx = 1; hist_idx < NBINS_fwd_cent+1; hist_idx++) {
        double tot_syst_PR = 0;
		double acc_PR = 0; double eff_PR = 0;
		double acc_NP = 0; double eff_NP = 0;
		//double bFrac_PR = 0; double bFrac_NP = 0;
        double sigPDF_PR = h_sigPDF_fwd_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_cent_PR->GetBinContent(hist_idx);
        double bFrac_PR = h_bFrac_fwd_cent_PR->GetBinContent(hist_idx);
        //double acc_PR = h_acc_fwd_cent_PR->GetBinContent(hist_idx);
        //double eff_PR = h_eff_fwd_cent_PR->GetBinContent(hist_idx);
        double HF_PR = h_HF_fwd_cent_PR->GetBinContent(hist_idx);
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                + Power(bFrac_PR,2) + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) );
        h_total_fwd_cent_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_cent_NP->GetBinContent(hist_idx);
        double bFrac_NP = h_bFrac_fwd_cent_NP->GetBinContent(hist_idx);
        //double acc_NP = h_acc_fwd_cent_NP->GetBinContent(hist_idx);
        //double eff_NP = h_eff_fwd_cent_NP->GetBinContent(hist_idx);
        double HF_NP = h_HF_fwd_cent_NP->GetBinContent(hist_idx);

        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                        + Power(bFrac_NP,2) + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) );
        h_total_fwd_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
        //cout << h_total_fwd_cent_NP->GetBinContent(hist_idx) << endl;
    }
    out_file.cd();
    out_file.Write();
    out_file.Close();
}
