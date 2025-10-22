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

// for ctau cut : remove the b fraction syst
// ========== End of Coment ========== //


#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>

using TMath::Power; using TMath::Sqrt;
using namespace std;

void compute_total_syst()
{
    
    // Start part 1-1: mid pT bins
    // Open files - pt
    auto in_sigPDF_pt = TFile::Open("syst_roots/syst_pt_sigPDF.root");
    auto in_sigPAR_pt = TFile::Open("syst_roots/syst_pt_sigPAR.root");
    auto in_bkgPDF_pt = TFile::Open("syst_roots/syst_pt_bkgPDF.root");
    //auto in_bFrac_pt = TFile::Open("syst_roots/syst_pt_bFrac.root");
    auto in_acc_pt = TFile::Open("syst_roots/syst_pt_acc.root");
    auto in_eff_pt = TFile::Open("syst_roots/syst_pt_Eff.root");
    auto in_HF_pt = TFile::Open("syst_roots/syst_pt_HF.root");
    auto in_TNP_pt = TFile::Open("syst_roots/syst_pt_TNP.root");
    auto in_TNPpp_pt = TFile::Open("syst_roots/syst_pt_TNP_pp.root");


    // Open files - cent
    auto in_sigPDF_cent = TFile::Open("syst_roots/syst_cent_sigPDF.root");
    auto in_sigPAR_cent = TFile::Open("syst_roots/syst_cent_sigPAR.root");
    auto in_bkgPDF_cent = TFile::Open("syst_roots/syst_cent_bkgPDF.root");
    //auto in_bFrac_cent = TFile::Open("syst_roots/syst_cent_bFrac.root");
    auto in_TNP_cent = TFile::Open("syst_roots/syst_cent_TNP.root");
    auto in_TNPpp_cent = TFile::Open("syst_roots/syst_cent_TNP_pp.root");
    auto in_acc_cent = TFile::Open("syst_roots/syst_cent_acc.root");
    auto in_eff_cent = TFile::Open("syst_roots/syst_cent_Eff.root");
    auto in_HF_cent = TFile::Open("syst_roots/syst_cent_HF.root");
    
    
    TFile out_file("syst_roots/total_syst.root", "recreate");
    
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
    //auto h_bFrac_mid_pt_PR = (TH1D*) in_bFrac_pt->Get("mid_PR");
    //auto h_bFrac_mid_pt_PR_pp = (TH1D*) in_bFrac_pt->Get("mid_PR_pp");
    //auto h_bFrac_mid_pt_PR_pb = (TH1D*) in_bFrac_pt->Get("mid_PR_pb");
    auto h_acc_mid_pt_PR = (TH1D*) in_acc_pt->Get("mid");
    auto h_acc_mid_pt_PR_pp = (TH1D*) in_acc_pt->Get("mid_pp");
    auto h_acc_mid_pt_PR_pb = (TH1D*) in_acc_pt->Get("mid_pb");
    auto h_eff_mid_pt_PR = (TH1D*) in_eff_pt->Get("mid_PR");
    auto h_eff_mid_pt_PR_pp = (TH1D*) in_eff_pt->Get("mid_PR_pp");
    auto h_eff_mid_pt_PR_pb = (TH1D*) in_eff_pt->Get("mid_PR_pb");
    auto h_HF_mid_pt_PR = (TH1D*) in_HF_pt->Get("mid_PR");
    //auto h_TNP_mid_pt_PR = (TH1D*) in_TNP_pt->Get("mid_PR");
    auto h_TNP_mid_pt_PR_pp = (TH1D*) in_TNPpp_pt->Get("mid_PR"); 
    auto h_TNP_mid_pt_PR_pb = (TH1D*) in_TNP_pt->Get("mid_PR"); 

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
    //auto h_bFrac_mid_pt_NP = (TH1D*) in_bFrac_pt->Get("mid_NP");
    //auto h_bFrac_mid_pt_NP_pp = (TH1D*) in_bFrac_pt->Get("mid_NP_pp");
    //auto h_bFrac_mid_pt_NP_pb = (TH1D*) in_bFrac_pt->Get("mid_NP_pb");
    auto h_acc_mid_pt_NP = (TH1D*) in_acc_pt->Get("mid");
    auto h_acc_mid_pt_NP_pp = (TH1D*) in_acc_pt->Get("mid_pp");
    auto h_acc_mid_pt_NP_pb = (TH1D*) in_acc_pt->Get("mid_pb");
    auto h_eff_mid_pt_NP = (TH1D*) in_eff_pt->Get("mid_NP");
    auto h_eff_mid_pt_NP_pp = (TH1D*) in_eff_pt->Get("mid_NP_pp");
    auto h_eff_mid_pt_NP_pb = (TH1D*) in_eff_pt->Get("mid_NP_pb");
    auto h_HF_mid_pt_NP = (TH1D*) in_HF_pt->Get("mid_NP");
    auto h_TNP_mid_pt_NP_pp = (TH1D*) in_TNPpp_pt->Get("mid_NP"); 
    auto h_TNP_mid_pt_NP_pb = (TH1D*) in_TNP_pt->Get("mid_NP"); 

    const int NBINS_mid_pt = 6;
    //double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 30, 50};
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    auto h_total_mid_pt_PR = new TH1D("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_PR_pp = new TH1D("mid_pt_PR_pp", "mid_PR_pp", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_PR_pb = new TH1D("mid_pt_PR_pb", "mid_PR_pb", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP = new TH1D("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP_pp = new TH1D("mid_pt_NP_pp", "mid_NP_pp", NBINS_mid_pt, edges_mid_pt);
    auto h_total_mid_pt_NP_pb = new TH1D("mid_pt_NP_pb", "mid_NP_pb", NBINS_mid_pt, edges_mid_pt);

    auto h_TNP_mid_pt_PR = new TH1D("TNP_mid_pt_PR", "TNP_mid_pt_PR", NBINS_mid_pt, edges_mid_pt);
    auto h_TNP_mid_pt_NP = new TH1D("TNP_mid_pt_NP", "TNP_mid_pt_NP", NBINS_mid_pt, edges_mid_pt);

    
    for (int hist_idx = 1; hist_idx < NBINS_mid_pt+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double TnP_mid_pt_PR = 0;
		//double acc_PR_pp = 0; double eff_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_mid_pt_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_mid_pt_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_mid_pt_PR_pp->GetBinContent(hist_idx);
        //double bFrac_PR_pp = h_bFrac_mid_pt_PR_pp->GetBinContent(hist_idx);
        double acc_PR_pp = h_acc_mid_pt_PR_pp->GetBinContent(hist_idx);
        double eff_PR_pp = h_eff_mid_pt_PR_pp->GetBinContent(hist_idx);
        double TNP_PR_pp = h_TNP_mid_pt_PR_pp->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                /*+ Power(bFrac_PR_pp,2)*/ + Power(acc_PR_pp,2) + Power(eff_PR_pp,2) + Power(TNP_PR_pp,2));
        
        h_total_mid_pt_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
		//double acc_PR_pb = 0; double eff_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_mid_pt_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_mid_pt_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_mid_pt_PR_pb->GetBinContent(hist_idx);
        //double bFrac_PR_pb = h_bFrac_mid_pt_PR_pb->GetBinContent(hist_idx);
        double acc_PR_pb = h_acc_mid_pt_PR_pb->GetBinContent(hist_idx);
        double eff_PR_pb = h_eff_mid_pt_PR_pb->GetBinContent(hist_idx);
        double HF_PR = h_HF_mid_pt_PR->GetBinContent(hist_idx);
        double TNP_PR_pb = h_TNP_mid_pt_PR_pb->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                /*+ Power(bFrac_PR_pb,2)*/ + Power(acc_PR_pb,2) + Power(eff_PR_pb,2) + Power(HF_PR,2) + Power(TNP_PR_pb,2));
        h_total_mid_pt_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
		//double acc_PR = 0; double eff_PR = 0;
        double sigPDF_PR = h_sigPDF_mid_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_pt_PR->GetBinContent(hist_idx);
        //double bFrac_PR = h_bFrac_mid_pt_PR->GetBinContent(hist_idx);
        double acc_PR = h_acc_mid_pt_PR->GetBinContent(hist_idx);
        double eff_PR = h_eff_mid_pt_PR->GetBinContent(hist_idx);
        double TNP_PR = Sqrt( Power(TNP_PR_pb,2) + Power(TNP_PR_pp,2) ); 
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                /*+ Power(bFrac_PR,2)*/ + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) + Power(TNP_PR,2));
        TnP_mid_pt_PR = Sqrt( Power(TNP_PR_pp,2) + Power(TNP_PR_pb,2) );
        h_TNP_mid_pt_PR->SetBinContent(hist_idx, TnP_mid_pt_PR);
        h_total_mid_pt_PR->SetBinContent(hist_idx, tot_syst_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double TnP_mid_pt_NP = 0;
		//double acc_NP_pp = 0; double eff_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_mid_pt_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_mid_pt_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_mid_pt_NP_pp->GetBinContent(hist_idx);
        //double bFrac_NP_pp = h_bFrac_mid_pt_NP_pp->GetBinContent(hist_idx);
        double acc_NP_pp = h_acc_mid_pt_NP_pp->GetBinContent(hist_idx);
        double eff_NP_pp = h_eff_mid_pt_NP_pp->GetBinContent(hist_idx);
        double TNP_NP_pp = h_TNP_mid_pt_NP_pp->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                /*+ Power(bFrac_NP_pp,2)*/ + Power(acc_NP_pp,2) + Power(eff_NP_pp,2) + Power(TNP_NP_pp,2));
        h_total_mid_pt_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
		//double acc_NP_pb = 0; double eff_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_mid_pt_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_mid_pt_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_mid_pt_NP_pb->GetBinContent(hist_idx);
        //double bFrac_NP_pb = h_bFrac_mid_pt_NP_pb->GetBinContent(hist_idx);
        double acc_NP_pb = h_acc_mid_pt_NP_pb->GetBinContent(hist_idx);
        double eff_NP_pb = h_eff_mid_pt_NP_pb->GetBinContent(hist_idx);
        double HF_NP = h_HF_mid_pt_NP->GetBinContent(hist_idx);
        double TNP_NP_pb = h_TNP_mid_pt_NP_pb->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                /*+ Power(bFrac_NP_pb,2)*/ + Power(acc_NP_pb,2) + Power(eff_NP_pb,2) + Power(HF_NP,2) + Power(TNP_NP_pb,2));
        TnP_mid_pt_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) );
        h_total_mid_pt_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        h_TNP_mid_pt_NP->SetBinContent(hist_idx, TnP_mid_pt_NP);
        //cout << h_total_mid_pt_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
		//double acc_NP = 0; double eff_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_pt_NP->GetBinContent(hist_idx);
        //double bFrac_NP = h_bFrac_mid_pt_NP->GetBinContent(hist_idx);
        double acc_NP = h_acc_mid_pt_NP->GetBinContent(hist_idx);
        double eff_NP = h_eff_mid_pt_NP->GetBinContent(hist_idx);
        double TNP_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) ); 
		cout << "TNP NP Tot : " << TNP_NP << endl;
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                /*+ Power(bFrac_NP,2)*/ + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) + Power(TNP_NP,2));
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
    //auto h_bFrac_fwd_pt_PR = (TH1D*) in_bFrac_pt->Get("fwd_PR");
    //auto h_bFrac_fwd_pt_PR_pp = (TH1D*) in_bFrac_pt->Get("fwd_PR_pp");
    //auto h_bFrac_fwd_pt_PR_pb = (TH1D*) in_bFrac_pt->Get("fwd_PR_pb");
    auto h_acc_fwd_pt_PR = (TH1D*) in_acc_pt->Get("mid");
    auto h_acc_fwd_pt_PR_pp = (TH1D*) in_acc_pt->Get("fwd_pp");
    auto h_acc_fwd_pt_PR_pb = (TH1D*) in_acc_pt->Get("fwd_pb");
    auto h_eff_fwd_pt_PR = (TH1D*) in_eff_pt->Get("fwd_PR");
    auto h_eff_fwd_pt_PR_pp = (TH1D*) in_eff_pt->Get("fwd_PR_pp");
    auto h_eff_fwd_pt_PR_pb = (TH1D*) in_eff_pt->Get("fwd_PR_pb");
    auto h_HF_fwd_pt_PR = (TH1D*) in_HF_pt->Get("fwd_PR");
    auto h_TNP_fwd_pt_PR_pp = (TH1D*) in_TNPpp_pt->Get("fwd_PR");
    auto h_TNP_fwd_pt_PR_pb = (TH1D*) in_TNP_pt->Get("fwd_PR");

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
    //auto h_bFrac_fwd_pt_NP = (TH1D*) in_bFrac_pt->Get("fwd_NP");
    //auto h_bFrac_fwd_pt_NP_pp = (TH1D*) in_bFrac_pt->Get("fwd_NP_pp");
    //auto h_bFrac_fwd_pt_NP_pb = (TH1D*) in_bFrac_pt->Get("fwd_NP_pb");
    auto h_acc_fwd_pt_NP = (TH1D*) in_acc_pt->Get("mid");
    auto h_acc_fwd_pt_NP_pp = (TH1D*) in_acc_pt->Get("fwd_pp");
    auto h_acc_fwd_pt_NP_pb = (TH1D*) in_acc_pt->Get("fwd_pb");
    auto h_eff_fwd_pt_NP = (TH1D*) in_eff_pt->Get("fwd_NP");
    auto h_eff_fwd_pt_NP_pp = (TH1D*) in_eff_pt->Get("fwd_NP_pp");
    auto h_eff_fwd_pt_NP_pb = (TH1D*) in_eff_pt->Get("fwd_NP_pb");
    auto h_HF_fwd_pt_NP = (TH1D*) in_HF_pt->Get("fwd_NP");
    auto h_TNP_fwd_pt_NP_pp = (TH1D*) in_TNPpp_pt->Get("fwd_NP");
    auto h_TNP_fwd_pt_NP_pb = (TH1D*) in_TNP_pt->Get("fwd_NP");

    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 40};
    auto h_total_fwd_pt_PR = new TH1D("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_PR_pp = new TH1D("fwd_pt_PR_pp", "fwd_PR_pp", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_PR_pb = new TH1D("fwd_pt_PR_pb", "fwd_PR_pb", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP = new TH1D("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP_pp = new TH1D("fwd_pt_NP_pp", "fwd_NP_pp", NBINS_fwd_pt, edges_fwd_pt);
    auto h_total_fwd_pt_NP_pb = new TH1D("fwd_pt_NP_pb", "fwd_NP_pb", NBINS_fwd_pt, edges_fwd_pt);

    auto h_TNP_fwd_pt_PR = new TH1D("TNP_fwd_pt_PR", "TNP_fwd_pt_PR", NBINS_mid_pt+1, edges_mid_pt);
    auto h_TNP_fwd_pt_NP = new TH1D("TNP_fwd_pt_NP", "TNP_fwd_pt_NP", NBINS_mid_pt+1, edges_mid_pt);
    
    for (int hist_idx = 1; hist_idx < NBINS_fwd_pt+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double TnP_fwd_pt_PR = 0;
		//double acc_PR_pp = 0; double eff_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_fwd_pt_PR_pp->GetBinContent(hist_idx);
        //double bFrac_PR_pp = h_bFrac_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double acc_PR_pp = h_acc_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double eff_PR_pp = h_eff_fwd_pt_PR_pp->GetBinContent(hist_idx);
        double TNP_PR_pp = h_TNP_fwd_pt_PR_pp->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                /*+ Power(bFrac_PR_pp,2)*/ + Power(acc_PR_pp,2) + Power(eff_PR_pp,2) + Power(TNP_PR_pp,2));
        h_total_fwd_pt_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
		//double acc_PR_pb = 0; double eff_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_fwd_pt_PR_pb->GetBinContent(hist_idx);
        //double bFrac_PR_pb = h_bFrac_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double acc_PR_pb = h_acc_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double eff_PR_pb = h_eff_fwd_pt_PR_pb->GetBinContent(hist_idx);
        double HF_PR = h_HF_fwd_pt_PR->GetBinContent(hist_idx);
        double TNP_PR_pb = h_TNP_fwd_pt_PR_pb->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                /*+ Power(bFrac_PR_pb,2)*/ + Power(acc_PR_pb,2) + Power(eff_PR_pb,2) + Power(HF_PR,2) + Power(TNP_PR_pb,2));
        h_total_fwd_pt_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);
        cout << "TnP PR Pb : " << TNP_PR_pb << "\tTotal syst PR Pb : " << tot_syst_PR_pb << endl;

        double tot_syst_PR = 0;
		//double acc_PR = 0; double eff_PR = 0;
        double sigPDF_PR = h_sigPDF_fwd_pt_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_pt_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_pt_PR->GetBinContent(hist_idx);
        //double bFrac_PR = h_bFrac_fwd_pt_PR->GetBinContent(hist_idx);
        double acc_PR = h_acc_fwd_pt_PR->GetBinContent(hist_idx);
        double eff_PR = h_eff_fwd_pt_PR->GetBinContent(hist_idx);
        double TNP_PR = Sqrt( Power(TNP_PR_pb,2) + Power(TNP_PR_pp,2) ); 
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                /*+ Power(bFrac_PR,2)*/ + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) + Power(TNP_PR,2));
        TnP_fwd_pt_PR = Sqrt( Power(TNP_PR_pp,2) + Power(TNP_PR_pb,2) );
        h_total_fwd_pt_PR->SetBinContent(hist_idx, tot_syst_PR);
        h_TNP_fwd_pt_PR->SetBinContent(hist_idx, TnP_fwd_pt_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double TnP_fwd_pt_NP = 0;
		//double acc_NP_pp = 0; double eff_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_fwd_pt_NP_pp->GetBinContent(hist_idx);
        //double bFrac_NP_pp = h_bFrac_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double acc_NP_pp = h_acc_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double eff_NP_pp = h_eff_fwd_pt_NP_pp->GetBinContent(hist_idx);
        double TNP_NP_pp = h_TNP_fwd_pt_NP_pp->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                /*+ Power(bFrac_NP_pp,2)*/ + Power(acc_NP_pp,2) + Power(eff_NP_pp,2) + Power(TNP_NP_pp,2));
        h_total_fwd_pt_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
		//double acc_NP_pb = 0; double eff_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_fwd_pt_NP_pb->GetBinContent(hist_idx);
        //double bFrac_NP_pb = h_bFrac_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double acc_NP_pb = h_acc_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double eff_NP_pb = h_eff_fwd_pt_NP_pb->GetBinContent(hist_idx);
        double HF_NP = h_HF_fwd_pt_NP->GetBinContent(hist_idx);
        double TNP_NP_pb = h_TNP_fwd_pt_NP_pb->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                /*+ Power(bFrac_NP_pb,2)*/ + Power(acc_NP_pb,2) + Power(eff_NP_pb,2) + Power(HF_NP,2) + Power(TNP_NP_pb,2));
        h_total_fwd_pt_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        TnP_fwd_pt_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) );
        //cout << h_total_fwd_pt_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
		//double acc_NP = 0; double eff_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_pt_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_pt_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_pt_NP->GetBinContent(hist_idx);
        //double bFrac_NP = h_bFrac_fwd_pt_NP->GetBinContent(hist_idx);
        double acc_NP = h_acc_fwd_pt_NP->GetBinContent(hist_idx);
        double eff_NP = h_eff_fwd_pt_NP->GetBinContent(hist_idx);
        double TNP_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) ); 
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                /*+ Power(bFrac_NP,2)*/ + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) + Power(TNP_NP,2));
        h_total_fwd_pt_NP->SetBinContent(hist_idx, tot_syst_NP);
        h_TNP_fwd_pt_NP->SetBinContent(hist_idx, TnP_fwd_pt_NP);
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
    //auto h_bFrac_mid_cent_PR = (TH1D*) in_bFrac_cent->Get("mid_PR");
    //auto h_bFrac_mid_cent_PR_pp = (TH1D*) in_bFrac_cent->Get("mid_PR_pp");
    //auto h_bFrac_mid_cent_PR_pb = (TH1D*) in_bFrac_cent->Get("mid_PR_pb");
    auto h_acc_mid_cent_PR = (TH1D*) in_acc_cent->Get("mid");
    auto h_acc_mid_cent_PR_pp = (TH1D*) in_acc_cent->Get("mid_pp");
    auto h_acc_mid_cent_PR_pb = (TH1D*) in_acc_cent->Get("mid_pb");
    auto h_eff_mid_cent_PR = (TH1D*) in_eff_cent->Get("mid_PR");
    auto h_eff_mid_cent_PR_pp = (TH1D*) in_eff_cent->Get("mid_PR_pp");
    auto h_eff_mid_cent_PR_pb = (TH1D*) in_eff_cent->Get("mid_PR_pb");
    auto h_HF_mid_cent_PR = (TH1D*) in_HF_cent->Get("mid_PR");
    //auto h_TNP_mid_cent_PR_pb = (TH1D*) in_TNP_cent->Get("mid_PR_pb");
    auto h_TNP_mid_cent_PR_pp = (TH1D*) in_TNPpp_cent->Get("mid_PR");
    auto h_TNP_mid_cent_PR_pb = (TH1D*) in_TNP_cent->Get("mid_PR");

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
    //auto h_bFrac_mid_cent_NP = (TH1D*) in_bFrac_cent->Get("mid_NP");
    //auto h_bFrac_mid_cent_NP_pp = (TH1D*) in_bFrac_cent->Get("mid_NP_pp");
    //auto h_bFrac_mid_cent_NP_pb = (TH1D*) in_bFrac_cent->Get("mid_NP_pb");
    auto h_acc_mid_cent_NP = (TH1D*) in_acc_cent->Get("mid");
    auto h_acc_mid_cent_NP_pp = (TH1D*) in_acc_cent->Get("mid_pp");
    auto h_acc_mid_cent_NP_pb = (TH1D*) in_acc_cent->Get("mid_pb");
    auto h_eff_mid_cent_NP = (TH1D*) in_eff_cent->Get("mid_NP");
    auto h_eff_mid_cent_NP_pp = (TH1D*) in_eff_cent->Get("mid_NP_pp");
    auto h_eff_mid_cent_NP_pb = (TH1D*) in_eff_cent->Get("mid_NP_pb");
    auto h_HF_mid_cent_NP = (TH1D*) in_HF_cent->Get("mid_NP");
    //auto h_TNP_mid_cent_NP_pb = (TH1D*) in_TNP_cent->Get("mid_NP_pb");
    auto h_TNP_mid_cent_NP_pp = (TH1D*) in_TNPpp_cent->Get("mid_NP");
    auto h_TNP_mid_cent_NP_pb = (TH1D*) in_TNP_cent->Get("mid_NP");

    

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0,10,20,30,40,50,90};
    auto h_total_mid_cent_PR = new TH1D("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_PR_pp = new TH1D("mid_cent_PR_pp", "mid_PR_pp", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_PR_pb = new TH1D("mid_cent_PR_pb", "mid_PR_pb", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP = new TH1D("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP_pp = new TH1D("mid_cent_NP_pp", "mid_NP_pp", NBINS_mid_cent, edges_mid_cent);
    auto h_total_mid_cent_NP_pb = new TH1D("mid_cent_NP_pb", "mid_NP_pb", NBINS_mid_cent, edges_mid_cent);
    auto h_TNP_mid_cent_PR = new TH1D("TNP_mid_cent_PR", "TNP_mid_cent_PR", NBINS_mid_cent+1, edges_mid_cent);
    auto h_TNP_mid_cent_NP = new TH1D("TNP_mid_cent_NP", "TNP_mid_cent_NP", NBINS_mid_cent+1, edges_mid_cent);

    for (int hist_idx = 1; hist_idx < NBINS_mid_cent+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double TnP_mid_cent_PR = 0;
        
		//double acc_PR_pp = 0; double eff_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_mid_cent_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_mid_cent_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_mid_cent_PR_pp->GetBinContent(hist_idx);
        
        //double bFrac_PR_pp = h_bFrac_mid_cent_PR_pp->GetBinContent(hist_idx);
        double acc_PR_pp = h_acc_mid_cent_PR_pp->GetBinContent(hist_idx);
        double eff_PR_pp = h_eff_mid_cent_PR_pp->GetBinContent(hist_idx);
        double TNP_PR_pp = h_TNP_mid_cent_PR_pp->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                /*+ Power(bFrac_PR_pp,2)*/ + Power(acc_PR_pp,2) + Power(eff_PR_pp,2) + Power(TNP_PR_pp,2));
        h_total_mid_cent_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
		//double acc_PR_pb = 0; double eff_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_mid_cent_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_mid_cent_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_mid_cent_PR_pb->GetBinContent(hist_idx);
        //double bFrac_PR_pb = h_bFrac_mid_cent_PR_pb->GetBinContent(hist_idx);
        double acc_PR_pb = h_acc_mid_cent_PR_pb->GetBinContent(hist_idx);
        double eff_PR_pb = h_eff_mid_cent_PR_pb->GetBinContent(hist_idx);
        double HF_PR = h_HF_mid_cent_PR->GetBinContent(hist_idx);
        double TNP_PR_pb = h_TNP_mid_cent_PR_pb->GetBinContent(hist_idx+1);
		cout << "HF PR : " << HF_PR << endl;
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                /*+ Power(bFrac_PR_pb,2)*/ + Power(acc_PR_pb,2) + Power(eff_PR_pb,2) + Power(HF_PR,2) + Power(TNP_PR_pb,2));
        h_total_mid_cent_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
		//double acc_PR = 0; double eff_PR = 0;
        double sigPDF_PR = h_sigPDF_mid_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_mid_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_mid_cent_PR->GetBinContent(hist_idx);
        //double bFrac_PR = h_bFrac_mid_cent_PR->GetBinContent(hist_idx);
        double acc_PR = h_acc_mid_cent_PR->GetBinContent(hist_idx);
        double eff_PR = h_eff_mid_cent_PR->GetBinContent(hist_idx);
        double TNP_PR = Sqrt( Power(TNP_PR_pb,2) + Power(TNP_PR_pp,2) ); 
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                /*+ Power(bFrac_PR,2)*/ + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) + Power(TNP_PR,2));
        TnP_mid_cent_PR = Sqrt( Power(TNP_PR_pp,2) + Power(TNP_PR_pb,2) );
        h_total_mid_cent_PR->SetBinContent(hist_idx, tot_syst_PR);
        h_TNP_mid_cent_PR->SetBinContent(hist_idx, TnP_mid_cent_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double TnP_mid_cent_NP = 0;
		//double acc_NP_pp = 0; double eff_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_mid_cent_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_mid_cent_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_mid_cent_NP_pp->GetBinContent(hist_idx);
        //double bFrac_NP_pp = h_bFrac_mid_cent_NP_pp->GetBinContent(hist_idx);
        double acc_NP_pp = h_acc_mid_cent_NP_pp->GetBinContent(hist_idx);
        double eff_NP_pp = h_eff_mid_cent_NP_pp->GetBinContent(hist_idx);
        double TNP_NP_pp = h_TNP_mid_cent_NP_pp->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                /*+ Power(bFrac_NP_pp,2)*/ + Power(acc_NP_pp,2) + Power(eff_NP_pp,2) + Power(TNP_NP_pp,2));
        h_total_mid_cent_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
		//double acc_NP_pb = 0; double eff_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_mid_cent_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_mid_cent_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_mid_cent_NP_pb->GetBinContent(hist_idx);
        //double bFrac_NP_pb = h_bFrac_mid_cent_NP_pb->GetBinContent(hist_idx);
        double acc_NP_pb = h_acc_mid_cent_NP_pb->GetBinContent(hist_idx);
        double eff_NP_pb = h_eff_mid_cent_NP_pb->GetBinContent(hist_idx);
        double HF_NP = h_HF_mid_cent_NP->GetBinContent(hist_idx);
        double TNP_NP_pb = h_TNP_mid_cent_NP_pb->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                /*+ Power(bFrac_NP_pb,2)*/ + Power(acc_NP_pb,2) + Power(eff_NP_pb,2) + Power(HF_NP,2) + Power(TNP_NP_pb,2));
        h_total_mid_cent_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        //cout << h_total_mid_cent_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
		//double acc_NP = 0; double eff_NP = 0;
        double sigPDF_NP = h_sigPDF_mid_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_mid_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_mid_cent_NP->GetBinContent(hist_idx);
        //double bFrac_NP = h_bFrac_mid_cent_NP->GetBinContent(hist_idx);
        double acc_NP = h_acc_mid_cent_NP->GetBinContent(hist_idx);
        double eff_NP = h_eff_mid_cent_NP->GetBinContent(hist_idx);
        double TNP_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) ); 
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                /*+ Power(bFrac_NP,2)*/ + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) + Power(TNP_NP,2));
        TnP_mid_cent_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) );
        h_total_mid_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
        h_TNP_mid_cent_NP->SetBinContent(hist_idx, TnP_mid_cent_NP);
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
    //auto h_bFrac_fwd_cent_PR = (TH1D*) in_bFrac_cent->Get("fwd_PR");
    //auto h_bFrac_fwd_cent_PR_pp = (TH1D*) in_bFrac_cent->Get("fwd_PR_pp");
    //auto h_bFrac_fwd_cent_PR_pb = (TH1D*) in_bFrac_cent->Get("fwd_PR_pb");
    auto h_acc_fwd_cent_PR = (TH1D*) in_acc_cent->Get("fwd");
    auto h_acc_fwd_cent_PR_pp = (TH1D*) in_acc_cent->Get("fwd_pp");
    auto h_acc_fwd_cent_PR_pb = (TH1D*) in_acc_cent->Get("fwd_pb");
    auto h_eff_fwd_cent_PR = (TH1D*) in_eff_cent->Get("fwd_PR");
    auto h_eff_fwd_cent_PR_pp = (TH1D*) in_eff_cent->Get("fwd_PR_pp");
    auto h_eff_fwd_cent_PR_pb = (TH1D*) in_eff_cent->Get("fwd_PR_pb");
    auto h_HF_fwd_cent_PR = (TH1D*) in_HF_cent->Get("fwd_PR");
    auto h_TNP_fwd_cent_PR_pp = (TH1D*) in_TNPpp_cent->Get("fwd_PR");
    auto h_TNP_fwd_cent_PR_pb = (TH1D*) in_TNP_cent->Get("fwd_PR");

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
    //auto h_bFrac_fwd_cent_NP = (TH1D*) in_bFrac_cent->Get("fwd_NP");
    //auto h_bFrac_fwd_cent_NP_pp = (TH1D*) in_bFrac_cent->Get("fwd_NP_pp");
    //auto h_bFrac_fwd_cent_NP_pb = (TH1D*) in_bFrac_cent->Get("fwd_NP_pb");
    auto h_acc_fwd_cent_NP = (TH1D*) in_acc_cent->Get("fwd");
    auto h_acc_fwd_cent_NP_pp = (TH1D*) in_acc_cent->Get("fwd_pp");
    auto h_acc_fwd_cent_NP_pb = (TH1D*) in_acc_cent->Get("fwd_pb");
    auto h_eff_fwd_cent_NP = (TH1D*) in_eff_cent->Get("fwd_NP");
    auto h_eff_fwd_cent_NP_pp = (TH1D*) in_eff_cent->Get("fwd_NP_pp");
    auto h_eff_fwd_cent_NP_pb = (TH1D*) in_eff_cent->Get("fwd_NP_pb");
    auto h_HF_fwd_cent_NP = (TH1D*) in_HF_cent->Get("fwd_NP");
    auto h_TNP_fwd_cent_NP_pp = (TH1D*) in_TNPpp_cent->Get("fwd_NP");
    auto h_TNP_fwd_cent_NP_pb = (TH1D*) in_TNP_cent->Get("fwd_NP");


    const int NBINS_fwd_cent = 4;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 30, 50, 90};
    auto h_total_fwd_cent_PR = new TH1D("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_PR_pp = new TH1D("fwd_cent_PR_pp", "fwd_PR_pp", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_PR_pb = new TH1D("fwd_cent_PR_pb", "fwd_PR_pb", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP = new TH1D("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP_pp = new TH1D("fwd_cent_NP_pp", "fwd_NP_pp", NBINS_fwd_cent, edges_fwd_cent);
    auto h_total_fwd_cent_NP_pb = new TH1D("fwd_cent_NP_pb", "fwd_NP_pb", NBINS_fwd_cent, edges_fwd_cent);
    auto h_TNP_fwd_cent_PR = new TH1D("TNP_fwd_cent_PR", "TNP_fwd_cent_PR", NBINS_fwd_cent+1, edges_fwd_cent);
    auto h_TNP_fwd_cent_NP = new TH1D("TNP_fwd_cent_NP", "TNP_fwd_cent_NP", NBINS_fwd_cent+1, edges_fwd_cent);
    
    for (int hist_idx = 1; hist_idx < NBINS_fwd_cent+1; hist_idx++) {
        double tot_syst_PR_pp = 0;
        double TnP_fwd_cent_PR = 0;
		//double acc_PR_pp = 0; double eff_PR_pp = 0;
        double sigPDF_PR_pp = h_sigPDF_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double sigPAR_PR_pp = h_sigPAR_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double bkgPDF_PR_pp = h_bkgPDF_fwd_cent_PR_pp->GetBinContent(hist_idx);
        //double bFrac_PR_pp = h_bFrac_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double acc_PR_pp = h_acc_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double eff_PR_pp = h_eff_fwd_cent_PR_pp->GetBinContent(hist_idx);
        double TNP_PR_pp = h_TNP_fwd_cent_PR_pp->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pp = Sqrt( Power(sigPDF_PR_pp,2) + Power(sigPAR_PR_pp,2) + Power(bkgPDF_PR_pp,2)
                /*+ Power(bFrac_PR_pp,2)*/ + Power(acc_PR_pp,2) + Power(eff_PR_pp,2) + Power(TNP_PR_pp,2));
        h_total_fwd_cent_PR_pp->SetBinContent(hist_idx, tot_syst_PR_pp);

        double tot_syst_PR_pb = 0;
		//double acc_PR_pb = 0; double eff_PR_pb = 0;
        double sigPDF_PR_pb = h_sigPDF_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double sigPAR_PR_pb = h_sigPAR_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double bkgPDF_PR_pb = h_bkgPDF_fwd_cent_PR_pb->GetBinContent(hist_idx);
        //double bFrac_PR_pb = h_bFrac_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double acc_PR_pb = h_acc_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double eff_PR_pb = h_eff_fwd_cent_PR_pb->GetBinContent(hist_idx);
        double HF_PR = h_HF_fwd_cent_PR->GetBinContent(hist_idx);
        double TNP_PR_pb = h_TNP_fwd_cent_PR_pb->GetBinContent(hist_idx+1);
        
        tot_syst_PR_pb = Sqrt( Power(sigPDF_PR_pb,2) + Power(sigPAR_PR_pb,2) + Power(bkgPDF_PR_pb,2)
                /*+ Power(bFrac_PR_pb,2)*/ + Power(acc_PR_pb,2) + Power(eff_PR_pb,2) + Power(HF_PR,2) + Power(TNP_PR_pb,2));
        h_total_fwd_cent_PR_pb->SetBinContent(hist_idx, tot_syst_PR_pb);

        double tot_syst_PR = 0;
		//double acc_PR = 0; double eff_PR = 0;
        double sigPDF_PR = h_sigPDF_fwd_cent_PR->GetBinContent(hist_idx);
        double sigPAR_PR = h_sigPAR_fwd_cent_PR->GetBinContent(hist_idx);
        double bkgPDF_PR = h_bkgPDF_fwd_cent_PR->GetBinContent(hist_idx);
        //double bFrac_PR = h_bFrac_fwd_cent_PR->GetBinContent(hist_idx);
        double acc_PR = h_acc_fwd_cent_PR->GetBinContent(hist_idx);
        double eff_PR = h_eff_fwd_cent_PR->GetBinContent(hist_idx);
        double TNP_PR = Sqrt( Power(TNP_PR_pb,2) + Power(TNP_PR_pp,2) ); 
        
        tot_syst_PR = Sqrt( Power(sigPDF_PR,2) + Power(sigPAR_PR,2) + Power(bkgPDF_PR,2)
                /*+ Power(bFrac_PR,2)*/ + Power(acc_PR,2) + Power(eff_PR,2) + Power(HF_PR,2) + Power(TNP_PR,2));

        TnP_fwd_cent_PR = Sqrt( Power(TNP_PR_pp,2) + Power(TNP_PR_pb,2) );
        h_total_fwd_cent_PR->SetBinContent(hist_idx, tot_syst_PR);
        h_TNP_fwd_cent_PR->SetBinContent(hist_idx, TnP_fwd_cent_PR);

        // NP
        double tot_syst_NP_pp = 0;
        double TnP_fwd_cent_NP = 0;
		//double acc_NP_pp = 0; double eff_NP_pp = 0;
        double sigPDF_NP_pp = h_sigPDF_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double sigPAR_NP_pp = h_sigPAR_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double bkgPDF_NP_pp = h_bkgPDF_fwd_cent_NP_pp->GetBinContent(hist_idx);
        //double bFrac_NP_pp = h_bFrac_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double acc_NP_pp = h_acc_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double eff_NP_pp = h_eff_fwd_cent_NP_pp->GetBinContent(hist_idx);
        double TNP_NP_pp = h_TNP_fwd_cent_NP_pp->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pp = Sqrt( Power(sigPDF_NP_pp,2) + Power(sigPAR_NP_pp,2) + Power(bkgPDF_NP_pp,2)
                /*+ Power(bFrac_NP_pp,2)*/ + Power(acc_NP_pp,2) + Power(eff_NP_pp,2) + Power(TNP_NP_pp,2));
        h_total_fwd_cent_NP_pp->SetBinContent(hist_idx, tot_syst_NP_pp);

        double tot_syst_NP_pb = 0;
		//double acc_NP_pb = 0; double eff_NP_pb = 0;
        double sigPDF_NP_pb = h_sigPDF_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double sigPAR_NP_pb = h_sigPAR_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double bkgPDF_NP_pb = h_bkgPDF_fwd_cent_NP_pb->GetBinContent(hist_idx);
        //double bFrac_NP_pb = h_bFrac_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double acc_NP_pb = h_acc_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double eff_NP_pb = h_eff_fwd_cent_NP_pb->GetBinContent(hist_idx);
        double HF_NP = h_HF_fwd_cent_NP->GetBinContent(hist_idx);
        double TNP_NP_pb = h_TNP_fwd_cent_NP_pb->GetBinContent(hist_idx+1);
        
        tot_syst_NP_pb = Sqrt( Power(sigPDF_NP_pb,2) + Power(sigPAR_NP_pb,2) + Power(bkgPDF_NP_pb,2)
                /*+ Power(bFrac_NP_pb,2)*/ + Power(acc_NP_pb,2) + Power(eff_NP_pb,2) + Power(HF_NP,2) + Power(TNP_NP_pb,2));
        h_total_fwd_cent_NP_pb->SetBinContent(hist_idx, tot_syst_NP_pb);
        //cout << h_total_fwd_cent_NP->GetBinContent(hist_idx) << endl;
		//
        double tot_syst_NP = 0;
		//double acc_NP = 0; double eff_NP = 0;
        double sigPDF_NP = h_sigPDF_fwd_cent_NP->GetBinContent(hist_idx);
        double sigPAR_NP = h_sigPAR_fwd_cent_NP->GetBinContent(hist_idx);
        double bkgPDF_NP = h_bkgPDF_fwd_cent_NP->GetBinContent(hist_idx);
        //double bFrac_NP = h_bFrac_fwd_cent_NP->GetBinContent(hist_idx);
        double acc_NP = h_acc_fwd_cent_NP->GetBinContent(hist_idx);
        double eff_NP = h_eff_fwd_cent_NP->GetBinContent(hist_idx);
        double TNP_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) ); 
        
        tot_syst_NP = Sqrt( Power(sigPDF_NP,2) + Power(sigPAR_NP,2) + Power(bkgPDF_NP,2)
                /*+ Power(bFrac_NP,2)*/ + Power(acc_NP,2) + Power(eff_NP,2) + Power(HF_NP,2) + Power(TNP_NP,2));
        TnP_fwd_cent_NP = Sqrt( Power(TNP_NP_pp,2) + Power(TNP_NP_pb,2) );
        h_total_fwd_cent_NP->SetBinContent(hist_idx, tot_syst_NP);
        h_TNP_fwd_cent_NP->SetBinContent(hist_idx, TnP_fwd_cent_NP);
    }
    out_file.cd();
    out_file.Write();
    out_file.Close();
}
