#include "TH1.h"
#include "TFile.h"
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst);
double compute_uncertainty_pp(double pp_nomi, double pp_syst);
double compute_uncertainty_pb(double pb_nomi, double pb_syst);

void syst_acc()
{
    gSystem->mkdir("syst_roots");
    string syst_type = "acc";

    // Path label
    /*
    TString path = "./Efficiency/20230607_EffNom_";
    TString fname_pp = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW%d_tnp1_20230602.root";
    TString fname_pb = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW%d_tnp1_20230602.root";
    TString pp_nomi_PR = "pp/roots/"+fname_pp;
    TString pp_nomi_NP = "pp/roots/"+fname_pp;
    TString pp_syst_PR = "pp_Off/roots/"+fname_pp;
    TString pp_syst_NP = "pp_Off/roots/"+fname_pp;
    TString pb_nomi_PR = "pbpb_On/roots/"+fname_pb;
    TString pb_nomi_NP = "pbpb_On_Npr/roots/"+fname_pb;
    TString pb_syst_PR = "pbpb_Off/roots/"+fname_pb;
    TString pb_syst_NP = "pbpb_Off_Npr/roots/"+fname_pb;
*/


    
    TString path = "../../Eff_Acc/roots/";
	//TString pb_nomi = "acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20230728.root";
	TString pb_nomi = "acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20240326_ppWfunc.root";
	//TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root";
	//TString pp_nomi = "acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240314.root";
	TString pp_nomi = "acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240319.root";
    //TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20231019.root";
    //TString pp_nomi_NP = "mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW1_tnp1_20230728.root";
    //TString pb_nomi_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20231019.root";
    //TString pb_nomi_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2s_PtW1_tnp1_new_20231019.root";

    TString pp_syst = "acceptance_Prompt_psi2s_GenOnly_wgt0_pp_SysUp0_20240319.root";
    TString pb_syst = "acceptance_Prompt_psi2s_GenOnly_wgt0_PbPb_SysUp0_20240314.root";
    

    TString h_pp_mid = "hAccPt_2021_midy";
    TString h_pb_mid = "hAccPt_2021_midy";
    

    // Define pointers for input files
    TFile *pp_nomi_input = new TFile(path+pp_nomi.Data());
    TFile *pb_nomi_input = new TFile(path+pb_nomi.Data());
    TFile *pp_syst_input = new TFile(path+pp_syst.Data());
    TFile *pb_syst_input = new TFile(path+pb_syst.Data());
    
    TH1D *h_pb_nomi = (TH1D *)pb_nomi_input->Get(h_pb_mid.Data());
    TH1D *h_pb_syst = (TH1D *)pb_syst_input->Get(h_pb_mid.Data());
    TH1D *h_pp_nomi = (TH1D *)pp_nomi_input->Get(h_pp_mid.Data());
    TH1D *h_pp_syst = (TH1D *)pp_syst_input->Get(h_pp_mid.Data());


    // Start loop
        // loop1 - mid_pt
        // loop2 - fwd_pt, results of loop 1 and 2 are saved into syst_pt_[something].root
        // loop3 - mid_cent
        // loop4 - fwd_cent, results of loop 3 and 4 are saved into syst_cent_[something].root
    
    // Start loop1 - mid_pt
    string out_name = "./syst_roots/syst_pt_" + syst_type + ".root";
    TFile out_pt(out_name.c_str(), "recreate");

    const int NBINS_mid_pt = 6;
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 50};
    TH1D mid_pt("mid_pt", "mid", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_pp("mid_pt_pp", "mid_pp", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_pb("mid_pt_pb", "mid_pb", NBINS_mid_pt, edges_mid_pt);

    for (int i = 0; i < NBINS_mid_pt; i++) {
        // Get number of PR and NP Jpsi
        // pp nominal
        double n_pp_nomi = h_pp_nomi->GetBinContent(i+1);
        
        // pp syst
        double n_pp_syst = h_pp_syst->GetBinContent(i+1);
        
        // Pb nominal
        double n_Pb_nomi = h_pb_nomi->GetBinContent(i+1);
        
        // Pb syst
        double n_Pb_syst = h_pb_syst->GetBinContent(i+1);


        // Compute uncertainty
        // Proppt
        double uncert = compute_uncertainty(n_pp_nomi, n_pp_syst, n_Pb_nomi, n_Pb_syst);
        double uncert_pp = compute_uncertainty_pp(n_pp_nomi, n_pp_syst);
        double uncert_pb = compute_uncertainty_pb(n_Pb_nomi, n_Pb_syst);
        
        
        // Non-prompt
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_pt.SetBinContent(i+1, uncert); // i starts from 0, hist elements starts from 1
        mid_pt_pp.SetBinContent(i+1, uncert_pp);
        mid_pt_pb.SetBinContent(i+1, uncert_pb);
    }

    // loop2 - fwd_pt
    TString h_pp_fwd = "hAccPt_2021_Fory";
    TString h_pb_fwd = "hAccPt_2021_Fory";

    h_pb_nomi = (TH1D *)pb_nomi_input->Get(h_pb_fwd.Data());
    h_pb_syst = (TH1D *)pb_syst_input->Get(h_pb_fwd.Data());
    h_pp_nomi = (TH1D *)pp_nomi_input->Get(h_pp_fwd.Data());
    h_pp_syst = (TH1D *)pp_syst_input->Get(h_pp_fwd.Data());


    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 50};
    TH1D fwd_pt("fwd_pt", "fwd", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_pp("fwd_pt_pp", "fwd_pp", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_pb("fwd_pt_pb", "fwd_pb", NBINS_fwd_pt, edges_fwd_pt);
    for (int i = 0; i < pp_fwd_pt.size(); i++) {
         // Open input files
        // Get number of PR and NP Jpsi
        // pp nominal
        double n_pp_nomi = h_pp_nomi->GetBinContent(i+1);
        
        // pp syst
        double n_pp_syst = h_pp_syst->GetBinContent(i+1);
        
        // Pb nominal
        double n_Pb_nomi = h_pb_nomi->GetBinContent(i+1);

        // Pb syst
        double n_Pb_syst = h_pb_syst->GetBinContent(i+1);


        // Compute uncertainty
        // Proppt
        double uncert = compute_uncertainty(n_pp_nomi, n_pp_syst, n_Pb_nomi, n_Pb_syst);
        double uncert_pp = compute_uncertainty_pp(n_pp_nomi, n_pp_syst);
        double uncert_pb = compute_uncertainty_pb(n_Pb_nomi, n_Pb_syst);
        
        
        // Non-prompt
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        fwd_pt.SetBinContent(i+1, uncert); // i starts from 0, hist elements starts from 1
        fwd_pt_pp.SetBinContent(i+1, uncert_pp);
        fwd_pt_pb.SetBinContent(i+1, uncert_pb);
    }
    // Save results of loop 1 and 2
    out_pt.cd();
    mid_pt.SetName("mid");
    mid_pt_pp.SetName("mid_pp");
    mid_pt_pb.SetName("mid_pb");

    fwd_pt.SetName("fwd");
    fwd_pt_pp.SetName("fwd_pp");
    fwd_pt_pb.SetName("fwd_pb");
    
    mid_pt.Write();
    mid_pt_pp.Write();
    mid_pt_pb.Write();
    
    fwd_pt.Write();
    fwd_pt_pp.Write();
    fwd_pt_pb.Write();
    
    mid_pt.SetName("mid_pt"); // Restore names
    mid_pt_pp.SetName("mid_pt_pp");
    mid_pt_pb.SetName("mid_pt_pb");

    fwd_pt.SetName("fwd_pt");
    fwd_pt_pp.SetName("fwd_pt_pp");
    fwd_pt_pb.SetName("fwd_pt_pb");
    out_pt.Close();
 

    // Start loop3 - mid_cent
    out_name = "./syst_roots/syst_cent_" + syst_type + ".root";
    TFile out_cent(out_name.c_str(), "recreate");

    h_pp_mid = "hAccPt_2021_midy_Int";
    h_pb_mid = "hAccPt_2021_midy_Int";
    h_pb_nomi = (TH1D *)pb_nomi_input->Get(h_pb_mid.Data());
    h_pb_syst = (TH1D *)pb_syst_input->Get(h_pb_mid.Data());
    h_pp_nomi = (TH1D *)pp_nomi_input->Get(h_pp_mid.Data());
    h_pp_syst = (TH1D *)pp_syst_input->Get(h_pp_mid.Data());

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0, 10, 20, 30, 40, 50, 90};
	TH1D mid_cent("mid_cent", "mid", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_pb("mid_cent_pb", "mid_pb", NBINS_mid_cent, edges_mid_cent);

	// pp nominal
	double n_pp_nomi = h_pp_nomi->GetBinContent(1);

	// pp syst
	double n_pp_syst = h_pp_syst->GetBinContent(1);

	// Pb nominal
	double n_Pb_nomi = h_pb_nomi->GetBinContent(1);

	// Pb syst
	double n_Pb_syst = h_pb_syst->GetBinContent(1);


	// Compute uncertainty
	// Proppt
	double uncert_mid = compute_uncertainty(n_pp_nomi, n_pp_syst, n_Pb_nomi, n_Pb_syst);
    double uncert_mid_pb = compute_uncertainty_pb(n_Pb_nomi, n_Pb_syst);
    

	// Non-prompt
	//cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

	// Fill histograms
	for(int i = 0; i < NBINS_mid_cent; i++){
		mid_cent.SetBinContent(i+1, uncert_mid); // i starts from 0, hist elements starts from 1
        mid_cent_pb.SetBinContent(i+1, uncert_mid_pb);
	}

	// Start loop4 - fwd_cent
	h_pp_fwd = "hAccPt_2021_Fory_Int"; 
	h_pb_fwd = "hAccPt_2021_Fory_Int"; 
    h_pb_nomi = (TH1D *)pb_nomi_input->Get(h_pb_fwd.Data());
    h_pb_syst = (TH1D *)pb_syst_input->Get(h_pb_fwd.Data());
    h_pp_nomi = (TH1D *)pp_nomi_input->Get(h_pp_fwd.Data());
    h_pp_syst = (TH1D *)pp_syst_input->Get(h_pp_fwd.Data());

    const int NBINS_fwd_cent = 6;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 20, 30, 40 , 50, 90};
    TH1D fwd_cent("fwd_cent", "fwd", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_pb("fwd_cent_pb", "fwd_pb", NBINS_fwd_cent, edges_fwd_cent);

	// pp nominal
	n_pp_nomi = h_pp_nomi->GetBinContent(1);

	// pp syst
	n_pp_syst = h_pp_syst->GetBinContent(1);

	// Pb nominal
	n_Pb_nomi = h_pb_nomi->GetBinContent(1);

	// Pb syst
	n_Pb_syst = h_pb_syst->GetBinContent(1);


	// Compute uncertainty
	// Proppt
	double uncert_fwd = compute_uncertainty(n_pp_nomi, n_pp_syst, n_Pb_nomi, n_Pb_syst);
    double uncert_fwd_pb = compute_uncertainty_pb(n_Pb_nomi, n_Pb_syst);
    

	// Non-prompt
	//cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

	// Fill histograms
	for(int i = 0; i < NBINS_fwd_cent; i++){
		fwd_cent.SetBinContent(i+1, uncert_fwd); // i starts from 0, hist elements starts from 1
        fwd_cent_pb.SetBinContent(i+1, uncert_fwd_pb);
	}

    // Save results of loop 3 and 4
    out_cent.cd();
    mid_cent.SetName("mid");
    mid_cent_pb.SetName("mid_pb");

    fwd_cent.SetName("fwd");
    fwd_cent_pb.SetName("fwd_pb");
    
    mid_cent.Write();
    mid_cent_pb.Write();
    
    fwd_cent.Write();
    fwd_cent_pb.Write();
    
    mid_cent.SetName("mid_cent"); // Restore names
    mid_cent_pb.SetName("mid_cent_pb");

    fwd_cent.SetName("fwd_cent");
    fwd_cent_pb.SetName("fwd_cent_pb");
    out_cent.Close();
}

double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst)
{
    // sqrt of (diff/n_PbPb)^2 + (diff/n_PP)^2
    return TMath::Sqrt(TMath::Power(((pb_nomi-pb_syst)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_syst)/pp_nomi), 2));
}

double compute_uncertainty_pp(double pp_nomi, double pp_syst)
{
    return TMath::Sqrt(TMath::Power(((pp_nomi-pp_syst)/pp_nomi), 2));
}

double compute_uncertainty_pb(double pb_nomi, double pb_syst)
{
    return TMath::Sqrt(TMath::Power(((pb_nomi-pb_syst)/pb_nomi), 2));
}

  
    // Open file
    // Read N_jpsi and b_fraction
    // Compute # of PR and NP
    // Compute difference btw nominal and syst.
    // Fill the containers (n_PR and n_NP)
        


    // Save or Print as root file
    // hist[0,1]
