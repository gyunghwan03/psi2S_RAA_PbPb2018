#include "TH1.h"
#include "TFile.h"
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
double compute_uncertainty(double PR_nomi, double PR_syst, double NP_nomi, double NP_syst);
double compute_uncertainty(double nominal_value, double syst_value);

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
	//TString PR_nomi = "acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20241127_fwdNew.root";
	TString PR_nomi = "acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20250123.root";
	//TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root";
	//TString pp_nomi = "acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240314.root";
	TString NP_nomi = "acceptance_NonPrompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20250123.root";
    //TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20231019.root";
    //TString pp_nomi_NP = "mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW1_tnp1_20230728.root";
    //TString pb_nomi_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20231019.root";
    //TString pb_nomi_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2s_PtW1_tnp1_new_20231019.root";

    TString PR_syst = "acceptance_Prompt_Jpsi_GenOnly_wgt0_pp_SysUp0_20250114.root";
    TString NP_syst = "acceptance_NonPrompt_Jpsi_GenOnly_wgt0_pp_SysUp0_20250114.root";
    

    TString h_PR_mid = "hAccPt_2021_midy";
    TString h_NP_mid = "hAccPt_2021_midy";
    

    // Define pointers for input files
    TFile *PR_nomi_input = new TFile(path+PR_nomi.Data());
    TFile *NP_nomi_input = new TFile(path+NP_nomi.Data());
    TFile *PR_syst_input = new TFile(path+PR_syst.Data());
    TFile *NP_syst_input = new TFile(path+NP_syst.Data());
    
    TH1D *h_NP_nomi = (TH1D *)NP_nomi_input->Get(h_NP_mid.Data());
    TH1D *h_NP_syst = (TH1D *)NP_syst_input->Get(h_NP_mid.Data());
    TH1D *h_PR_nomi = (TH1D *)PR_nomi_input->Get(h_PR_mid.Data());
    TH1D *h_PR_syst = (TH1D *)PR_syst_input->Get(h_PR_mid.Data());


    // Start loop
        // loop1 - mid_pt
        // loop2 - fwd_pt, results of loop 1 and 2 are saved into syst_pt_[something].root
        // loop3 - mid_cent
        // loop4 - fwd_cent, results of loop 3 and 4 are saved into syst_cent_[something].root
    
    // Start loop1 - mid_pt
    string out_name = "./syst_roots/syst_pt_" + syst_type + ".root";
    TFile out_pt(out_name.c_str(), "recreate");

    const int NBINS_mid_pt = 6;
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    TH1D mid_pt("mid_pt", "mid", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_PR("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);

    for (int i = 0; i < NBINS_mid_pt; i++) {
        // Get number of PR and NP Jpsi
        // PR nominal
        double n_PR_nomi = h_PR_nomi->GetBinContent(i+1);
        
        // PR syst
        double n_PR_syst = h_PR_syst->GetBinContent(i+1);
        
        // Pb nominal
        double n_NP_nomi = h_NP_nomi->GetBinContent(i+1);
        
        // Pb syst
        double n_NP_syst = h_NP_syst->GetBinContent(i+1);


        // Compute uncertainty
        // ProPRt
        double uncert = compute_uncertainty(n_PR_nomi, n_PR_syst, n_NP_nomi, n_NP_syst);
        double uncert_PR = compute_uncertainty(n_PR_nomi, n_PR_syst);
        double uncert_NP = compute_uncertainty(n_NP_nomi, n_NP_syst);
        
        
        // Non-prompt
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_pt.SetBinContent(i+1, uncert); // i starts from 0, hist elements starts from 1
        mid_pt_PR.SetBinContent(i+1, uncert_PR);
        mid_pt_NP.SetBinContent(i+1, uncert_NP);
    }

    // loop2 - fwd_pt
    TString h_PR_fwd = "hAccPt_2021_Fory";
    TString h_NP_fwd = "hAccPt_2021_Fory";

    h_NP_nomi = (TH1D *)NP_nomi_input->Get(h_NP_fwd.Data());
    h_NP_syst = (TH1D *)NP_syst_input->Get(h_NP_fwd.Data());
    h_PR_nomi = (TH1D *)PR_nomi_input->Get(h_PR_fwd.Data());
    h_PR_syst = (TH1D *)PR_syst_input->Get(h_PR_fwd.Data());


    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 50};
    TH1D fwd_pt("fwd_pt", "fwd", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_PR("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    for (int i = 0; i < NBINS_fwd_pt; i++) {
         // Open input files
        // Get number of PR and NP Jpsi
        // PR nominal
        double n_PR_nomi = h_PR_nomi->GetBinContent(i+1);
        
        // PR syst
        double n_PR_syst = h_PR_syst->GetBinContent(i+1);
        
        // Pb nominal
        double n_Pb_nomi = h_NP_nomi->GetBinContent(i+1);

        // Pb syst
        double n_Pb_syst = h_NP_syst->GetBinContent(i+1);


        // Compute uncertainty
        // ProPRt
        double uncert = compute_uncertainty(n_PR_nomi, n_PR_syst, n_Pb_nomi, n_Pb_syst);
        double uncert_PR = compute_uncertainty(n_PR_nomi, n_PR_syst);
        double uncert_NP = compute_uncertainty(n_Pb_nomi, n_Pb_syst);
        
        
        // Non-prompt
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        fwd_pt.SetBinContent(i+1, uncert); // i starts from 0, hist elements starts from 1
        fwd_pt_PR.SetBinContent(i+1, uncert_PR);
        fwd_pt_NP.SetBinContent(i+1, uncert_NP);
    }
    // Save results of loop 1 and 2
    out_pt.cd();
    mid_pt.SetName("mid");
    mid_pt_PR.SetName("mid_PR");
    mid_pt_NP.SetName("mid_NP");

    fwd_pt.SetName("fwd");
    fwd_pt_PR.SetName("fwd_PR");
    fwd_pt_NP.SetName("fwd_NP");
    
    mid_pt.Write();
    mid_pt_PR.Write();
    mid_pt_NP.Write();
    
    fwd_pt.Write();
    fwd_pt_PR.Write();
    fwd_pt_NP.Write();
    
    mid_pt.SetName("mid_pt"); // Restore names
    mid_pt_PR.SetName("mid_pt_PR");
    mid_pt_NP.SetName("mid_pt_NP");

    fwd_pt.SetName("fwd_pt");
    fwd_pt_PR.SetName("fwd_pt_PR");
    fwd_pt_NP.SetName("fwd_pt_NP");
    out_pt.Close();
 

    // Start loop3 - mid_cent
    out_name = "./syst_roots/syst_cent_" + syst_type + ".root";
    TFile out_cent(out_name.c_str(), "recreate");

    h_PR_mid = "hAccPt_2021_midy_Int";
    h_NP_mid = "hAccPt_2021_midy_Int";
    h_NP_nomi = (TH1D *)NP_nomi_input->Get(h_NP_mid.Data());
    h_NP_syst = (TH1D *)NP_syst_input->Get(h_NP_mid.Data());
    h_PR_nomi = (TH1D *)PR_nomi_input->Get(h_PR_mid.Data());
    h_PR_syst = (TH1D *)PR_syst_input->Get(h_PR_mid.Data());

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0,10,20,30,40,50,90};
	TH1D mid_cent("mid_cent", "mid", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_PR("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);

	// PR nominal
	double n_PR_nomi = h_PR_nomi->GetBinContent(1);

	// PR syst
	double n_PR_syst = h_PR_syst->GetBinContent(1);

	// Pb nominal
	double n_Pb_nomi = h_NP_nomi->GetBinContent(1);

	// Pb syst
	double n_Pb_syst = h_NP_syst->GetBinContent(1);


	// Compute uncertainty
	// ProPRt
	double uncert_mid = compute_uncertainty(n_PR_nomi, n_PR_syst, n_Pb_nomi, n_Pb_syst);
    double uncert_mid_PR = compute_uncertainty(n_PR_nomi, n_PR_syst);
    double uncert_mid_NP = compute_uncertainty(n_Pb_nomi, n_Pb_syst);
    

	// Non-prompt
	//cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

	// Fill histograms
	for(int i = 0; i < NBINS_mid_cent; i++){
		mid_cent.SetBinContent(i+1, uncert_mid); // i starts from 0, hist elements starts from 1
        mid_cent_PR.SetBinContent(i+1, uncert_mid_PR);
        mid_cent_NP.SetBinContent(i+1, uncert_mid_NP);
	}

	// Start loop4 - fwd_cent
	h_PR_fwd = "hAccPt_2021_Fory_Int"; 
	h_NP_fwd = "hAccPt_2021_Fory_Int"; 
    h_NP_nomi = (TH1D *)NP_nomi_input->Get(h_NP_fwd.Data());
    h_NP_syst = (TH1D *)NP_syst_input->Get(h_NP_fwd.Data());
    h_PR_nomi = (TH1D *)PR_nomi_input->Get(h_PR_fwd.Data());
    h_PR_syst = (TH1D *)PR_syst_input->Get(h_PR_fwd.Data());

    const int NBINS_fwd_cent = 4;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 30, 50, 90};
    TH1D fwd_cent("fwd_cent", "fwd", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_PR("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);

	// PR nominal
	n_PR_nomi = h_PR_nomi->GetBinContent(1);

	// PR syst
	n_PR_syst = h_PR_syst->GetBinContent(1);

	// Pb nominal
	n_Pb_nomi = h_NP_nomi->GetBinContent(1);

	// Pb syst
	n_Pb_syst = h_NP_syst->GetBinContent(1);


	// Compute uncertainty
	// ProPRt
	double uncert_fwd = compute_uncertainty(n_PR_nomi, n_PR_syst, n_Pb_nomi, n_Pb_syst);
    double uncert_fwd_PR = compute_uncertainty(n_PR_nomi, n_PR_syst);
    double uncert_fwd_NP = compute_uncertainty(n_Pb_nomi, n_Pb_syst);
    

	// Non-prompt
	//cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

	// Fill histograms
	for(int i = 0; i < NBINS_fwd_cent; i++){
		fwd_cent.SetBinContent(i+1, uncert_fwd); // i starts from 0, hist elements starts from 1
        fwd_cent_PR.SetBinContent(i+1, uncert_fwd_PR);
        fwd_cent_NP.SetBinContent(i+1, uncert_fwd_NP);
	}

    // Save results of loop 3 and 4
    out_cent.cd();
    mid_cent.SetName("mid");
    mid_cent_PR.SetName("mid_PR");
    mid_cent_NP.SetName("mid_NP");

    fwd_cent.SetName("fwd");
    fwd_cent_PR.SetName("fwd_PR");
    fwd_cent_NP.SetName("fwd_NP");
    
    mid_cent.Write();
    mid_cent_PR.Write();
    mid_cent_NP.Write();
    
    fwd_cent.Write();
    fwd_cent_PR.Write();
    fwd_cent_NP.Write();
    
    mid_cent.SetName("mid_cent"); // Restore names
    mid_cent_PR.SetName("mid_cent_PR");
    mid_cent_NP.SetName("mid_cent_NP");

    fwd_cent.SetName("fwd_cent");
    fwd_cent_PR.SetName("fwd_cent_PR");
    fwd_cent_NP.SetName("fwd_cent_NP");
    out_cent.Close();
}

double compute_uncertainty(double PR_nomi, double PR_syst, double NP_nomi, double NP_syst)
{
    // sqrt of (diff/n_PbPb)^2 + (diff/n_PP)^2
    return TMath::Sqrt(TMath::Power(((NP_nomi-NP_syst)/NP_nomi), 2) + TMath::Power(((PR_nomi-PR_syst)/PR_nomi), 2));
}

double compute_uncertainty(double nominal_value, double syst_value)
{
    return TMath::Sqrt(TMath::Power(((nominal_value-syst_value)/nominal_value), 2));
}
  
    // Open file
    // Read N_jpsi and b_fraction
    // Compute # of PR and NP
    // Compute difference btw nominal and syst.
    // Fill the containers (n_PR and n_NP)
        


    // Save or Print as root file
    // hist[0,1]
