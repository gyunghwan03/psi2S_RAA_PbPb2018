#include "TH1.h"
#include "TFile.h"
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst);


void syst_eff()
{
    gSystem->mkdir("syst_roots");
    string syst_type = "Eff";

    // Path label
    TString path = "../../Eff_Acc/roots/";
    TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root";
    TString pp_nomi_NP = "mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW1_tnp1_20230728.root";
    TString pb_nomi_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20230729.root";
    TString pb_nomi_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2s_PtW1_tnp1_20230729.root";

    TString pp_syst_PR = "mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW0_tnp1_20231002.root";
    TString pp_syst_NP = "mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW0_tnp1_20231002.root";
    TString pb_syst_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW0_tnp1_20231002.root";
    TString pb_syst_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2s_PtW0_tnp1_new_20231002.root";

    TString h_pp_mid = "mc_eff_vs_pt_TnP1_PtW%d_absy0_1p6";
    TString h_pb_mid = "mc_eff_vs_pt_TnP1_PtW%d_cent_0_to_180_absy0_1p6";

    // Define pointers for input files
    TFile *PR_pp_nomi_input =  new TFile(path+pp_nomi_PR.Data());
    TFile *PR_pb_nomi_input = new TFile(path+pb_nomi_PR.Data());
    TFile *PR_pp_syst_input = new TFile(path+pp_syst_PR.Data());
    TFile *PR_pb_syst_input = new TFile(path+pb_syst_PR.Data());
    
    TFile *NP_pp_nomi_input =  new TFile(path+pp_nomi_NP.Data());
    TFile *NP_pb_nomi_input = new TFile(path+pb_nomi_NP.Data());
    TFile *NP_pp_syst_input = new TFile(path+pp_syst_NP.Data());
    TFile *NP_pb_syst_input = new TFile(path+pb_syst_NP.Data());

    TH1D *h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    TH1D *h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    TH1D *h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    TH1D *h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_mid.Data(),0));

    TH1D *h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    TH1D *h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    TH1D *h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    TH1D *h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_mid.Data(),0));

    // Start loop
        // loop1 - mid_pt
        // loop2 - fwd_pt, results of loop 1 and 2 are saved into syst_pt_[something].root
        // loop3 - mid_cent
        // loop4 - fwd_cent, results of loop 3 and 4 are saved into syst_cent_[something].root
    
    // Start loop1 - mid_pt
    string out_name = "./syst_roots/syst_pt_" + syst_type + ".root";
    TFile out_pt(out_name.c_str(), "recreate");

    const int NBINS_mid_pt = 6;
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25,  50};
    TH1D mid_pt_PR("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);

    for (int i = 0; i < NBINS_mid_pt; i++) {
        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi = h_pp_nomi_PR->GetBinContent(i+2);
        double n_NP_pp_nomi = h_pp_nomi_NP->GetBinContent(i+2);
        
        // pp syst
        double n_PR_pp_syst = h_pp_syst_PR->GetBinContent(i+2);
        double n_NP_pp_syst = h_pp_syst_NP->GetBinContent(i+2);
        
        // Pb nominal
        double n_PR_Pb_nomi = h_pb_nomi_PR->GetBinContent(i+2);
        double n_NP_Pb_nomi = h_pb_nomi_NP->GetBinContent(i+2);
        
        // Pb syst
        double n_PR_Pb_syst= h_pb_syst_PR->GetBinContent(i+2);;
        double n_NP_Pb_syst= h_pb_syst_NP->GetBinContent(i+2);;


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_pt_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        mid_pt_NP.SetBinContent(i+1, NP_uncert); 
    }

    // loop2 - fwd_pt
    TString h_pp_fwd = "mc_eff_vs_pt_TnP1_PtW%d_absy1p6_2p4";
    TString h_pb_fwd = "mc_eff_vs_pt_TnP1_PtW%d_cent_0_to_180_absy1p6_2p4;";

    h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 5, 6.5, 12, 50};
    TH1D fwd_pt_PR("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    for (int i = 0; i < pp_fwd_pt.size(); i++) {
         // Open input files
        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi = h_pp_nomi_PR->GetBinContent(i+2);
        double n_NP_pp_nomi = h_pp_nomi_NP->GetBinContent(i+2);
        
        // pp syst
        double n_PR_pp_syst = h_pp_syst_PR->GetBinContent(i+2);
        double n_NP_pp_syst = h_pp_syst_NP->GetBinContent(i+2);
        
        // Pb nominal
        double n_PR_Pb_nomi = h_pb_nomi_PR->GetBinContent(i+2);
        double n_NP_Pb_nomi = h_pb_nomi_NP->GetBinContent(i+2);
        
        // Pb syst
        double n_PR_Pb_syst= h_pb_syst_PR->GetBinContent(i+2);;
        double n_NP_Pb_syst= h_pb_syst_NP->GetBinContent(i+2);;


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

        // Fill histograms
        fwd_pt_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        fwd_pt_NP.SetBinContent(i+1, NP_uncert);
    }
    // Save results of loop 1 and 2
    out_pt.cd();
    mid_pt_PR.SetName("mid_PR");
    mid_pt_NP.SetName("mid_NP");
    fwd_pt_PR.SetName("fwd_PR");
    fwd_pt_NP.SetName("fwd_NP");
    mid_pt_PR.Write();
    mid_pt_NP.Write();
    fwd_pt_PR.Write();
    fwd_pt_NP.Write();
    mid_pt_PR.SetName("mid_pt_PR"); // Restore names
    mid_pt_NP.SetName("mid_pt_NP");
    fwd_pt_PR.SetName("fwd_pt_PR");
    fwd_pt_NP.SetName("fwd_pt_NP");
    out_pt.Close();
 

    // Start loop3 - mid_cent
    out_name = "./syst_roots/syst_cent_" + syst_type + ".root";
    TFile out_cent(out_name.c_str(), "recreate");

    h_pp_mid = "mc_eff_Integrated_TnP1_PtW%d_absy0_1p6";
    h_pb_mid = "mc_eff_vs_cent_TnP1_PtW%d_pt_6p5_to_50_absy0_1p6";
    h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_mid.Data(),0));

    h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_mid.Data(),0));
    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    TH1D mid_cent_PR("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);

    for (int i = 0; i < NBINS_mid_cent; i++) {
        // pp nominal
        double n_PR_pp_nomi = h_pp_nomi_PR->GetBinContent(1);
        double n_NP_pp_nomi = h_pp_nomi_NP->GetBinContent(1);
        
        // pp syst
        double n_PR_pp_syst = h_pp_syst_PR->GetBinContent(1);
        double n_NP_pp_syst = h_pp_syst_NP->GetBinContent(1);
        
        // Pb nominal
        double n_PR_Pb_nomi = h_pb_nomi_PR->GetBinContent(i+1);
        double n_NP_Pb_nomi = h_pb_nomi_NP->GetBinContent(i+1);
        
        // Pb syst
        double n_PR_Pb_syst = h_pb_syst_PR->GetBinContent(i+1);
        double n_NP_Pb_syst = h_pb_syst_NP->GetBinContent(i+1);


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_cent_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        mid_cent_NP.SetBinContent(i+1, NP_uncert);
    }

    // Start loop4 - fwd_cent
    h_pp_fwd = "mc_eff_Integrated_TnP1_PtW%d_absy1p6_2p4";
    h_pb_fwd = "mc_eff_vs_cent_TnP1_PtW%d_pt_3_to_50_absy1p6_2p4";
    h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    const int NBINS_fwd_cent = 6;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    TH1D fwd_cent_PR("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);

    for (int i = 0; i < NBINS_fwd_cent; i++) {
        // pp nominal
        double n_PR_pp_nomi = h_pp_nomi_PR->GetBinContent(1);
        double n_NP_pp_nomi = h_pp_nomi_NP->GetBinContent(1);
        
        // pp syst
        double n_PR_pp_syst = h_pp_syst_PR->GetBinContent(1);
        double n_NP_pp_syst = h_pp_syst_NP->GetBinContent(1);
        
        // Pb nominal
        double n_PR_Pb_nomi = h_pb_nomi_PR->GetBinContent(i+1);
        double n_NP_Pb_nomi = h_pb_nomi_NP->GetBinContent(i+1);
        
        // Pb syst
        double n_PR_Pb_syst = h_pb_syst_PR->GetBinContent(i+1);
        double n_NP_Pb_syst = h_pb_syst_NP->GetBinContent(i+1);


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

        // Fill histograms
        fwd_cent_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        fwd_cent_NP.SetBinContent(i+1, NP_uncert);
    }
    // Save results of loop 3 and 4
    out_cent.cd();
    mid_cent_PR.SetName("mid_PR");
    mid_cent_NP.SetName("mid_NP");
    fwd_cent_PR.SetName("fwd_PR");
    fwd_cent_NP.SetName("fwd_NP");
    mid_cent_PR.Write();
    mid_cent_NP.Write();
    fwd_cent_PR.Write();
    fwd_cent_NP.Write();
    mid_cent_PR.SetName("mid_cent_PR"); // Restore names
    mid_cent_NP.SetName("mid_cent_NP");
    fwd_cent_PR.SetName("fwd_cent_PR");
    fwd_cent_NP.SetName("fwd_cent_NP");
    out_cent.Close();
}

double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst)
{
    // sqrt of (diff/n_PbPb)^2 + (diff/n_PP)^2
    return TMath::Sqrt(TMath::Power(((pb_nomi-pb_syst)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_syst)/pp_nomi), 2));
}


  
    // Open file
    // Read N_jpsi and b_fraction
    // Compute # of PR and NP
    // Compute difference btw nominal and syst.
    // Fill the containers (n_PR and n_NP)
        


    // Save or Print as root file
    // hist[0,1]
