#include "TH1.h"
#include "TFile.h"
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst);
double compute_uncertainty(double pp_nomi, double pp_syst);


void syst_eff()
{
    gSystem->mkdir("syst_roots");
    string syst_type = "Eff";

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
	//TString pb_nomi_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW1_tnp1_20241121_ppWfunc.root";
	//TString pb_nomi_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW1_tnp1_20241121_ppWfunc.root";
	TString pb_nomi_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW1_tnp1_20250114_ppWfunc.root";
	TString pb_nomi_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW1_tnp1_20250114_ppWfunc.root";
	TString pp_nomi_PR = "mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW1_tnp1_20241011.root";
	TString pp_nomi_NP = "mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtW1_tnp1_20241011.root";

    TString pp_syst_PR = "mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW0_tnp1_20241119.root";
    TString pp_syst_NP = "mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtW0_tnp1_20241119.root";
    TString pb_syst_PR = "mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW0_tnp1_20250114_ppWfunc.root";
    TString pb_syst_NP = "mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW0_tnp1_20250114_ppWfunc.root";
    

    TString h_pp_mid = "mc_eff_vs_pt_TnP1_PtW%d_absy0_1p6";
    TString h_pb_mid = "mc_eff_vs_pt_TnP1_PtW%d_cent_0_to_180_absy0_1p6";
    

    // Define pointers for input files
    TFile *PR_pp_nomi_input = new TFile(Form(path+pp_nomi_PR.Data(),1));
    TFile *PR_pb_nomi_input = new TFile(Form(path+pb_nomi_PR.Data(),1));
    TFile *PR_pp_syst_input = new TFile(Form(path+pp_syst_PR.Data(),0));
    TFile *PR_pb_syst_input = new TFile(Form(path+pb_syst_PR.Data(),0));
    
    TFile *NP_pp_nomi_input = new TFile(Form(path+pp_nomi_NP.Data(),1));
    TFile *NP_pb_nomi_input = new TFile(Form(path+pb_nomi_NP.Data(),1));
    TFile *NP_pp_syst_input = new TFile(Form(path+pp_syst_NP.Data(),0));
    TFile *NP_pb_syst_input = new TFile(Form(path+pb_syst_NP.Data(),0));

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
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    TH1D mid_pt_PR("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_PR_pp("mid_pt_PR_pp", "mid_PR_pp", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_PR_pb("mid_pt_PR_pb", "mid_PR_pb", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP_pp("mid_pt_NP_pp", "mid_NP_pp", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP_pb("mid_pt_NP_pb", "mid_NP_pb", NBINS_mid_pt, edges_mid_pt);

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
		//cout << "#### Nominal ####" << " PR PbPb : " << n_PR_Pb_nomi << ", NP PbPb : " << n_NP_Pb_nomi << endl;
        
        // Pb syst
        double n_PR_Pb_syst= h_pb_syst_PR->GetBinContent(i+2);;
        double n_NP_Pb_syst= h_pb_syst_NP->GetBinContent(i+2);;
		//cout << "##### Syst. #####" << " PR PbPb : " << n_PR_Pb_syst << ", NP PbPb : " << n_NP_Pb_syst << endl;


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_pt_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        mid_pt_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        mid_pt_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        mid_pt_NP.SetBinContent(i+1, NP_uncert);
        mid_pt_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        mid_pt_NP_pb.SetBinContent(i+1, NP_uncert_pb);
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
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 40};
    TH1D fwd_pt_PR("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_PR_pp("fwd_pt_PR_pp", "fwd_PR_pp", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_PR_pb("fwd_pt_PR_pb", "fwd_PR_pb", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP_pp("fwd_pt_NP_pp", "fwd_NP_pp", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP_pb("fwd_pt_NP_pb", "fwd_NP_pb", NBINS_fwd_pt, edges_fwd_pt);
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
		//cout << "#### Nominal ####" << " PR PbPb : " << n_PR_Pb_nomi << ", NP PbPb : " << n_NP_Pb_nomi << endl;
        
        // Pb syst
        double n_PR_Pb_syst= h_pb_syst_PR->GetBinContent(i+2);;
        double n_NP_Pb_syst= h_pb_syst_NP->GetBinContent(i+2);;
		//cout << "##### Syst. #####" << " PR PbPb : " << n_PR_Pb_syst << ", NP PbPb : " << n_NP_Pb_syst << endl;


        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        fwd_pt_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        fwd_pt_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        fwd_pt_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        fwd_pt_NP.SetBinContent(i+1, NP_uncert);
        fwd_pt_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        fwd_pt_NP_pb.SetBinContent(i+1, NP_uncert_pb);
    }
    // Save results of loop 1 and 2
    out_pt.cd();
    mid_pt_PR.SetName("mid_PR");
    mid_pt_PR_pp.SetName("mid_PR_pp");
    mid_pt_PR_pb.SetName("mid_PR_pb");
    mid_pt_NP.SetName("mid_NP");
    mid_pt_NP_pp.SetName("mid_NP_pp");
    mid_pt_NP_pb.SetName("mid_NP_pb");

    fwd_pt_PR.SetName("fwd_PR");
    fwd_pt_PR_pp.SetName("fwd_PR_pp");
    fwd_pt_PR_pb.SetName("fwd_PR_pb");
    fwd_pt_NP.SetName("fwd_NP");
    fwd_pt_NP_pp.SetName("fwd_NP_pp");
    fwd_pt_NP_pb.SetName("fwd_NP_pb");
    
    mid_pt_PR.Write();
    mid_pt_PR_pp.Write();
    mid_pt_PR_pb.Write();
    mid_pt_NP.Write();
    mid_pt_NP_pp.Write();
    mid_pt_NP_pb.Write();
    
    fwd_pt_PR.Write();
    fwd_pt_PR_pp.Write();
    fwd_pt_PR_pb.Write();
    fwd_pt_NP.Write();
    fwd_pt_NP_pp.Write();
    fwd_pt_NP_pb.Write();
    
    mid_pt_PR.SetName("mid_pt_PR"); // Restore names
    mid_pt_PR_pp.SetName("mid_pt_PR_pp");
    mid_pt_PR_pb.SetName("mid_pt_PR_pb");

    mid_pt_NP.SetName("mid_pt_NP");
    mid_pt_NP_pp.SetName("mid_pt_NP_pp");
    mid_pt_NP_pb.SetName("mid_pt_NP_pb");

    fwd_pt_PR.SetName("fwd_pt_PR");
    fwd_pt_PR_pp.SetName("fwd_pt_PR_pp");
    fwd_pt_PR_pb.SetName("fwd_pt_PR_pb");

    fwd_pt_NP.SetName("fwd_pt_NP");
    fwd_pt_NP_pp.SetName("fwd_pt_NP_pp");
    fwd_pt_NP_pb.SetName("fwd_pt_NP_pb");
    out_pt.Close();
 

    // Start loop3 - mid_cent
    out_name = "./syst_roots/syst_cent_" + syst_type + ".root";
    TFile out_cent(out_name.c_str(), "recreate");

    h_pp_mid = "mc_eff_Integrated_TnP1_PtW%d_absy0_1p6";
    h_pb_mid = "mc_eff_vs_cent_TnP1_PtW%d_pt_6p5_to_40_absy0_1p6";
    h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_mid.Data(),0));

    h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_mid.Data(),1));
    h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_mid.Data(),0));
    h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_mid.Data(),1));
    h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_mid.Data(),0));
    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0,10,20,30,40,50,90};
    TH1D mid_cent_PR("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_PR_pp("mid_cent_PR_pp", "mid_PR_pp", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_PR_pb("mid_cent_PR_pb", "mid_PR_pb", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP_pp("mid_cent_NP_pp", "mid_NP_pp", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP_pb("mid_cent_NP_pb", "mid_NP_pb", NBINS_mid_cent, edges_mid_cent);

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
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_cent_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        mid_cent_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        mid_cent_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        mid_cent_NP.SetBinContent(i+1, NP_uncert);
        mid_cent_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        mid_cent_NP_pb.SetBinContent(i+1, NP_uncert_pb);
    }

    // Start loop4 - fwd_cent
    h_pp_fwd = "mc_eff_Integrated_TnP1_PtW%d_absy1p6_2p4";
    h_pb_fwd = "mc_eff_vs_cent_TnP1_PtW%d_pt_3_to_40_absy1p6_2p4";
    h_pb_nomi_PR = (TH1D *)PR_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_PR = (TH1D *)PR_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_PR = (TH1D *)PR_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_PR = (TH1D *)PR_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    h_pb_nomi_NP = (TH1D *)NP_pb_nomi_input->Get(Form(h_pb_fwd.Data(),1));
    h_pb_syst_NP = (TH1D *)NP_pb_syst_input->Get(Form(h_pb_fwd.Data(),0));
    h_pp_nomi_NP = (TH1D *)NP_pp_nomi_input->Get(Form(h_pp_fwd.Data(),1));
    h_pp_syst_NP = (TH1D *)NP_pp_syst_input->Get(Form(h_pp_fwd.Data(),0));

    const int NBINS_fwd_cent = 4;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 30, 50, 90};
    TH1D fwd_cent_PR("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_PR_pp("fwd_cent_PR_pp", "fwd_PR_pp", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_PR_pb("fwd_cent_PR_pb", "fwd_PR_pb", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP_pp("fwd_cent_NP_pp", "fwd_NP_pp", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP_pb("fwd_cent_NP_pb", "fwd_NP_pb", NBINS_fwd_cent, edges_fwd_cent);

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

		cout << "#### Nominal ####" << " PR PbPb : " << n_PR_Pb_nomi << ", NP PbPb : " << n_NP_Pb_nomi << endl;
		cout << "##### Syst. #####" << " PR PbPb : " << n_PR_Pb_syst << ", NP PbPb : " << n_NP_Pb_syst << endl;

        // Compute uncertainty
        // Proppt
        double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst, n_PR_Pb_nomi, n_PR_Pb_syst);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_syst);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_syst);
        
        // Non-prompt
        double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst, n_NP_Pb_nomi, n_NP_Pb_syst);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_syst);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_syst);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        fwd_cent_PR.SetBinContent(i+1, PR_uncert); // i starts from 0, hist elements starts from 1
        fwd_cent_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        fwd_cent_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        fwd_cent_NP.SetBinContent(i+1, NP_uncert);
        fwd_cent_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        fwd_cent_NP_pb.SetBinContent(i+1, NP_uncert_pb);
    }
    // Save results of loop 3 and 4
    out_cent.cd();
    mid_cent_PR.SetName("mid_PR");
    mid_cent_PR_pp.SetName("mid_PR_pp");
    mid_cent_PR_pb.SetName("mid_PR_pb");
    mid_cent_NP.SetName("mid_NP");
    mid_cent_NP_pp.SetName("mid_NP_pp");
    mid_cent_NP_pb.SetName("mid_NP_pb");

    fwd_cent_PR.SetName("fwd_PR");
    fwd_cent_PR_pp.SetName("fwd_PR_pp");
    fwd_cent_PR_pb.SetName("fwd_PR_pb");
    fwd_cent_NP.SetName("fwd_NP");
    fwd_cent_NP_pp.SetName("fwd_NP_pp");
    fwd_cent_NP_pb.SetName("fwd_NP_pb");
    
    mid_cent_PR.Write();
    mid_cent_PR_pp.Write();
    mid_cent_PR_pb.Write();
    mid_cent_NP.Write();
    mid_cent_NP_pp.Write();
    mid_cent_NP_pb.Write();
    
    fwd_cent_PR.Write();
    fwd_cent_PR_pp.Write();
    fwd_cent_PR_pb.Write();
    fwd_cent_NP.Write();
    fwd_cent_NP_pp.Write();
    fwd_cent_NP_pb.Write();
    
    mid_cent_PR.SetName("mid_cent_PR"); // Restore names
    mid_cent_PR_pp.SetName("mid_cent_PR_pp");
    mid_cent_PR_pb.SetName("mid_cent_PR_pb");

    mid_cent_NP.SetName("mid_cent_NP");
    mid_cent_NP_pp.SetName("mid_cent_NP_pp");
    mid_cent_NP_pb.SetName("mid_cent_NP_pb");

    fwd_cent_PR.SetName("fwd_cent_PR");
    fwd_cent_PR_pp.SetName("fwd_cent_PR_pp");
    fwd_cent_PR_pb.SetName("fwd_cent_PR_pb");

    fwd_cent_NP.SetName("fwd_cent_NP");
    fwd_cent_NP_pp.SetName("fwd_cent_NP_pp");
    fwd_cent_NP_pb.SetName("fwd_cent_NP_pb");
    out_cent.Close();
}

double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst)
{
    // sqrt of (diff/n_PbPb)^2 + (diff/n_PP)^2
    return TMath::Sqrt(TMath::Power(((pb_nomi-pb_syst)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_syst)/pp_nomi), 2));
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
