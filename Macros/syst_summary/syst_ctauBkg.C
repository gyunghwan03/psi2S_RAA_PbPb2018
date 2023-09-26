#include <vector>
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
double compute_uncertainty(double pp_nomi, double pp_syst, double pb_nomi, double pb_syst);


void syst_ctauBkg()
{
    gSystem->mkdir("syst_roots");
    string syst_type = "ctauBkg";

    // Path label
    string nominal_path_pp = "../pp_psi2S_230512/roots/2DFit_No_Weight/Final/";
    string nominal_path_pb = "../psi2S_230512/roots/2DFit_No_Weight/Final/";
    string syst_path_pp = "../pp_psi2S_230512/systematic/" + syst_type + "/roots/2DFit_No_Weight/Final/";
    string syst_path_pb = "../psi2S_230512/systematics/" + syst_type + "/roots/2DFit_No_Weight/Final/";

    // Define pointers for input files
    TFile *pp_nominal_input = nullptr;
    TFile *pb_nominal_input = nullptr;
    TFile *pp_syst_input = nullptr;
    TFile *pb_syst_input = nullptr;    
    

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
    TH1D mid_pt_PR("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);

    for (int i = 0; i < pp_mid_pt.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_mid_pt[i].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = nominal_path_pb + pb_mid_pt[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp + pp_mid_pt[i];
        pp_syst_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pb + pb_mid_pt[i];
        pb_syst_input = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_syst;
        double n_NP_pp_syst;
        compute_n_jpsi(pp_syst_input, n_PR_pp_syst, n_NP_pp_syst);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_syst;
        double n_NP_Pb_syst;
        compute_n_jpsi(pb_syst_input, n_PR_Pb_syst, n_NP_Pb_syst);


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
    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 5, 6.5, 12, 50};
    TH1D fwd_pt_PR("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    for (int i = 0; i < pp_fwd_pt.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_fwd_pt[i].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = nominal_path_pb + pb_fwd_pt[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp + pp_fwd_pt[i];
        pp_syst_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pb + pb_fwd_pt[i];
        pb_syst_input = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_syst;
        double n_NP_pp_syst;
        compute_n_jpsi(pp_syst_input, n_PR_pp_syst, n_NP_pp_syst);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_syst;
        double n_NP_Pb_syst;
        compute_n_jpsi(pb_syst_input, n_PR_Pb_syst, n_NP_Pb_syst);


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

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    TH1D mid_cent_PR("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);

    for (int i = 0; i < pb_mid_cent.size(); i++) {
        // Open input files
        // pp_mid_cent has only one elements
        string temp_input_path = nominal_path_pp + pp_mid_cent[0].c_str(); 
        pp_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = nominal_path_pb + pb_mid_cent[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp + pp_mid_cent[0];
        pp_syst_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pb + pb_mid_cent[i];
        pb_syst_input = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_syst;
        double n_NP_pp_syst;
        compute_n_jpsi(pp_syst_input, n_PR_pp_syst, n_NP_pp_syst);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_syst;
        double n_NP_Pb_syst;
        compute_n_jpsi(pb_syst_input, n_PR_Pb_syst, n_NP_Pb_syst);


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
    const int NBINS_fwd_cent = 6;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    TH1D fwd_cent_PR("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);

    for (int i = 0; i < pb_fwd_cent.size(); i++) {
        // Open input files
        // pp_fwd_cent has only one elements
        string temp_input_path = nominal_path_pp + pp_fwd_cent[0].c_str(); 
        pp_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = nominal_path_pb + pb_fwd_cent[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp + pp_fwd_cent[0];
        pp_syst_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pb + pb_fwd_cent[i];
        pb_syst_input = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_syst;
        double n_NP_pp_syst;
        compute_n_jpsi(pp_syst_input, n_PR_pp_syst, n_NP_pp_syst);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_syst;
        double n_NP_Pb_syst;
        compute_n_jpsi(pb_syst_input, n_PR_Pb_syst, n_NP_Pb_syst);


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



void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP)
{
    //Read N_jpsi, b_fraction
    auto fit_result = (RooFitResult*)my_file->Get("fitresult_pdfCTAUMASS_Tot_dsToFit");
    auto hist_frac = (TH1D*)my_file->Get("2DfitResults"); // b_fraction -> only first bin is used
    auto RooReal_n_jpsi = (RooRealVar*)fit_result->constPars().find("N_Jpsi");
    double n_jpsi = RooReal_n_jpsi->getVal();
    double b_frac = hist_frac->GetBinContent(1);

    // Compute # of PR and NP
    n_NP = n_jpsi * b_frac;
    n_PR = n_jpsi - n_NP;
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
