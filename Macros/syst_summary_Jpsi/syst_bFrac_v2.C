#include <vector>
#include "input_list.h" // Input file list
using namespace std;


void compute_n_jpsi(TFile *my_file, double &n_PR, double &n_NP);
//double compute_uncertainty(double pp_nomi, double pp_alpha, double pp_n, double pp_f, double pp_x, double pb_nomi, double pb_alpha, double pb_n, double pb_f, double pb_x);
double compute_uncertainty(double nomi, double alpha, double n_value, double f_value, double x_value);
double compute_uncertainty(double pp_syst, double pb_syst);


void syst_bFrac_v2()
{
    gSystem->mkdir("syst_roots");
    string syst_type = "bFrac";

    
    // Path label
    string nominal_path_pp = "../pp_Jpsi/roots/2DFit_No_Weight/Final/";
    string nominal_path_pb = "../Jpsi/roots/2DFit_No_Weight/Final/";
    string syst_path_pp_Err = "../pp_Jpsi/systematics/ctauERR/roots/2DFit_No_Weight/Final/";
    string syst_path_pb_Err = "../Jpsi/systematics/ctauERR/roots/2DFit_No_Weight/Final/";
    string syst_path_pp_Res = "../pp_Jpsi/systematics/ctauRES/roots/2DFit_No_Weight/Final/";
    string syst_path_pb_Res = "../Jpsi/systematics/ctauRES/roots/2DFit_No_Weight/Final/";
    string syst_path_pp_Bkg = "../pp_Jpsi/systematics/ctauBKG/roots/2DFit_No_Weight/Final/";
    string syst_path_pb_Bkg = "../Jpsi/systematics/ctauBKG/roots/2DFit_No_Weight/Final/";
    string syst_path_pp_True = "../pp_Jpsi/systematics/ctauTRUE/roots/2DFit_No_Weight/Final/";
    string syst_path_pb_True = "../Jpsi/systematics/ctauTRUE/roots/2DFit_No_Weight/Final/";
    // ============== CAUTION ============== //
    // Variables below have suffixes such that alpha, n, f, and x
    // But this code compute b-fraction uncertainty (Changed input pathes above)
    // (alpha -> Err, n -> Res, f -> Bkg, x -> True)
    // ============== End CAUTION ============== //


    // Define pointers for input files
    TFile *pp_nominal_input = nullptr;
    TFile *pb_nominal_input = nullptr;
    TFile *pp_syst_input_alpha = nullptr;
    TFile *pb_syst_input_alpha = nullptr;
    TFile *pp_syst_input_n = nullptr;
    TFile *pb_syst_input_n = nullptr;
    TFile *pp_syst_input_f = nullptr;
    TFile *pb_syst_input_f = nullptr;
    TFile *pp_syst_input_x = nullptr;
    TFile *pb_syst_input_x = nullptr;
    TFile *pb_nominal_int_input = nullptr;
    TFile *pb_syst_int_input_alpha = nullptr;
    TFile *pb_syst_int_input_n = nullptr;
    TFile *pb_syst_int_input_f = nullptr;
    TFile *pb_syst_int_input_x = nullptr;


    // Start loop
        // loop1 - mid_pt
        // loop2 - fwd_pt, results of loop 1 and 2 are saved into syst_pt_[something].root
        // loop3 - mid_cent
        // loop4 - fwd_cent, results of loop 3 and 4 are saved into syst_cent_[something].root
    
    // Start loop1 - mid_pt
    string out_name = "./syst_roots/syst_pt_" + syst_type + ".root";
    TFile out_pt(out_name.c_str(), "recreate");
    cout << "Mid pt" << endl;

    const int NBINS_mid_pt = 6;
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    TH1D mid_pt_PR("mid_pt_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_PR_pp("mid_pt_PR_pp", "mid_PR_pp", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_PR_pb("mid_pt_PR_pb", "mid_PR_pb", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP("mid_pt_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP_pp("mid_pt_NP_pp", "mid_NP_pp", NBINS_mid_pt, edges_mid_pt);
    TH1D mid_pt_NP_pb("mid_pt_NP_pb", "mid_NP_pb", NBINS_mid_pt, edges_mid_pt);

	TH1D mid_pt_Err("mid_pt_Err" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_Res("mid_pt_Res" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_Bkg("mid_pt_Bkg" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_True("mid_pt_True" , "", NBINS_mid_pt, edges_mid_pt);

	TH1D mid_pt_PR_Err_pp("mid_pt_PR_Err_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_Res_pp("mid_pt_PR_Res_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_Bkg_pp("mid_pt_PR_Bkg_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_True_pp("mid_pt_PR_True_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_Err_pb("mid_pt_Err_PR_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_Res_pb("mid_pt_Res_PR_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_Bkg_pb("mid_pt_Bkg_PR_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_PR_True_pb("mid_pt_True_PR_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Err_pp("mid_pt_Err_NP_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Res_pp("mid_pt_Res_NP_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Bkg_pp("mid_pt_Bkg_NP_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_True_pp("mid_pt_True_NP_pp" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Err_pb("mid_pt_Err_NP_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Res_pb("mid_pt_Res_NP_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_Bkg_pb("mid_pt_Bkg_NP_pb" , "", NBINS_mid_pt, edges_mid_pt);
	TH1D mid_pt_NP_True_pb("mid_pt_True_NP_pb" , "", NBINS_mid_pt, edges_mid_pt);

    for (int i = 0; i < pp_mid_pt.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_mid_pt[i].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_mid_pt[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Err + pp_mid_pt[i];
        pp_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_mid_pt[i];
        pb_syst_input_alpha = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Res + pp_mid_pt[i];
        pp_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_mid_pt[i];
        pb_syst_input_n = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Bkg + pp_mid_pt[i];
        pp_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_mid_pt[i];
        pb_syst_input_f = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_True + pp_mid_pt[i];
        pp_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_mid_pt[i];
        pb_syst_input_x = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_alpha;
        double n_NP_pp_alpha;
        double n_PR_pp_n;
        double n_NP_pp_n;
        double n_PR_pp_f;
        double n_NP_pp_f;
        double n_PR_pp_x;
        double n_NP_pp_x;
        compute_n_jpsi(pp_syst_input_alpha, n_PR_pp_alpha, n_NP_pp_alpha);
        compute_n_jpsi(pp_syst_input_n, n_PR_pp_n, n_NP_pp_n);
        compute_n_jpsi(pp_syst_input_f, n_PR_pp_f, n_NP_pp_f);
        compute_n_jpsi(pp_syst_input_x, n_PR_pp_x, n_NP_pp_x);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_alpha;
        double n_NP_Pb_alpha;
        double n_PR_Pb_n;
        double n_NP_Pb_n;
        double n_PR_Pb_f;
        double n_NP_Pb_f;
        double n_PR_Pb_x;
        double n_NP_Pb_x;
        compute_n_jpsi(pb_syst_input_alpha, n_PR_Pb_alpha, n_NP_Pb_alpha);
        compute_n_jpsi(pb_syst_input_n, n_PR_Pb_n, n_NP_Pb_n);
        compute_n_jpsi(pb_syst_input_f, n_PR_Pb_f, n_NP_Pb_f);
        compute_n_jpsi(pb_syst_input_x, n_PR_Pb_x, n_NP_Pb_x);


        // Compute uncertainty
        // Proppt
        //double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert = compute_uncertainty(PR_uncert_pp, PR_uncert_pb);

        double PR_uncertErr_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertErr_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertRes_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_n, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertRes_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_n, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertBkg_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_f, n_PR_pp_nomi);
        double PR_uncertBkg_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_f, n_PR_Pb_nomi);
        double PR_uncertTrue_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_x);
        double PR_uncertTrue_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_x);
        
        // Non-prompt
        //double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert = compute_uncertainty(NP_uncert_pp, NP_uncert_pb);
        
        double NP_uncertErr_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertErr_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertRes_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_n, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertRes_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_n, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertBkg_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_f, n_NP_pp_nomi);
        double NP_uncertBkg_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_f, n_NP_Pb_nomi);
        double NP_uncertTrue_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_x);
        double NP_uncertTrue_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_x);
        
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
		//cout << "pp nomi : " << n_PR_pp_nomi << "\tPb nomi : " << n_PR_Pb_nomi  << "\tpp Err Syst : " << n_PR_pp_alpha << "\tPb Err Syst: " << n_PR_Pb_alpha << "\tPR Err pp : " << PR_uncertErr_pp << "\tPR Err Pb :" << PR_uncertErr_pb << endl;
        //printf("pt %.1f - %.1f\n", edges_mid_pt[i], edges_mid_pt[i + 1]);
        //cout << fixed << setw(5) << "Pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << n_NP_Pb_alpha << "\tRes Syst : " << n_NP_Pb_n << "\tPb Bkg Syst : " << n_NP_Pb_f << "\tNP True Pb : " << n_NP_Pb_x << "\tNP Uncert : " << NP_uncert_pb << endl;
        //cout << "pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << NP_uncertErr_pb << "\tRes Syst : " << NP_uncertRes_pb << "\tpb Bkg Syst : " << NP_uncertBkg_pb << "\tNP True pb : " << NP_uncertTrue_pb << "\tNP Uncert : " << NP_uncert_pb << endl;

        // Fill histograms
        mid_pt_PR.SetBinContent(i+1, PR_uncert); // The i starts from 0, hist elements start from 1
        mid_pt_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        mid_pt_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        mid_pt_NP.SetBinContent(i+1, NP_uncert);
        mid_pt_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        mid_pt_NP_pb.SetBinContent(i+1, NP_uncert_pb);

		mid_pt_PR_Err_pp.SetBinContent(i+1, PR_uncertErr_pp);
		mid_pt_PR_Res_pp.SetBinContent(i+1, PR_uncertRes_pp);
		mid_pt_PR_Bkg_pp.SetBinContent(i+1, PR_uncertBkg_pp);
		mid_pt_PR_True_pp.SetBinContent(i+1, PR_uncertTrue_pp);
		mid_pt_PR_Err_pb.SetBinContent(i+1, PR_uncertErr_pb);
		mid_pt_PR_Res_pb.SetBinContent(i+1, PR_uncertRes_pb);
		mid_pt_PR_Bkg_pb.SetBinContent(i+1, PR_uncertBkg_pb);
		mid_pt_PR_True_pb.SetBinContent(i+1, PR_uncertTrue_pb);
		mid_pt_NP_Err_pp.SetBinContent(i+1, NP_uncertErr_pp);
		mid_pt_NP_Res_pp.SetBinContent(i+1, NP_uncertRes_pp);
		mid_pt_NP_Bkg_pp.SetBinContent(i+1, NP_uncertBkg_pp);
		mid_pt_NP_True_pp.SetBinContent(i+1, NP_uncertTrue_pp);
		mid_pt_NP_Err_pb.SetBinContent(i+1, NP_uncertErr_pb);
		mid_pt_NP_Res_pb.SetBinContent(i+1, NP_uncertRes_pb);
		mid_pt_NP_Bkg_pb.SetBinContent(i+1, NP_uncertBkg_pb);
		mid_pt_NP_True_pb.SetBinContent(i+1, NP_uncertTrue_pb);

    }

    cout << "Fwd pt" << endl;

    // loop2 - fwd_pt
    const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 40};
    TH1D fwd_pt_PR("fwd_pt_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_PR_pp("fwd_pt_PR_pp", "fwd_PR_pp", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_PR_pb("fwd_pt_PR_pb", "fwd_PR_pb", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP("fwd_pt_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP_pp("fwd_pt_NP_pp", "fwd_NP_pp", NBINS_fwd_pt, edges_fwd_pt);
    TH1D fwd_pt_NP_pb("fwd_pt_NP_pb", "fwd_NP_pb", NBINS_fwd_pt, edges_fwd_pt);
    
	TH1D fwd_pt_PR_Err_pp("fwd_pt_PR_Err_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_Res_pp("fwd_pt_PR_Res_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_Bkg_pp("fwd_pt_PR_Bkg_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_True_pp("fwd_pt_PR_True_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_Err_pb("fwd_pt_Err_PR_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_Res_pb("fwd_pt_Res_PR_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_Bkg_pb("fwd_pt_Bkg_PR_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_PR_True_pb("fwd_pt_True_PR_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Err_pp("fwd_pt_Err_NP_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Res_pp("fwd_pt_Res_NP_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Bkg_pp("fwd_pt_Bkg_NP_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_True_pp("fwd_pt_True_NP_pp" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Err_pb("fwd_pt_Err_NP_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Res_pb("fwd_pt_Res_NP_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_Bkg_pb("fwd_pt_Bkg_NP_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
	TH1D fwd_pt_NP_True_pb("fwd_pt_True_NP_pb" , "", NBINS_fwd_pt, edges_fwd_pt);
    for (int i = 0; i < pp_fwd_pt.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_fwd_pt[i].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_fwd_pt[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Err + pp_fwd_pt[i];
        pp_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_fwd_pt[i];
        pb_syst_input_alpha = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Res + pp_fwd_pt[i];
        pp_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_fwd_pt[i];
        pb_syst_input_n = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Bkg + pp_fwd_pt[i];
        pp_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_fwd_pt[i];
        pb_syst_input_f = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_True + pp_fwd_pt[i];
        pp_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_fwd_pt[i];
        pb_syst_input_x = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_alpha;
        double n_NP_pp_alpha;
        double n_PR_pp_n;
        double n_NP_pp_n;
        double n_PR_pp_f;
        double n_NP_pp_f;
        double n_PR_pp_x;
        double n_NP_pp_x;
        compute_n_jpsi(pp_syst_input_alpha, n_PR_pp_alpha, n_NP_pp_alpha);
        compute_n_jpsi(pp_syst_input_n, n_PR_pp_n, n_NP_pp_n);
        compute_n_jpsi(pp_syst_input_f, n_PR_pp_f, n_NP_pp_f);
        compute_n_jpsi(pp_syst_input_x, n_PR_pp_x, n_NP_pp_x);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        
        // Pb syst
        double n_PR_Pb_alpha;
        double n_NP_Pb_alpha;
        double n_PR_Pb_n;
        double n_NP_Pb_n;
        double n_PR_Pb_f;
        double n_NP_Pb_f;
        double n_PR_Pb_x;
        double n_NP_Pb_x;
        compute_n_jpsi(pb_syst_input_alpha, n_PR_Pb_alpha, n_NP_Pb_alpha);
        compute_n_jpsi(pb_syst_input_n, n_PR_Pb_n, n_NP_Pb_n);
        compute_n_jpsi(pb_syst_input_f, n_PR_Pb_f, n_NP_Pb_f);
        compute_n_jpsi(pb_syst_input_x, n_PR_Pb_x, n_NP_Pb_x);        

        // Compute uncertainty
        // Proppt
        //double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert = compute_uncertainty(PR_uncert_pp, PR_uncert_pb);

        double PR_uncertErr_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertErr_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertRes_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_n, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertRes_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_n, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertBkg_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_f, n_PR_pp_nomi);
        double PR_uncertBkg_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_f, n_PR_Pb_nomi);
        double PR_uncertTrue_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_x);
        double PR_uncertTrue_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_x);

        
        // Non-prompt
        //double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert = compute_uncertainty(NP_uncert_pp, NP_uncert_pb);

        double NP_uncertErr_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertErr_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertRes_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_n, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertRes_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_n, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertBkg_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_f, n_NP_pp_nomi);
        double NP_uncertBkg_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_f, n_NP_Pb_nomi);
        double NP_uncertTrue_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_x);
        double NP_uncertTrue_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_x);
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
		//cout << "pp nomi : " << n_PR_pp_nomi << "\tPb nomi : " << n_PR_Pb_nomi  << "\tpp Err Syst : " << n_PR_pp_alpha << "\tPb Err Syst: " << n_PR_Pb_alpha << "\tPR Err pp : " << PR_uncertErr_pp << "\tPR Err Pb :" << PR_uncertErr_pb << endl;

        //printf("pt %.1f - %.1f\n", edges_fwd_pt[i], edges_fwd_pt[i + 1]);
        //cout << fixed << setw(5) << "Pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << n_NP_Pb_alpha << "\tRes Syst : " << n_NP_Pb_n << "\tPb Bkg Syst : " << n_NP_Pb_f << "\tNP True Pb : " << n_NP_Pb_x << "\tNP Uncert : " << NP_uncert_pb << endl;
        //cout << "pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << NP_uncertErr_pb << "\tRes Syst : " << NP_uncertRes_pb << "\tpb Bkg Syst : " << NP_uncertBkg_pb << "\tNP True pb : " << NP_uncertTrue_pb << "\tNP Uncert : " << NP_uncert_pb << endl;
        //cout << fixed << setw(5) << "pp nomi : " << n_NP_pp_nomi << "\tErr Syst : " << n_NP_pp_alpha << "\tRes Syst : " << n_NP_pp_n << "\tpp Bkg Syst : " << n_NP_pp_f << "\tNP True pp : " << n_NP_pp_x << "\tNP Uncert : " << NP_uncert_pp <<  endl;
		//cout << "pp nomi : " << n_NP_pp_nomi << "\tErr Syst : " << NP_uncertErr_pp << "\tRes Syst : " << NP_uncertRes_pp << "\tpp Bkg Syst : " << NP_uncertBkg_pp << "\tNP True pp : " << NP_uncertTrue_pp << "\tNP Uncert : " << NP_uncert_pp <<  endl;

        // Fill histograms
        fwd_pt_PR.SetBinContent(i+1, PR_uncert); // The i starts from 0, hist elements start from 1
        fwd_pt_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        fwd_pt_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        fwd_pt_NP.SetBinContent(i+1, NP_uncert);
        fwd_pt_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        fwd_pt_NP_pb.SetBinContent(i+1, NP_uncert_pb);

		fwd_pt_PR_Err_pp.SetBinContent(i+1, PR_uncertErr_pp);
		fwd_pt_PR_Res_pp.SetBinContent(i+1, PR_uncertRes_pp);
		fwd_pt_PR_Bkg_pp.SetBinContent(i+1, PR_uncertBkg_pp);
		fwd_pt_PR_True_pp.SetBinContent(i+1, PR_uncertTrue_pp);
		fwd_pt_PR_Err_pb.SetBinContent(i+1, PR_uncertErr_pb);
		fwd_pt_PR_Res_pb.SetBinContent(i+1, PR_uncertRes_pb);
		fwd_pt_PR_Bkg_pb.SetBinContent(i+1, PR_uncertBkg_pb);
		fwd_pt_PR_True_pb.SetBinContent(i+1, PR_uncertTrue_pb);
		fwd_pt_NP_Err_pp.SetBinContent(i+1, NP_uncertErr_pp);
		fwd_pt_NP_Res_pp.SetBinContent(i+1, NP_uncertRes_pp);
		fwd_pt_NP_Bkg_pp.SetBinContent(i+1, NP_uncertBkg_pp);
		fwd_pt_NP_True_pp.SetBinContent(i+1, NP_uncertTrue_pp);
		fwd_pt_NP_Err_pb.SetBinContent(i+1, NP_uncertErr_pb);
		fwd_pt_NP_Res_pb.SetBinContent(i+1, NP_uncertRes_pb);
		fwd_pt_NP_Bkg_pb.SetBinContent(i+1, NP_uncertBkg_pb);
		fwd_pt_NP_True_pb.SetBinContent(i+1, NP_uncertTrue_pb);
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

	fwd_pt_PR_Err_pp.Write();
	fwd_pt_PR_Res_pp.Write();
	fwd_pt_PR_Bkg_pp.Write();
	fwd_pt_PR_True_pp.Write();
	fwd_pt_PR_Err_pb.Write();
	fwd_pt_PR_Res_pb.Write();
	fwd_pt_PR_Bkg_pb.Write();
	fwd_pt_PR_True_pb.Write();
	fwd_pt_NP_Err_pp.Write();
	fwd_pt_NP_Res_pp.Write();
	fwd_pt_NP_Bkg_pp.Write();
	fwd_pt_NP_True_pp.Write();
	fwd_pt_NP_Err_pb.Write();
	fwd_pt_NP_Res_pb.Write();
	fwd_pt_NP_Bkg_pb.Write();
	fwd_pt_NP_True_pb.Write();

	mid_pt_PR_Err_pp.Write();
	mid_pt_PR_Res_pp.Write();
	mid_pt_PR_Bkg_pp.Write();
	mid_pt_PR_True_pp.Write();
	mid_pt_PR_Err_pb.Write();
	mid_pt_PR_Res_pb.Write();
	mid_pt_PR_Bkg_pb.Write();
	mid_pt_PR_True_pb.Write();
	mid_pt_NP_Err_pp.Write();
	mid_pt_NP_Res_pp.Write();
	mid_pt_NP_Bkg_pp.Write();
	mid_pt_NP_True_pp.Write();
	mid_pt_NP_Err_pb.Write();
	mid_pt_NP_Res_pb.Write();
	mid_pt_NP_Bkg_pb.Write();
	mid_pt_NP_True_pb.Write();
    
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

    out_name = "./syst_roots/syst_int_" + syst_type + ".root";
    TFile out_int(out_name.c_str(), "recreate");

    const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0,10,20,30,40,50,90};
    TH1D mid_cent_PR("mid_cent_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_PR_pp("mid_cent_PR_pp", "mid_PR_pp", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_PR_pb("mid_cent_PR_pb", "mid_PR_pb", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP("mid_cent_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP_pp("mid_cent_NP_pp", "mid_NP_pp", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_cent_NP_pb("mid_cent_NP_pb", "mid_NP_pb", NBINS_mid_cent, edges_mid_cent);
	
	TH1D mid_cent_PR_Err_pp("mid_cent_PR_Err_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_Res_pp("mid_cent_PR_Res_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_Bkg_pp("mid_cent_PR_Bkg_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_True_pp("mid_cent_PR_True_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_Err_pb("mid_cent_Err_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_Res_pb("mid_cent_Res_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_Bkg_pb("mid_cent_Bkg_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_PR_True_pb("mid_cent_True_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Err_pp("mid_cent_Err_NP_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Res_pp("mid_cent_Res_NP_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Bkg_pp("mid_cent_Bkg_NP_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_True_pp("mid_cent_True_NP_pp" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Err_pb("mid_cent_Err_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Res_pb("mid_cent_Res_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_Bkg_pb("mid_cent_Bkg_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
	TH1D mid_cent_NP_True_pb("mid_cent_True_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);

    TH1D mid_int_PR("mid_int_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_PR_pb("mid_int_PR_pb", "mid_PR_pb", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_NP("mid_int_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_NP_pb("mid_int_NP_pb", "mid_NP_pb", NBINS_mid_cent, edges_mid_cent);

    TH1D mid_int_PR_Err_pb("mid_int_Err_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_PR_Res_pb("mid_int_Res_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_PR_Bkg_pb("mid_int_Bkg_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_PR_True_pb("mid_int_True_PR_pb" , "", NBINS_mid_cent, edges_mid_cent);

    TH1D mid_int_NP_Err_pb("mid_int_Err_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_NP_Res_pb("mid_int_Res_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_NP_Bkg_pb("mid_int_Bkg_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);
    TH1D mid_int_NP_True_pb("mid_int_True_NP_pb" , "", NBINS_mid_cent, edges_mid_cent);

    cout << "Mid Cent" << endl;

    for (int i = 0; i < pb_mid_cent.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_mid_cent[0].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_mid_cent[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_mid_int[0].c_str();
        pb_nominal_int_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Err + pp_mid_cent[0];
        pp_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_mid_cent[i];
        pb_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_mid_int[0];
        pb_syst_int_input_alpha = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Res + pp_mid_cent[0];
        pp_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_mid_cent[i];
        pb_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_mid_int[0];
        pb_syst_int_input_n = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Bkg + pp_mid_cent[0];
        pp_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_mid_cent[i];
        pb_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_mid_int[0];
        pb_syst_int_input_f = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_True + pp_mid_cent[0];
        pp_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_mid_cent[i];
        pb_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_mid_int[0];
        pb_syst_int_input_x = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_alpha;
        double n_NP_pp_alpha;
        double n_PR_pp_n;
        double n_NP_pp_n;
        double n_PR_pp_f;
        double n_NP_pp_f;
        double n_PR_pp_x;
        double n_NP_pp_x;
        compute_n_jpsi(pp_syst_input_alpha, n_PR_pp_alpha, n_NP_pp_alpha);
        compute_n_jpsi(pp_syst_input_n, n_PR_pp_n, n_NP_pp_n);
        compute_n_jpsi(pp_syst_input_f, n_PR_pp_f, n_NP_pp_f);
        compute_n_jpsi(pp_syst_input_x, n_PR_pp_x, n_NP_pp_x);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        double n_PR_Pb_int_nomi;
        double n_NP_Pb_int_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);
        compute_n_jpsi(pb_nominal_int_input, n_PR_Pb_int_nomi, n_NP_Pb_int_nomi);
        
        // Pb syst
        double n_PR_Pb_alpha;
        double n_NP_Pb_alpha;
        double n_PR_Pb_n;
        double n_NP_Pb_n;
        double n_PR_Pb_f;
        double n_NP_Pb_f;
        double n_PR_Pb_x;
        double n_NP_Pb_x;

        double n_PR_Pb_int_alpha;
        double n_NP_Pb_int_alpha;
        double n_PR_Pb_int_n;
        double n_NP_Pb_int_n;
        double n_PR_Pb_int_f;
        double n_NP_Pb_int_f;
        double n_PR_Pb_int_x;
        double n_NP_Pb_int_x;

        compute_n_jpsi(pb_syst_input_alpha, n_PR_Pb_alpha, n_NP_Pb_alpha);
        compute_n_jpsi(pb_syst_input_n, n_PR_Pb_n, n_NP_Pb_n);
        compute_n_jpsi(pb_syst_input_f, n_PR_Pb_f, n_NP_Pb_f);
        compute_n_jpsi(pb_syst_input_x, n_PR_Pb_x, n_NP_Pb_x);

        compute_n_jpsi(pb_syst_int_input_alpha, n_PR_Pb_int_alpha, n_NP_Pb_int_alpha);
        compute_n_jpsi(pb_syst_int_input_n, n_PR_Pb_int_n, n_NP_Pb_int_n);
        compute_n_jpsi(pb_syst_int_input_f, n_PR_Pb_int_f, n_NP_Pb_int_f);
        compute_n_jpsi(pb_syst_int_input_x, n_PR_Pb_int_x, n_NP_Pb_int_x);




        // Compute uncertainty
        // Proppt
        //double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert = compute_uncertainty(PR_uncert_pp, PR_uncert_pb);
        cout << "PR_uncert: " << PR_uncert << endl;

        //double PR_uncert_int = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_n, n_PR_Pb_int_f, n_PR_Pb_int_x);
        double PR_uncert_int_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_n, n_PR_Pb_int_f, n_PR_Pb_int_x);
        double PR_uncert_int = compute_uncertainty(PR_uncert_pp, PR_uncert_int_pb);

        double PR_uncertErr_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertErr_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertRes_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_n, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertRes_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_n, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertBkg_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_f, n_PR_pp_nomi);
        double PR_uncertBkg_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_f, n_PR_Pb_nomi);
        double PR_uncertTrue_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_x);
        double PR_uncertTrue_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_x);

        double PR_uncert_int_Err_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi);
        double PR_uncert_int_Res_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_n, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi);
        double PR_uncert_int_Bkg_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_f, n_PR_Pb_int_nomi);
        double PR_uncert_int_True_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_x);
        
        // Non-prompt
        //double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert = compute_uncertainty(NP_uncert_pp, NP_uncert_pb);

        //double NP_uncert_int = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_n, n_NP_Pb_int_f, n_NP_Pb_int_x);
        double NP_uncert_int_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_n, n_NP_Pb_int_f, n_NP_Pb_int_x);
        double NP_uncert_int = compute_uncertainty(NP_uncert_pp, NP_uncert_int_pb);

        double NP_uncertErr_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertErr_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertRes_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_n, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertRes_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_n, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertBkg_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_f, n_NP_pp_nomi);
        double NP_uncertBkg_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_f, n_NP_Pb_nomi);
        double NP_uncertTrue_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_x);
        double NP_uncertTrue_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_x);

        double NP_uncert_int_Err_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi);
        double NP_uncert_int_Res_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_n, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi);
        double NP_uncert_int_Bkg_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_f, n_NP_Pb_int_nomi);
        double NP_uncert_int_True_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_x);
        

        //printf("Cent %.f - %.f\n", edges_mid_cent[i] ,edges_mid_cent[i+1]);
		//cout << fixed << setw(5) << "Pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << n_NP_Pb_alpha << "\tRes Syst : " << n_NP_Pb_n << "\tPb Bkg Syst : " << n_NP_Pb_f << "\tNP True Pb : " << n_NP_Pb_x << "\tNP Uncert : " << NP_uncert_pb <<  endl;
		//cout << "pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << NP_uncertErr_pb << "\tRes Syst : " << NP_uncertRes_pb << "\tpb Bkg Syst : " << NP_uncertBkg_pb << "\tNP True pb : " << NP_uncertTrue_pb << "\tNP Uncert : " << NP_uncert_pb <<  endl;
        
        cout << "pb int nomi : " << n_PR_Pb_int_nomi << "\tErr Syst : " << n_PR_Pb_int_alpha << "\tRes Syst : " << n_PR_Pb_int_n << "\tPb Bkg Syst : " << n_PR_Pb_int_f << "\tPR True Pb : " << n_PR_Pb_int_x << "\tPR Uncert : " << PR_uncert_int_pb << endl;
        cout << "pb int nomi : " << n_NP_Pb_int_nomi << "\tErr Syst : " << n_NP_Pb_int_alpha << "\tRes Syst : " << n_NP_Pb_int_n << "\tPb Bkg Syst : " << n_NP_Pb_int_f << "\tNP True Pb : " << n_NP_Pb_int_x << "\tNP Uncert : " << NP_uncert_int_pb << endl;
        //cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;
        
        // Fill histograms
        mid_cent_PR.SetBinContent(i+1, PR_uncert); // The i starts from 0, hist elements start from 1
        mid_cent_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        mid_cent_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        mid_cent_NP.SetBinContent(i+1, NP_uncert);
        mid_cent_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        mid_cent_NP_pb.SetBinContent(i+1, NP_uncert_pb);

		mid_cent_PR_Err_pp.SetBinContent(i+1, PR_uncertErr_pp);
		mid_cent_PR_Res_pp.SetBinContent(i+1, PR_uncertRes_pp);
		mid_cent_PR_Bkg_pp.SetBinContent(i+1, PR_uncertBkg_pp);
		mid_cent_PR_True_pp.SetBinContent(i+1, PR_uncertTrue_pp);
		mid_cent_PR_Err_pb.SetBinContent(i+1, PR_uncertErr_pb);
		mid_cent_PR_Res_pb.SetBinContent(i+1, PR_uncertRes_pb);
		mid_cent_PR_Bkg_pb.SetBinContent(i+1, PR_uncertBkg_pb);
		mid_cent_PR_True_pb.SetBinContent(i+1, PR_uncertTrue_pb);
		mid_cent_NP_Err_pp.SetBinContent(i+1, NP_uncertErr_pp);
		mid_cent_NP_Res_pp.SetBinContent(i+1, NP_uncertRes_pp);
		mid_cent_NP_Bkg_pp.SetBinContent(i+1, NP_uncertBkg_pp);
		mid_cent_NP_True_pp.SetBinContent(i+1, NP_uncertTrue_pp);
		mid_cent_NP_Err_pb.SetBinContent(i+1, NP_uncertErr_pb);
		mid_cent_NP_Res_pb.SetBinContent(i+1, NP_uncertRes_pb);
		mid_cent_NP_Bkg_pb.SetBinContent(i+1, NP_uncertBkg_pb);
		mid_cent_NP_True_pb.SetBinContent(i+1, NP_uncertTrue_pb);

        mid_int_PR.SetBinContent(i+1, PR_uncert_int);
        mid_int_PR_pb.SetBinContent(i+1, PR_uncert_int_pb);
        mid_int_NP.SetBinContent(i+1, NP_uncert_int);
        mid_int_NP_pb.SetBinContent(i+1, NP_uncert_int_pb);

        mid_int_PR_Err_pb.SetBinContent(i+1, PR_uncert_int_Err_pb);
        mid_int_PR_Res_pb.SetBinContent(i+1, PR_uncert_int_Res_pb);
        mid_int_PR_Bkg_pb.SetBinContent(i+1, PR_uncert_int_Bkg_pb);
        mid_int_PR_True_pb.SetBinContent(i+1, PR_uncert_int_True_pb);
        mid_int_NP_Err_pb.SetBinContent(i+1, NP_uncert_int_Err_pb);
        mid_int_NP_Res_pb.SetBinContent(i+1, NP_uncert_int_Res_pb);
        mid_int_NP_Bkg_pb.SetBinContent(i+1, NP_uncert_int_Bkg_pb);
        mid_int_NP_True_pb.SetBinContent(i+1, NP_uncert_int_True_pb);

    }

    // Start loop4 - fwd_cent
    const int NBINS_fwd_cent = 4;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 30, 50, 90};
    TH1D fwd_cent_PR("fwd_cent_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_PR_pp("fwd_cent_PR_pp", "fwd_PR_pp", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_PR_pb("fwd_cent_PR_pb", "fwd_PR_pb", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP("fwd_cent_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP_pp("fwd_cent_NP_pp", "fwd_NP_pp", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_cent_NP_pb("fwd_cent_NP_pb", "fwd_NP_pb", NBINS_fwd_cent, edges_fwd_cent);

	TH1D fwd_cent_PR_Err_pp("fwd_cent_PR_Err_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_Res_pp("fwd_cent_PR_Res_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_Bkg_pp("fwd_cent_PR_Bkg_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_True_pp("fwd_cent_PR_True_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_Err_pb("fwd_cent_Err_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_Res_pb("fwd_cent_Res_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_Bkg_pb("fwd_cent_Bkg_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_PR_True_pb("fwd_cent_True_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Err_pp("fwd_cent_Err_NP_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Res_pp("fwd_cent_Res_NP_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Bkg_pp("fwd_cent_Bkg_NP_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_True_pp("fwd_cent_True_NP_pp" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Err_pb("fwd_cent_Err_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Res_pb("fwd_cent_Res_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_Bkg_pb("fwd_cent_Bkg_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
	TH1D fwd_cent_NP_True_pb("fwd_cent_True_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);

    TH1D fwd_int_PR("fwd_int_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_PR_pb("fwd_int_PR_pb", "fwd_PR_pb", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_NP("fwd_int_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_NP_pb("fwd_int_NP_pb", "fwd_NP_pb", NBINS_fwd_cent, edges_fwd_cent);

    TH1D fwd_int_PR_Err_pb("fwd_int_Err_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_PR_Res_pb("fwd_int_Res_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_PR_Bkg_pb("fwd_int_Bkg_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_PR_True_pb("fwd_int_True_PR_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    
    TH1D fwd_int_NP_Err_pb("fwd_int_Err_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_NP_Res_pb("fwd_int_Res_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_NP_Bkg_pb("fwd_int_Bkg_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);
    TH1D fwd_int_NP_True_pb("fwd_int_True_NP_pb" , "", NBINS_fwd_cent, edges_fwd_cent);

    cout << "Fwd Cent" << endl;

    for (int i = 0; i < pb_fwd_cent.size(); i++) {
         // Open input files
        string temp_input_path = nominal_path_pp + pp_fwd_cent[0].c_str();
        pp_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_fwd_cent[i].c_str();
        pb_nominal_input = TFile::Open(temp_input_path.c_str());
        temp_input_path = nominal_path_pb + pb_fwd_int[0].c_str();
        pb_nominal_int_input = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Err + pp_fwd_cent[0];
        pp_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_fwd_cent[i];
        pb_syst_input_alpha = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Err + pb_fwd_int[0];
        pb_syst_int_input_alpha = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Res + pp_fwd_cent[0];
        pp_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_fwd_cent[i];
        pb_syst_input_n = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Res + pb_fwd_int[0];
        pb_syst_int_input_n = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_Bkg + pp_fwd_cent[0];
        pp_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_fwd_cent[i];
        pb_syst_input_f = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_Bkg + pb_fwd_int[0];
        pb_syst_int_input_f = TFile::Open(temp_input_path.c_str());

        temp_input_path = syst_path_pp_True + pp_fwd_cent[0];
        pp_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_fwd_cent[i];
        pb_syst_input_x = TFile::Open(temp_input_path.c_str());
        temp_input_path = syst_path_pb_True + pb_fwd_int[0];
        pb_syst_int_input_x = TFile::Open(temp_input_path.c_str());


        // Get number of PR and NP Jpsi
        // pp nominal
        double n_PR_pp_nomi;
        double n_NP_pp_nomi;
        compute_n_jpsi(pp_nominal_input, n_PR_pp_nomi, n_NP_pp_nomi);
        //cout << "n_PR_pp_nomi: " << n_PR_pp_nomi << "\tn_NP_pp_nomi: " << n_NP_pp_nomi << endl;
        
        // pp syst
        double n_PR_pp_alpha;
        double n_NP_pp_alpha;
        double n_PR_pp_n;
        double n_NP_pp_n;
        double n_PR_pp_f;
        double n_NP_pp_f;
        double n_PR_pp_x;
        double n_NP_pp_x;
        compute_n_jpsi(pp_syst_input_alpha, n_PR_pp_alpha, n_NP_pp_alpha);
        compute_n_jpsi(pp_syst_input_n, n_PR_pp_n, n_NP_pp_n);
        compute_n_jpsi(pp_syst_input_f, n_PR_pp_f, n_NP_pp_f);
        compute_n_jpsi(pp_syst_input_x, n_PR_pp_x, n_NP_pp_x);
        
        // Pb nominal
        double n_PR_Pb_nomi;
        double n_NP_Pb_nomi;
        compute_n_jpsi(pb_nominal_input, n_PR_Pb_nomi, n_NP_Pb_nomi);

        double n_PR_Pb_int_nomi;
        double n_NP_Pb_int_nomi;
        compute_n_jpsi(pb_nominal_int_input, n_PR_Pb_int_nomi, n_NP_Pb_int_nomi);
        
        // Pb syst
        double n_PR_Pb_alpha;
        double n_NP_Pb_alpha;
        double n_PR_Pb_n;
        double n_NP_Pb_n;
        double n_PR_Pb_f;
        double n_NP_Pb_f;
        double n_PR_Pb_x;
        double n_NP_Pb_x;
        compute_n_jpsi(pb_syst_input_alpha, n_PR_Pb_alpha, n_NP_Pb_alpha);
        compute_n_jpsi(pb_syst_input_n, n_PR_Pb_n, n_NP_Pb_n);
        compute_n_jpsi(pb_syst_input_f, n_PR_Pb_f, n_NP_Pb_f);
        compute_n_jpsi(pb_syst_input_x, n_PR_Pb_x, n_NP_Pb_x);

        double n_PR_Pb_int_alpha;
        double n_NP_Pb_int_alpha;
        double n_PR_Pb_int_n;
        double n_NP_Pb_int_n;
        double n_PR_Pb_int_f;
        double n_NP_Pb_int_f;
        double n_PR_Pb_int_x;
        double n_NP_Pb_int_x;
        compute_n_jpsi(pb_syst_int_input_alpha, n_PR_Pb_int_alpha, n_NP_Pb_int_alpha);
        compute_n_jpsi(pb_syst_int_input_n, n_PR_Pb_int_n, n_NP_Pb_int_n);
        compute_n_jpsi(pb_syst_int_input_f, n_PR_Pb_int_f, n_NP_Pb_int_f);
        compute_n_jpsi(pb_syst_int_input_x, n_PR_Pb_int_x, n_NP_Pb_int_x);


        // Compute uncertainty
        // Prompt
        //double PR_uncert = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x);
        double PR_uncert_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_n, n_PR_Pb_f, n_PR_Pb_x);
        double PR_uncert = compute_uncertainty(PR_uncert_pp, PR_uncert_pb);
        double PR_uncertErr_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertErr_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_alpha, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertRes_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_n, n_PR_pp_nomi, n_PR_pp_nomi);
        double PR_uncertRes_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_n, n_PR_Pb_nomi, n_PR_Pb_nomi);
        double PR_uncertBkg_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_f, n_PR_pp_nomi);
        double PR_uncertBkg_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_f, n_PR_Pb_nomi);
        double PR_uncertTrue_pp = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_nomi, n_PR_pp_x);
        double PR_uncertTrue_pb = compute_uncertainty(n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_nomi, n_PR_Pb_x);

        //double PR_uncert_int = compute_uncertainty(n_PR_pp_nomi, n_PR_pp_alpha, n_PR_pp_n, n_PR_pp_f, n_PR_pp_x, n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_n, n_PR_Pb_int_f, n_PR_Pb_int_x);
        double PR_uncert_int_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_n, n_PR_Pb_int_f, n_PR_Pb_int_x);
        double PR_uncert_int = compute_uncertainty(PR_uncert_pp, PR_uncert_int_pb);
        double PR_uncert_int_Err_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_alpha, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi);
        double PR_uncert_int_Res_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_n, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi);
        double PR_uncert_int_Bkg_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_f, n_PR_Pb_int_nomi);
        double PR_uncert_int_True_pb = compute_uncertainty(n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_nomi, n_PR_Pb_int_x);
        
        // Non-prompt
        //double NP_uncert = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x);
        double NP_uncert_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_n, n_NP_Pb_f, n_NP_Pb_x);
        double NP_uncert = compute_uncertainty(NP_uncert_pp, NP_uncert_pb);
        double NP_uncertErr_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertErr_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_alpha, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertRes_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_n, n_NP_pp_nomi, n_NP_pp_nomi);
        double NP_uncertRes_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_n, n_NP_Pb_nomi, n_NP_Pb_nomi);
        double NP_uncertBkg_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_f, n_NP_pp_nomi);
        double NP_uncertBkg_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_f, n_NP_Pb_nomi);
        double NP_uncertTrue_pp = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_nomi, n_NP_pp_x);
        double NP_uncertTrue_pb = compute_uncertainty(n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_nomi, n_NP_Pb_x);

        //double NP_uncert_int = compute_uncertainty(n_NP_pp_nomi, n_NP_pp_alpha, n_NP_pp_n, n_NP_pp_f, n_NP_pp_x, n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_n, n_NP_Pb_int_f, n_NP_Pb_int_x);
        double NP_uncert_int_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_n, n_NP_Pb_int_f, n_NP_Pb_int_x);
        double NP_uncert_int = compute_uncertainty(NP_uncert_pp, NP_uncert_int_pb);
        double NP_uncert_int_Err_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_alpha, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi);
        double NP_uncert_int_Res_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_n, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi);
        double NP_uncert_int_Bkg_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_f, n_NP_Pb_int_nomi);
        double NP_uncert_int_True_pb = compute_uncertainty(n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_nomi, n_NP_Pb_int_x);

        // cout << "pp nomi : " << n_PR_pp_nomi << "\tPb nomi : " << n_PR_Pb_nomi  << "\tpp Err Syst : " << n_PR_pp_alpha << "\tPb Err Syst: " << n_PR_Pb_alpha << "\tPR Err pp : " << PR_uncertErr_pp << "\tPR Err Pb :" << PR_uncertErr_pb << endl;
        //cout << "pp nomi : " << n_PR_pp_nomi << "\tPb nomi : " << n_PR_Pb_nomi  << "\tpp Err Syst : " << n_PR_pp_x << "\tPb True Syst: " << n_PR_Pb_x << "\tPR True pp : " << PR_uncertTrue_pp << "\tPR True Pb :" << PR_uncertTrue_pb << endl;
        cout << "int nomi : " << n_NP_Pb_int_nomi << "\tint Err Syst : " << n_NP_Pb_int_alpha << "\tint Res Syst : " << n_NP_Pb_int_n << "\tint Bkg Syst : " << n_NP_Pb_int_f <<  "\tint True Syst: " << n_NP_Pb_int_x << "\tPb NP Uncert : "<< NP_uncert_int_pb <<  "\tNP Syst : " << NP_uncert_int << endl;
        // printf("Cent %.f - %.f\n", edges_fwd_cent[i] ,edges_fwd_cent[i+1]);
        // cout << fixed << setw(5) << "Pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << n_NP_Pb_alpha << "\tRes Syst : " << n_NP_Pb_n << "\tPb Bkg Syst : " << n_NP_Pb_f << "\tNP True Pb : " << n_NP_Pb_x << "\tNP Uncert : " << NP_uncert_pb <<  endl;
        // cout << "pb nomi : " << n_NP_Pb_nomi << "\tErr Syst : " << NP_uncertErr_pb << "\tRes Syst : " << NP_uncertRes_pb << "\tpb Bkg Syst : " << NP_uncertBkg_pb << "\tNP True pb : " << NP_uncertTrue_pb << "\tNP Uncert : " << NP_uncert_pb <<  endl;
        // cout << "PR_uncert: " << PR_uncert << "\tNP_uncert: " << NP_uncert << endl;

        // Fill histograms
        fwd_cent_PR.SetBinContent(i+1, PR_uncert); // The i starts from 0, hist elements start from 1
        fwd_cent_PR_pp.SetBinContent(i+1, PR_uncert_pp);
        fwd_cent_PR_pb.SetBinContent(i+1, PR_uncert_pb);
        fwd_cent_NP.SetBinContent(i+1, NP_uncert);
        fwd_cent_NP_pp.SetBinContent(i+1, NP_uncert_pp);
        fwd_cent_NP_pb.SetBinContent(i+1, NP_uncert_pb);

		fwd_cent_PR_Err_pp.SetBinContent(i+1, PR_uncertErr_pp);
		fwd_cent_PR_Res_pp.SetBinContent(i+1, PR_uncertRes_pp);
		fwd_cent_PR_Bkg_pp.SetBinContent(i+1, PR_uncertBkg_pp);
		fwd_cent_PR_True_pp.SetBinContent(i+1, PR_uncertTrue_pp);
		fwd_cent_PR_Err_pb.SetBinContent(i+1, PR_uncertErr_pb);
		fwd_cent_PR_Res_pb.SetBinContent(i+1, PR_uncertRes_pb);
		fwd_cent_PR_Bkg_pb.SetBinContent(i+1, PR_uncertBkg_pb);
		fwd_cent_PR_True_pb.SetBinContent(i+1, PR_uncertTrue_pb);
		fwd_cent_NP_Err_pp.SetBinContent(i+1, NP_uncertErr_pp);
		fwd_cent_NP_Res_pp.SetBinContent(i+1, NP_uncertRes_pp);
		fwd_cent_NP_Bkg_pp.SetBinContent(i+1, NP_uncertBkg_pp);
		fwd_cent_NP_True_pp.SetBinContent(i+1, NP_uncertTrue_pp);
		fwd_cent_NP_Err_pb.SetBinContent(i+1, NP_uncertErr_pb);
		fwd_cent_NP_Res_pb.SetBinContent(i+1, NP_uncertRes_pb);
		fwd_cent_NP_Bkg_pb.SetBinContent(i+1, NP_uncertBkg_pb);
		fwd_cent_NP_True_pb.SetBinContent(i+1, NP_uncertTrue_pb);

        fwd_int_PR.SetBinContent(i+1, PR_uncert_int);
        fwd_int_PR_pb.SetBinContent(i+1, PR_uncert_int_pb);
        fwd_int_NP.SetBinContent(i+1, NP_uncert_int);
        fwd_int_NP_pb.SetBinContent(i+1, NP_uncert_int_pb);

        fwd_int_PR_Err_pb.SetBinContent(i+1, PR_uncert_int_Err_pb);
        fwd_int_PR_Res_pb.SetBinContent(i+1, PR_uncert_int_Res_pb);
        fwd_int_PR_Bkg_pb.SetBinContent(i+1, PR_uncert_int_Bkg_pb);
        fwd_int_PR_True_pb.SetBinContent(i+1, PR_uncert_int_True_pb);
        fwd_int_NP_Err_pb.SetBinContent(i+1, NP_uncert_int_Err_pb);
        fwd_int_NP_Res_pb.SetBinContent(i+1, NP_uncert_int_Res_pb);
        fwd_int_NP_Bkg_pb.SetBinContent(i+1, NP_uncert_int_Bkg_pb);
        fwd_int_NP_True_pb.SetBinContent(i+1, NP_uncert_int_True_pb);
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

	fwd_cent_PR_Err_pp.Write();
	fwd_cent_PR_Res_pp.Write();
	fwd_cent_PR_Bkg_pp.Write();
	fwd_cent_PR_True_pp.Write();
	fwd_cent_PR_Err_pb.Write();
	fwd_cent_PR_Res_pb.Write();
	fwd_cent_PR_Bkg_pb.Write();
	fwd_cent_PR_True_pb.Write();
	fwd_cent_NP_Err_pp.Write();
	fwd_cent_NP_Res_pp.Write();
	fwd_cent_NP_Bkg_pp.Write();
	fwd_cent_NP_True_pp.Write();
	fwd_cent_NP_Err_pb.Write();
	fwd_cent_NP_Res_pb.Write();
	fwd_cent_NP_Bkg_pb.Write();
	fwd_cent_NP_True_pb.Write();

	mid_cent_PR_Err_pp.Write();
	mid_cent_PR_Res_pp.Write();
	mid_cent_PR_Bkg_pp.Write();
	mid_cent_PR_True_pp.Write();
	mid_cent_PR_Err_pb.Write();
	mid_cent_PR_Res_pb.Write();
	mid_cent_PR_Bkg_pb.Write();
	mid_cent_PR_True_pb.Write();
	mid_cent_NP_Err_pp.Write();
	mid_cent_NP_Res_pp.Write();
	mid_cent_NP_Bkg_pp.Write();
	mid_cent_NP_True_pp.Write();
	mid_cent_NP_Err_pb.Write();
	mid_cent_NP_Res_pb.Write();
	mid_cent_NP_Bkg_pb.Write();
	mid_cent_NP_True_pb.Write();
    
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

    out_int.cd();

    mid_int_PR.SetName("mid_int_PR");
    mid_int_PR_pb.SetName("mid_int_PR_pb");
    mid_int_NP.SetName("mid_int_NP");
    mid_int_NP_pb.SetName("mid_int_NP_pb");

    mid_int_PR_Err_pb.SetName("mid_int_Err_PR_pb");
    mid_int_PR_Res_pb.SetName("mid_int_Res_PR_pb");
    mid_int_PR_Bkg_pb.SetName("mid_int_Bkg_PR_pb");
    mid_int_PR_True_pb.SetName("mid_int_True_PR_pb");

    mid_int_NP_Err_pb.SetName("mid_int_Err_NP_pb");
    mid_int_NP_Res_pb.SetName("mid_int_Res_NP_pb");
    mid_int_NP_Bkg_pb.SetName("mid_int_Bkg_NP_pb");
    mid_int_NP_True_pb.SetName("mid_int_True_NP_pb");

    fwd_int_PR.SetName("fwd_int_PR");
    fwd_int_PR_pb.SetName("fwd_int_PR_pb");
    fwd_int_NP.SetName("fwd_int_NP");
    fwd_int_NP_pb.SetName("fwd_int_NP_pb");

    fwd_int_PR_Err_pb.SetName("fwd_int_Err_PR_pb");
    fwd_int_PR_Res_pb.SetName("fwd_int_Res_PR_pb");
    fwd_int_PR_Bkg_pb.SetName("fwd_int_Bkg_PR_pb");
    fwd_int_PR_True_pb.SetName("fwd_int_True_PR_pb");

    fwd_int_NP_Err_pb.SetName("fwd_int_Err_NP_pb");
    fwd_int_NP_Res_pb.SetName("fwd_int_Res_NP_pb");
    fwd_int_NP_Bkg_pb.SetName("fwd_int_Bkg_NP_pb");
    fwd_int_NP_True_pb.SetName("fwd_int_True_NP_pb");

    mid_int_PR.Write();
    mid_int_PR_pb.Write();
    mid_int_NP.Write();
    mid_int_NP_pb.Write();

    mid_int_PR_Err_pb.Write();
    mid_int_PR_Res_pb.Write();
    mid_int_PR_Bkg_pb.Write();
    mid_int_PR_True_pb.Write();

    mid_int_NP_Err_pb.Write();
    mid_int_NP_Res_pb.Write();
    mid_int_NP_Bkg_pb.Write();
    mid_int_NP_True_pb.Write();

    fwd_int_PR.Write();
    fwd_int_PR_pb.Write();
    fwd_int_NP.Write();
    fwd_int_NP_pb.Write();

    fwd_int_PR_Err_pb.Write();
    fwd_int_PR_Res_pb.Write();
    fwd_int_PR_Bkg_pb.Write();
    fwd_int_PR_True_pb.Write();

    fwd_int_NP_Err_pb.Write();
    fwd_int_NP_Res_pb.Write();
    fwd_int_NP_Bkg_pb.Write();
    fwd_int_NP_True_pb.Write();



    out_int.Close();

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
    n_NP = b_frac;
    n_PR = 1-b_frac;
}

//double compute_uncertainty(double pp_nomi, double pp_alpha, double pp_n, double pp_f, double pp_x, double pb_nomi, double pb_alpha, double pb_n, double pb_f, double pb_x)
double compute_uncertainty(double pp_syst, double pb_syst)
{
    // sqrt of (diff/n_PbPb)^2 + (diff/n_PP)^2
    return TMath::Sqrt((
        TMath::Power(pp_syst,2)+TMath::Power(pb_syst,2)
    ));
    //return TMath::Sqrt(
    //    TMath::Power(((pb_nomi-pb_alpha)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_alpha)/pp_nomi), 2)
    //    + TMath::Power(((pb_nomi-pb_n)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_n)/pp_nomi), 2)
    //    + TMath::Power(((pb_nomi-pb_f)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_f)/pp_nomi), 2)
    //    + TMath::Power(((pb_nomi-pb_x)/pb_nomi), 2) + TMath::Power(((pp_nomi-pp_x)/pp_nomi), 2)
    //    );
}

double compute_uncertainty(double nomi, double alpha, double n_value, double f_value, double x_value)
{
    return TMath::Sqrt((
        TMath::Power(((nomi-alpha)), 2)
        + TMath::Power(((nomi-n_value)), 2)
        + TMath::Power(((nomi-f_value)), 2)
        + TMath::Power(((nomi-x_value)), 2)
        )/4);
    //return TMath::Sqrt(
    //    TMath::Power(((nomi-alpha)/nomi), 2)
    //    + TMath::Power(((nomi-n_value)/nomi), 2)
    //    + TMath::Power(((nomi-f_value)/nomi), 2)
    //    + TMath::Power(((nomi-x_value)/nomi), 2)
    //    );
}
