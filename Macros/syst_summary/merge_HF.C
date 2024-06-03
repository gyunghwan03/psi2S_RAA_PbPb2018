#include "TMath.h"

using TMath::Power; using TMath::Sqrt;

void merge_HF()
{
	auto HF_up_pt = TFile::Open("syst_roots/syst_pt_HF_Up.root");
	auto HF_down_pt = TFile::Open("syst_roots/syst_pt_HF_Down.root");

	auto HF_up_cent = TFile::Open("syst_roots/syst_cent_HF_Up.root");
	auto HF_down_cent = TFile::Open("syst_roots/syst_cent_HF_Down.root");

	auto h_up_mid_pt_PR = (TH1D*) HF_up_pt->Get("mid_PR");
	auto h_up_mid_pt_NP = (TH1D*) HF_up_pt->Get("mid_NP");
	auto h_down_mid_pt_PR = (TH1D*) HF_down_pt->Get("mid_PR");
	auto h_down_mid_pt_NP = (TH1D*) HF_down_pt->Get("mid_NP");
	auto h_up_fwd_pt_PR = (TH1D*) HF_up_pt->Get("fwd_PR");
	auto h_up_fwd_pt_NP = (TH1D*) HF_up_pt->Get("fwd_NP");
	auto h_down_fwd_pt_PR = (TH1D*) HF_down_pt->Get("fwd_PR");
	auto h_down_fwd_pt_NP = (TH1D*) HF_down_pt->Get("fwd_NP");

	auto h_up_mid_cent_PR = (TH1D*) HF_up_cent->Get("mid_PR");
	auto h_up_mid_cent_NP = (TH1D*) HF_up_cent->Get("mid_NP");
	auto h_down_mid_cent_PR = (TH1D*) HF_down_cent->Get("mid_PR");
	auto h_down_mid_cent_NP = (TH1D*) HF_down_cent->Get("mid_NP");
	auto h_up_fwd_cent_PR = (TH1D*) HF_up_cent->Get("fwd_PR");
	auto h_up_fwd_cent_NP = (TH1D*) HF_up_cent->Get("fwd_NP");
	auto h_down_fwd_cent_PR = (TH1D*) HF_down_cent->Get("fwd_PR");
	auto h_down_fwd_cent_NP = (TH1D*) HF_down_cent->Get("fwd_NP");

	const int NBINS_mid_pt = 6;
    //double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 30, 50};
    double edges_mid_pt[NBINS_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 50};
    auto h_hf_mid_pt_PR = new TH1D("mid_PR", "mid_PR", NBINS_mid_pt, edges_mid_pt);
    auto h_hf_mid_pt_NP = new TH1D("mid_NP", "mid_NP", NBINS_mid_pt, edges_mid_pt);

	for (int hist_idx = 1; hist_idx < NBINS_mid_pt+1; hist_idx++) {
        double hf_mid_PR = 0;
		double up_mid_PR = h_up_mid_pt_PR->GetBinContent(hist_idx);
		double down_mid_PR = h_down_mid_pt_PR->GetBinContent(hist_idx);
		
		hf_mid_PR = Sqrt(Power(up_mid_PR,2) + Power(down_mid_PR,2));
		h_hf_mid_pt_PR -> SetBinContent(hist_idx,hf_mid_PR);

        double hf_mid_NP = 0;
		double up_mid_NP = h_up_mid_pt_NP->GetBinContent(hist_idx);
		double down_mid_NP = h_down_mid_pt_NP->GetBinContent(hist_idx);
		
		hf_mid_NP = Sqrt(Power(up_mid_NP,2) + Power(down_mid_NP,2));
		h_hf_mid_pt_NP -> SetBinContent(hist_idx,hf_mid_NP);
	}

	const int NBINS_fwd_pt = 4;
    double edges_fwd_pt[NBINS_fwd_pt+1] = {3.5, 6.5, 9, 12, 50};
    auto h_hf_fwd_pt_PR = new TH1D("fwd_PR", "fwd_PR", NBINS_fwd_pt, edges_fwd_pt);
    auto h_hf_fwd_pt_NP = new TH1D("fwd_NP", "fwd_NP", NBINS_fwd_pt, edges_fwd_pt);

	for (int hist_idx = 1; hist_idx < NBINS_fwd_pt+1; hist_idx++) {
        double hf_fwd_PR = 0;
		double up_fwd_PR = h_up_fwd_pt_PR->GetBinContent(hist_idx);
		double down_fwd_PR = h_down_fwd_pt_PR->GetBinContent(hist_idx);
		
		hf_fwd_PR = Sqrt(Power(up_fwd_PR,2) + Power(down_fwd_PR,2));
		h_hf_fwd_pt_PR -> SetBinContent(hist_idx,hf_fwd_PR);

        double hf_fwd_NP = 0;
		double up_fwd_NP = h_up_fwd_pt_NP->GetBinContent(hist_idx);
		double down_fwd_NP = h_down_fwd_pt_NP->GetBinContent(hist_idx);
		
		hf_fwd_NP = Sqrt(Power(up_fwd_NP,2) + Power(down_fwd_NP,2));
		h_hf_fwd_pt_NP -> SetBinContent(hist_idx,hf_fwd_NP);
	}
	TFile out_pt("syst_roots/syst_pt_HF.root", "recreate");
	out_pt.cd();
	h_hf_mid_pt_PR->Write();
	h_hf_fwd_pt_PR->Write();
	h_hf_mid_pt_NP->Write();
	h_hf_fwd_pt_NP->Write();

	const int NBINS_mid_cent = 6;
    double edges_mid_cent[NBINS_mid_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    auto h_hf_mid_cent_PR = new TH1D("mid_PR", "mid_PR", NBINS_mid_cent, edges_mid_cent);
    auto h_hf_mid_cent_NP = new TH1D("mid_NP", "mid_NP", NBINS_mid_cent, edges_mid_cent);

	for (int hist_idx = 1; hist_idx < NBINS_mid_cent+1; hist_idx++) {
        double hf_mid_PR = 0;
		double up_mid_PR = h_up_mid_cent_PR->GetBinContent(hist_idx);
		double down_mid_PR = h_down_mid_cent_PR->GetBinContent(hist_idx);
		
		hf_mid_PR = Sqrt(Power(up_mid_PR,2) + Power(down_mid_PR,2));
		h_hf_mid_cent_PR -> SetBinContent(hist_idx,hf_mid_PR);

        double hf_mid_NP = 0;
		double up_mid_NP = h_up_mid_cent_NP->GetBinContent(hist_idx);
		double down_mid_NP = h_down_mid_cent_NP->GetBinContent(hist_idx);
		
		hf_mid_NP = Sqrt(Power(up_mid_NP,2) + Power(down_mid_NP,2));
		h_hf_mid_cent_NP -> SetBinContent(hist_idx,hf_mid_NP);
	}

	const int NBINS_fwd_cent = 6;
    double edges_fwd_cent[NBINS_fwd_cent+1] = {0, 10, 20, 30, 40, 50, 90};
    auto h_hf_fwd_cent_PR = new TH1D("fwd_PR", "fwd_PR", NBINS_fwd_cent, edges_fwd_cent);
    auto h_hf_fwd_cent_NP = new TH1D("fwd_NP", "fwd_NP", NBINS_fwd_cent, edges_fwd_cent);

	for (int hist_idx = 1; hist_idx < NBINS_fwd_cent+1; hist_idx++) {
        double hf_fwd_PR = 0;
		double up_fwd_PR = h_up_fwd_cent_PR->GetBinContent(hist_idx);
		double down_fwd_PR = h_down_fwd_cent_PR->GetBinContent(hist_idx);
		
		hf_fwd_PR = Sqrt(Power(up_fwd_PR,2) + Power(down_fwd_PR,2));
		h_hf_fwd_cent_PR -> SetBinContent(hist_idx,hf_fwd_PR);

        double hf_fwd_NP = 0;
		double up_fwd_NP = h_up_fwd_cent_NP->GetBinContent(hist_idx);
		double down_fwd_NP = h_down_fwd_cent_NP->GetBinContent(hist_idx);
		
		hf_fwd_NP = Sqrt(Power(up_fwd_NP,2) + Power(down_fwd_NP,2));
		h_hf_fwd_cent_NP -> SetBinContent(hist_idx,hf_fwd_NP);
	}

	TFile out_cent("syst_roots/syst_cent_HF.root", "recreate");
	out_cent.cd();
	h_hf_mid_cent_PR->Write();
	h_hf_fwd_cent_PR->Write();
	h_hf_mid_cent_NP->Write();
	h_hf_fwd_cent_NP->Write();
}
