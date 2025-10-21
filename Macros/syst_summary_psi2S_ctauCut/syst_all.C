// syst_all.C
// - HFup/HFdown: apply variation on Pb side only (pp uses nominal)
// - sigPAR: combine {alpha, f, n, x} by quadrature (bin-by-bin)
// - Also save per-subvariation histograms for sigPAR: alpha, f, n, x
// - Suppress TH1 auto-registration & name collisions; free temps to avoid leaks.

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"   // TH1::AddDirectory
#include "TH1D.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"

// Provide these from your env (input_list.h):
//   pp_mid_pt, pb_mid_pt, pp_fwd_pt, pb_fwd_pt,
//   pp_mid_cent, pb_mid_cent, pp_fwd_cent, pb_fwd_cent
#include "input_list.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// ------------ small utils ------------
static string sanitize(const string& s) {
    string t = s;
    for (auto& c : t) if (c=='/' || c==' ' || c=='\\' || c==':') c = '_';
    return t;
}
static string sigpar_label_from_source(const string& source) {
    // Turn "sigPAR/sigPAR_alpha" -> "sigPAR_alpha", etc.
    const string key = "SigPAR_";
    size_t pos = source.find(key);
    if (pos != string::npos) {
        string sub = source.substr(pos + key.size()); // e.g. "alpha"
        return "sigPAR_" + sub;
    }
    // Fallback: sanitize
    return sanitize(source);
}

// Read yields (PR/NP) from directory layout:
//   file_dir + "PRMC/" + type[i]   and   file_dir + "NPMC/" + type[i]
static void get_n_jpsi_all(const string& file_dir,
                           const vector<string>& type,
                           int i,
                           double &n_PR,
                           double &n_NP)
{
    n_PR = 0; n_NP = 0;

    const string fPR = file_dir + "PRMC/" + type[i];
    const string fNP = file_dir + "NPMC/" + type[i];

    TFile *PR = TFile::Open(fPR.c_str());
    if (!PR || PR->IsZombie()) {
        cerr << "[WARN] Fail open PR: " << fPR << endl;
    } else {
        if (auto* h = dynamic_cast<TH1D*>(PR->Get("fitResults"))) n_PR = h->GetBinContent(1);
        else cerr << "[WARN] Missing hist 'fitResults' in " << fPR << endl;
        PR->Close();
    }

    TFile *NP = TFile::Open(fNP.c_str());
    if (!NP || NP->IsZombie()) {
        cerr << "[WARN] Fail open NP: " << fNP << endl;
    } else {
        if (auto* h = dynamic_cast<TH1D*>(NP->Get("fitResults"))) n_NP = h->GetBinContent(1);
        else cerr << "[WARN] Missing hist 'fitResults' in " << fNP << endl;
        NP->Close();
    }
}

static inline double frac_unc(double nomi, double syst)
{
    if (nomi == 0) return 0;
    return TMath::Abs(nomi - syst) / TMath::Abs(nomi);
}
static inline double RAA_like_unc(double pp_nomi, double pp_syst,
                                  double pb_nomi, double pb_syst)
{
    if (pp_nomi == 0 || pb_nomi == 0) return 0;
    const double a = (pb_nomi - pb_syst) / pb_nomi;
    const double b = (pp_nomi - pp_syst) / pp_nomi;
    return TMath::Sqrt(a*a + b*b);
}

// ------------ container for 6 hists ------------
struct SixHists {
    TH1D *PR, *PR_pp, *PR_pb;
    TH1D *NP, *NP_pp, *NP_pb;
};

static TH1D* mk1(const string& name, int nbins, const double* edges) {
    TH1D* h = new TH1D(name.c_str(), name.c_str(), nbins, edges);
    h->SetDirectory(nullptr);   // prevent auto registration to gDirectory
    return h;
}
static SixHists make_six(const string& tag, int nbins, const double* edges)
{
    SixHists S;
    S.PR    = mk1(tag + "_PR",    nbins, edges);
    S.PR_pp = mk1(tag + "_PR_pp", nbins, edges);
    S.PR_pb = mk1(tag + "_PR_pb", nbins, edges);
    S.NP    = mk1(tag + "_NP",    nbins, edges);
    S.NP_pp = mk1(tag + "_NP_pp", nbins, edges);
    S.NP_pb = mk1(tag + "_NP_pb", nbins, edges);
    return S;
}
static void free_six(SixHists& S){
    delete S.PR; delete S.PR_pp; delete S.PR_pb;
    delete S.NP; delete S.NP_pp; delete S.NP_pb;
}

static void quad_add(SixHists& dst, const SixHists& src, bool first)
{
    auto qadd = [first](TH1D* d, TH1D* s){
        const int nb = d->GetNbinsX();
        for (int i=1;i<=nb;++i){
            const double v2 = s->GetBinContent(i) * s->GetBinContent(i);
            if (first) d->SetBinContent(i, v2);
            else       d->SetBinContent(i, d->GetBinContent(i) + v2);
        }
    };
    qadd(dst.PR,    src.PR);
    qadd(dst.PR_pp, src.PR_pp);
    qadd(dst.PR_pb, src.PR_pb);
    qadd(dst.NP,    src.NP);
    qadd(dst.NP_pp, src.NP_pp);
    qadd(dst.NP_pb, src.NP_pb);
}
static void sqrt_finalize(SixHists& H)
{
    auto sq = [](TH1D* h){
        const int nb = h->GetNbinsX();
        for (int i=1;i<=nb;++i){
            const double x = h->GetBinContent(i);
            h->SetBinContent(i, TMath::Sqrt(std::max(0.0, x)));
        }
    };
    sq(H.PR); sq(H.PR_pp); sq(H.PR_pb);
    sq(H.NP); sq(H.NP_pp); sq(H.NP_pb);
}

// ------------ writer helper ------------
static void write_outputs(const string& out_pt_name,
                          const string& out_cent_name,
                          SixHists& mid_pt, SixHists& fwd_pt,
                          SixHists& mid_cent, SixHists& fwd_cent)
{
    // pt file
    {
        TFile f(out_pt_name.c_str(), "recreate");
        mid_pt.PR->SetName("mid_PR");         mid_pt.PR->Write();
        mid_pt.PR_pp->SetName("mid_PR_pp");   mid_pt.PR_pp->Write();
        mid_pt.PR_pb->SetName("mid_PR_pb");   mid_pt.PR_pb->Write();
        mid_pt.NP->SetName("mid_NP");         mid_pt.NP->Write();
        mid_pt.NP_pp->SetName("mid_NP_pp");   mid_pt.NP_pp->Write();
        mid_pt.NP_pb->SetName("mid_NP_pb");   mid_pt.NP_pb->Write();

        fwd_pt.PR->SetName("fwd_PR");         fwd_pt.PR->Write();
        fwd_pt.PR_pp->SetName("fwd_PR_pp");   fwd_pt.PR_pp->Write();
        fwd_pt.PR_pb->SetName("fwd_PR_pb");   fwd_pt.PR_pb->Write();
        fwd_pt.NP->SetName("fwd_NP");         fwd_pt.NP->Write();
        fwd_pt.NP_pp->SetName("fwd_NP_pp");   fwd_pt.NP_pp->Write();
        fwd_pt.NP_pb->SetName("fwd_NP_pb");   fwd_pt.NP_pb->Write();
    }
    // cent file
    {
        TFile f(out_cent_name.c_str(), "recreate");
        mid_cent.PR->SetName("mid_PR");       mid_cent.PR->Write();
        mid_cent.PR_pp->SetName("mid_PR_pp"); mid_cent.PR_pp->Write();
        mid_cent.PR_pb->SetName("mid_PR_pb"); mid_cent.PR_pb->Write();
        mid_cent.NP->SetName("mid_NP");       mid_cent.NP->Write();
        mid_cent.NP_pp->SetName("mid_NP_pp"); mid_cent.NP_pp->Write();
        mid_cent.NP_pb->SetName("mid_NP_pb"); mid_cent.NP_pb->Write();

        fwd_cent.PR->SetName("fwd_PR");       fwd_cent.PR->Write();
        fwd_cent.PR_pp->SetName("fwd_PR_pp"); fwd_cent.PR_pp->Write();
        fwd_cent.PR_pb->SetName("fwd_PR_pb"); fwd_cent.PR_pb->Write();
        fwd_cent.NP->SetName("fwd_NP");       fwd_cent.NP->Write();
        fwd_cent.NP_pp->SetName("fwd_NP_pp"); fwd_cent.NP_pp->Write();
        fwd_cent.NP_pb->SetName("fwd_NP_pb"); fwd_cent.NP_pb->Write();
    }
}

// ------------ core runner ------------
static void run_one_type_or_combined(const string& syst_type,
                                     const vector<string>& subtypes_to_combine,
                                     bool save_subcomponents = false)
{
    gSystem->mkdir("syst_roots", true);

    const string nominal_path_pp = "../psi2S_L_cut_250427/roots_2S_pp/";
    const string nominal_path_pb = "../psi2S_L_cut_250427/roots_2S_Pb/";
    const string syst_base_pp    = "../psi2S_L_cut_250427/systematic_pp/";
    const string syst_base_pb    = "../psi2S_L_cut_250427/systematic_Pb/";

    const bool do_combine = !subtypes_to_combine.empty();

    // ---- binning ----
    const int NB_mid_pt = 6; double E_mid_pt[NB_mid_pt+1] = {6.5, 9, 12, 15, 20, 25, 40};
    const int NB_fwd_pt = 4; double E_fwd_pt[NB_fwd_pt+1] = {3.5, 6.5, 9, 12, 40};
    const int NB_mid_cent = 6; double E_mid_cent[NB_mid_cent+1] = {0,10,20,30,40,50,90};
    const int NB_fwd_cent = 4; double E_fwd_cent[NB_fwd_cent+1] = {0,10,30,50,90};

    // ---- aggregated containers ----
    SixHists mid_pt  = make_six("mid_pt",  NB_mid_pt,  E_mid_pt);
    SixHists fwd_pt  = make_six("fwd_pt",  NB_fwd_pt,  E_fwd_pt);
    SixHists mid_cent= make_six("mid_cent",NB_mid_cent,E_mid_cent);
    SixHists fwd_cent= make_six("fwd_cent",NB_fwd_cent,E_fwd_cent);

    // ---- lambdas to compute one source ----
    auto fill_one_source_midpt = [&](const string& source)->SixHists{
        const string tag = "tmp_mid_pt_" + sanitize(source);
        SixHists H = make_six(tag, NB_mid_pt, E_mid_pt);

        const bool isHF = (source=="HFup" || source=="HFdown");
        const string pp_syst_dir = isHF ? nominal_path_pp : (syst_base_pp + source + "/roots/");

        for (int i=0;i<(int)pp_mid_pt.size();++i){
            double ppPRn, ppNPn, ppPRs, ppNPs;
            double pbPRn, pbNPn, pbPRs, pbNPs;

            get_n_jpsi_all(nominal_path_pp, pp_mid_pt, i, ppPRn, ppNPn);
            get_n_jpsi_all(pp_syst_dir,     pp_mid_pt, i, ppPRs, ppNPs);

            get_n_jpsi_all(nominal_path_pb, pb_mid_pt, i, pbPRn, pbNPn);
            get_n_jpsi_all(syst_base_pb + source + "/roots/", pb_mid_pt, i, pbPRs, pbNPs);

            const double PR_unc    = RAA_like_unc(ppPRn, ppPRs, pbPRn, pbPRs);
            const double PR_unc_pp = frac_unc(ppPRn,     ppPRs);
            const double PR_unc_pb = frac_unc(pbPRn,     pbPRs);

            const double NP_unc    = RAA_like_unc(ppNPn, ppNPs, pbNPn, pbNPs);
            const double NP_unc_pp = frac_unc(ppNPn,     ppNPs);
            const double NP_unc_pb = frac_unc(pbNPn,     pbNPs);

            H.PR->SetBinContent(i+1, PR_unc);
            H.PR_pp->SetBinContent(i+1, PR_unc_pp);
            H.PR_pb->SetBinContent(i+1, PR_unc_pb);
            H.NP->SetBinContent(i+1, NP_unc);
            H.NP_pp->SetBinContent(i+1, NP_unc_pp);
            H.NP_pb->SetBinContent(i+1, NP_unc_pb);
        }
        return H;
    };
    auto fill_one_source_fwdpt = [&](const string& source)->SixHists{
        const string tag = "tmp_fwd_pt_" + sanitize(source);
        SixHists H = make_six(tag, NB_fwd_pt, E_fwd_pt);

        const bool isHF = (source=="HFup" || source=="HFdown");
        const string pp_syst_dir = isHF ? nominal_path_pp : (syst_base_pp + source + "/roots/");

        for (int i=0;i<(int)pp_fwd_pt.size();++i){
            double ppPRn, ppNPn, ppPRs, ppNPs;
            double pbPRn, pbNPn, pbPRs, pbNPs;

            get_n_jpsi_all(nominal_path_pp, pp_fwd_pt, i, ppPRn, ppNPn);
            get_n_jpsi_all(pp_syst_dir,     pp_fwd_pt, i, ppPRs, ppNPs);

            get_n_jpsi_all(nominal_path_pb, pb_fwd_pt, i, pbPRn, pbNPn);
            get_n_jpsi_all(syst_base_pb + source + "/roots/", pb_fwd_pt, i, pbPRs, pbNPs);

            const double PR_unc    = RAA_like_unc(ppPRn, ppPRs, pbPRn, pbPRs);
            const double PR_unc_pp = frac_unc(ppPRn,     ppPRs);
            const double PR_unc_pb = frac_unc(pbPRn,     pbPRs);

            const double NP_unc    = RAA_like_unc(ppNPn, ppNPs, pbNPn, pbNPs);
            const double NP_unc_pp = frac_unc(ppNPn,     ppNPs);
            const double NP_unc_pb = frac_unc(pbNPn,     pbNPs);

            H.PR->SetBinContent(i+1, PR_unc);
            H.PR_pp->SetBinContent(i+1, PR_unc_pp);
            H.PR_pb->SetBinContent(i+1, PR_unc_pb);
            H.NP->SetBinContent(i+1, NP_unc);
            H.NP_pp->SetBinContent(i+1, NP_unc_pp);
            H.NP_pb->SetBinContent(i+1, NP_unc_pb);
        }
        return H;
    };
    auto fill_one_source_midcent = [&](const string& source)->SixHists{
        const string tag = "tmp_mid_cent_" + sanitize(source);
        SixHists H = make_six(tag, NB_mid_cent, E_mid_cent);

        const bool isHF = (source=="HFup" || source=="HFdown");
        const string pp_syst_dir = isHF ? nominal_path_pp : (syst_base_pp + source + "/roots/");

        for (int i=0;i<(int)pb_mid_cent.size();++i){
            double ppPRn, ppNPn, ppPRs, ppNPs;
            double pbPRn, pbNPn, pbPRs, pbNPs;

            // pp has single cent bin (index 0)
            get_n_jpsi_all(nominal_path_pp, pp_mid_cent, 0, ppPRn, ppNPn);
            get_n_jpsi_all(pp_syst_dir,     pp_mid_cent, 0, ppPRs, ppNPs);

            get_n_jpsi_all(nominal_path_pb, pb_mid_cent, i, pbPRn, pbNPn);
            get_n_jpsi_all(syst_base_pb + source + "/roots/", pb_mid_cent, i, pbPRs, pbNPs);

            const double PR_unc    = RAA_like_unc(ppPRn, ppPRs, pbPRn, pbPRs);
            const double PR_unc_pp = frac_unc(ppPRn,     ppPRs);
            const double PR_unc_pb = frac_unc(pbPRn,     pbPRs);

            const double NP_unc    = RAA_like_unc(ppNPn, ppNPs, pbNPn, pbNPs);
            const double NP_unc_pp = frac_unc(ppNPn,     ppNPs);
            const double NP_unc_pb = frac_unc(pbNPn,     pbNPs);

            H.PR->SetBinContent(i+1, PR_unc);
            H.PR_pp->SetBinContent(i+1, NP_unc_pp); // keep PR_pp/NP_pp consistent
            H.PR_pb->SetBinContent(i+1, PR_unc_pb);
            H.NP->SetBinContent(i+1, NP_unc);
            H.NP_pp->SetBinContent(i+1, NP_unc_pp);
            H.NP_pb->SetBinContent(i+1, NP_unc_pb);
        }
        return H;
    };
    auto fill_one_source_fwdcent = [&](const string& source)->SixHists{
        const string tag = "tmp_fwd_cent_" + sanitize(source);
        SixHists H = make_six(tag, NB_fwd_cent, E_fwd_cent);

        const bool isHF = (source=="HFup" || source=="HFdown");
        const string pp_syst_dir = isHF ? nominal_path_pp : (syst_base_pp + source + "/roots/");

        for (int i=0;i<(int)pb_fwd_cent.size();++i){
            double ppPRn, ppNPn, ppPRs, ppNPs;
            double pbPRn, pbNPn, pbPRs, pbNPs;

            // pp has single cent bin (index 0)
            get_n_jpsi_all(nominal_path_pp, pp_fwd_cent, 0, ppPRn, ppNPn);
            get_n_jpsi_all(pp_syst_dir,     pp_fwd_cent, 0, ppPRs, ppNPs);

            get_n_jpsi_all(nominal_path_pb, pb_fwd_cent, i, pbPRn, pbNPn);
            get_n_jpsi_all(syst_base_pb + source + "/roots/", pb_fwd_cent, i, pbPRs, pbNPs);

            const double PR_unc    = RAA_like_unc(ppPRn, ppPRs, pbPRn, pbPRs);
            const double PR_unc_pp = frac_unc(ppPRn,     ppPRs);
            const double PR_unc_pb = frac_unc(pbPRn,     pbPRs);

            const double NP_unc    = RAA_like_unc(ppNPn, ppNPs, pbNPn, pbNPs);
            const double NP_unc_pp = frac_unc(ppNPn,     ppNPs);
            const double NP_unc_pb = frac_unc(pbNPn,     pbNPs);

            H.PR->SetBinContent(i+1, PR_unc);
            H.PR_pp->SetBinContent(i+1, PR_unc_pp);
            H.PR_pb->SetBinContent(i+1, PR_unc_pb);
            H.NP->SetBinContent(i+1, NP_unc);
            H.NP_pp->SetBinContent(i+1, NP_unc_pp);
            H.NP_pb->SetBinContent(i+1, NP_unc_pb);
        }
        return H;
    };

    // ---- compute ----
    if (!do_combine) {
        // single source mode
        auto H1 = fill_one_source_midpt (syst_type);
        auto H2 = fill_one_source_fwdpt (syst_type);
        auto H3 = fill_one_source_midcent(syst_type);
        auto H4 = fill_one_source_fwdcent(syst_type);

        quad_add(mid_pt,  H1, /*first=*/true);
        quad_add(fwd_pt,  H2, /*first=*/true);
        quad_add(mid_cent,H3, /*first=*/true);
        quad_add(fwd_cent,H4, /*first=*/true);

        // free temps
        free_six(H1); free_six(H2); free_six(H3); free_six(H4);

        // finalize aggregated
        sqrt_finalize(mid_pt);
        sqrt_finalize(fwd_pt);
        sqrt_finalize(mid_cent);
        sqrt_finalize(fwd_cent);

        // write aggregated
        const string out_pt_name   = "./syst_roots/syst_pt_"   + syst_type + ".root";
        const string out_cent_name = "./syst_roots/syst_cent_" + syst_type + ".root";
        write_outputs(out_pt_name, out_cent_name, mid_pt, fwd_pt, mid_cent, fwd_cent);
        cout << "[DONE] " << syst_type << " -> " << out_pt_name << " , " << out_cent_name << endl;
    } else {
        // combined mode (e.g., sigPAR)
        vector<SixHists> subs_mid_pt, subs_fwd_pt, subs_mid_cent, subs_fwd_cent;
        bool first = true;

        for (const auto& sub : subtypes_to_combine) {
            auto H1 = fill_one_source_midpt (sub);
            auto H2 = fill_one_source_fwdpt (sub);
            auto H3 = fill_one_source_midcent(sub);
            auto H4 = fill_one_source_fwdcent(sub);

            quad_add(mid_pt,  H1, first);
            quad_add(fwd_pt,  H2, first);
            quad_add(mid_cent,H3, first);
            quad_add(fwd_cent,H4, first);
            first = false;

            if (save_subcomponents) {
                subs_mid_pt.push_back(H1);
                subs_fwd_pt.push_back(H2);
                subs_mid_cent.push_back(H3);
                subs_fwd_cent.push_back(H4);
            } else {
                // not saving, free immediately
                free_six(H1); free_six(H2); free_six(H3); free_six(H4);
            }
        }

        // finalize aggregated
        sqrt_finalize(mid_pt);
        sqrt_finalize(fwd_pt);
        sqrt_finalize(mid_cent);
        sqrt_finalize(fwd_cent);

        // write aggregated (e.g., sigPAR combined)
        const string out_pt_name   = "./syst_roots/syst_pt_"   + syst_type + ".root";
        const string out_cent_name = "./syst_roots/syst_cent_" + syst_type + ".root";
        write_outputs(out_pt_name, out_cent_name, mid_pt, fwd_pt, mid_cent, fwd_cent);
        cout << "[DONE] " << syst_type << " (combined) -> " << out_pt_name << " , " << out_cent_name << endl;

        // write and free subcomponents if requested
        if (save_subcomponents) {
            for (size_t i=0; i<subtypes_to_combine.size(); ++i) {
                const string label = sigpar_label_from_source(subtypes_to_combine[i]); // e.g. "sigPAR_alpha"
                const string out_pt_sub   = "./syst_roots/syst_pt_"   + label + ".root";
                const string out_cent_sub = "./syst_roots/syst_cent_" + label + ".root";
                write_outputs(out_pt_sub, out_cent_sub,
                              subs_mid_pt[i], subs_fwd_pt[i], subs_mid_cent[i], subs_fwd_cent[i]);
                cout << "  └─ saved subcomponent: " << label
                     << " -> " << out_pt_sub << " , " << out_cent_sub << endl;

                // free each subcomponent's hists
                free_six(subs_mid_pt[i]);
                free_six(subs_fwd_pt[i]);
                free_six(subs_mid_cent[i]);
                free_six(subs_fwd_cent[i]);
            }
        }
    }

    // free aggregated hists (avoid growth across calls)
    free_six(mid_pt);
    free_six(fwd_pt);
    free_six(mid_cent);
    free_six(fwd_cent);
}

// ------------ CSV parser & entry ------------
static vector<string> parse_csv_types(const string& csv)
{
    vector<string> out;
    std::istringstream iss(csv);
    string item;
    while (std::getline(iss, item, ',')) {
        size_t b = item.find_first_not_of(" \t");
        size_t e = item.find_last_not_of(" \t");
        if (b != string::npos) out.emplace_back(item.substr(b, e - b + 1));
    }
    return out;
}

void syst_all(const char* types_csv = "bkgPDF,HFdown,HFup,sigPDF,sigPAR")
{
    // turn off TH1 auto-registration to avoid name collisions
    const bool prevAddDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    const auto types = parse_csv_types(types_csv ? string(types_csv) : string());
    if (types.empty()) { cerr << "[ERROR] No systematic types provided.\n"; TH1::AddDirectory(prevAddDir); return; }

    for (const auto& t : types) {
        if (t == "sigPAR") {
            // Combine sub-variations in quadrature:
            // expected folder layout under systematic_{pp,Pb}/:
            //   sigPAR/sigPAR_{alpha,f,n,x}
            const vector<string> sub = {
                "sigPAR/sigPAR_alpha",
                "sigPAR/sigPAR_f",
                "sigPAR/sigPAR_n",
                "sigPAR/sigPAR_x"
            };
            cout << "\n========== Running combined sigPAR from {alpha,f,n,x} ==========\n";
            // save_subcomponents = true to also dump each variation's hists
            run_one_type_or_combined("sigPAR", sub, /*save_subcomponents=*/true);
        } else {
            cout << "\n========== Running syst type: " << t << " ==========\n";
            run_one_type_or_combined(t, /*subtypes_to_combine=*/{}, /*save_subcomponents=*/false);
        }
    }

    // restore previous AddDirectory status
    TH1::AddDirectory(prevAddDir);

    cout << "\nAll requested systematics finished.\n";
}
