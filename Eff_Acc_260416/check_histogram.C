#include <iostream>
#include <TFile.h>
#include <TH1D.h>

void check_histogram() {
    // Check if the root file exists
    TFile *f = TFile::Open("roots_2S_pp/ctau3D_cut_ptBin_PRMC_y0.0-2.4.root");
    if (!f || f->IsZombie()) {
        cout << "Cannot open file!" << endl;
        return;
    }
    
    cout << "File opened successfully" << endl;
    cout << "\nChecking mid rapidity (0-1.6) histograms:" << endl;
    
    // Check each pt bin for mid rapidity
    for (int ipt = 0; ipt < 7; ipt++) {
        TH1D* h_pr = (TH1D*)f->Get(Form("h_decay_pr_pt_mid_%d", ipt));
        TH1D* h_np = (TH1D*)f->Get(Form("h_decay_np_pt_mid_%d", ipt));
        
        if (h_pr && h_np) {
            double entries_pr = h_pr->GetEntries();
            double entries_np = h_np->GetEntries();
            double integral_pr = h_pr->Integral();
            double integral_np = h_np->Integral();
            
            cout << Form("  ipt=%d: PR entries=%.0f, integral=%.0f | NP entries=%.0f, integral=%.0f", 
                         ipt, entries_pr, integral_pr, entries_np, integral_np) << endl;
        } else {
            cout << Form("  ipt=%d: Histogram not found!", ipt) << endl;
        }
    }
    
    f->Close();
}
