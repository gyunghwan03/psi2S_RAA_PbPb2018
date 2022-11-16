void drawAcc(){
    gROOT->Macro("~/rootlogon.C");
    TFile *in1 = new TFile("acceptance_Prompt_GenOnly_wgt1_2021Psi2Sv2_20210601.root","READ");

    TH1F *hAcc1_pt1 = (TH1F*)in1->Get("hAccPt_2021_ally");
    TH1F *hAcc1_pt2 = (TH1F*)in1->Get("hAccPt_2021_Fory");

    hAcc1_pt1->SetMarkerStyle(20);
    hAcc1_pt1->SetMarkerColor(kRed+1);
    hAcc1_pt1->SetLineColor(kRed+1);
    hAcc1_pt2->SetMarkerStyle(24);
    hAcc1_pt2->SetMarkerColor(kBlue+1);
    hAcc1_pt2->SetLineColor(kBlue+1);

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();

    TLegend *lg1 = new TLegend(0.55,0.36,0.73,0.50);
    lg1->SetFillStyle(0);
    lg1->SetFillColor(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.035);
    lg1->AddEntry(hAcc1_pt1,"|y|: 0-2.4","lep");
    lg1->AddEntry(hAcc1_pt2,"|y|: 1.6-2.4","lep");

    TCanvas *c1 = new TCanvas("c1","",660,600);
    c1->cd();

    TH1F *hPad = new TH1F("hPad",";;",10,1.0,50.0);

    hPad->GetYaxis()->SetRangeUser(0.0,1.2);
    hPad->GetYaxis()->SetTitle("Acceptance");
    hPad->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPad->GetYaxis()->CenterTitle();
    hPad->GetXaxis()->CenterTitle();

    hAcc1_pt2->GetXaxis()->SetRangeUser(1.0,50.0);
    hAcc1_pt2->GetYaxis()->SetRangeUser(0.0,1.2);
    hAcc1_pt2->GetYaxis()->SetTitle("Acceptance");
    hAcc1_pt2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hAcc1_pt2->GetYaxis()->CenterTitle();
    hAcc1_pt2->GetXaxis()->CenterTitle();

    hPad->Draw();
    hAcc1_pt2->Draw("E same");
    hAcc1_pt1->Draw("E same");

    lg1->Draw("same");
    //lt1->DrawLatex(0.20,0.87,"Trigger Matching");
    lt1->DrawLatex(0.20,0.87,"CMS Simulation");
    lt1->DrawLatex(0.20,0.80,"Prompt #psi(2S)");
    lt1->DrawLatex(0.60,0.87,"PbPb #sqrt{s_{NN}} = 5.02 TeV");

    c1->SaveAs("acc_total_psi2s_pt.png");
    c1->SaveAs("acc_total_psi2s_pt.pdf");

    TFile *out = new TFile("acc_psi2s_pbpb_ptwgt1.root","RECREATE");
    out->cd();
    
    hAcc1_pt1->SetName("hAcc_pt1");
    hAcc1_pt2->SetName("hAcc_pt2");

    hAcc1_pt1->Write();
    hAcc1_pt2->Write();
    out->Close();


}
