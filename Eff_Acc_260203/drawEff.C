void drawEff(){
    gROOT->Macro("~/rootlogon.C");
    gStyle->SetOptFit(0);
    //TFile *in = new TFile("mc_eff_vs_pt_cent_0_to_200_rap_prompt_pbpb_Jpsi_PtW0_tnp0_total.root","READ");
    TFile *in1 = new TFile("mc_eff_vs_pt_cent_20_to_120_rap_prompt_pbpb_psi2s_PtW1_tnp1_20220125.root","READ");

    TH1F *hy1_eff = (TH1F*)in1->Get("hy_reco");
    TH1F *hy1_gen = (TH1F*)in1->Get("hy_gen");

    hy1_eff->Divide(hy1_gen);

    hy1_eff->SetMarkerStyle(20);
    hy1_eff->SetMarkerColor(kBlue+1);
    hy1_eff->SetLineColor(kBlue+1);

    TH1F *hEff1_pt1 = (TH1F*)in1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_20_to_120_absy0_2p4");
    TH1F *hEff1_pt2 = (TH1F*)in1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_20_to_120_absy0_1p2");
    TH1F *hEff1_pt3 = (TH1F*)in1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_20_to_120_absy1p2_1p6");
    TH1F *hEff1_pt4 = (TH1F*)in1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_20_to_120_absy1p6_1p8");
    TH1F *hEff1_pt5 = (TH1F*)in1->Get("mc_eff_vs_pt_TnP1_PtW1_cent_20_to_120_absy1p8_2p4");
    TH1F *hEff1_pt6 = (TH1F*)in1->Get("hpt_reco_6");


    /*
    TH1F *hEff1_pt1 = (TH1F*)in1->Get("hpt_reco_1");
    TH1F *hEff1_pt2 = (TH1F*)in1->Get("hpt_reco_2");
    TH1F *hEff1_pt3 = (TH1F*)in1->Get("hpt_reco_3");
    TH1F *hEff1_pt4 = (TH1F*)in1->Get("hpt_reco_4");
    TH1F *hEff1_pt5 = (TH1F*)in1->Get("hpt_reco_5");
    TH1F *hEff1_pt6 = (TH1F*)in1->Get("hpt_reco_6");


    TH1F *hGen1_pt1 = (TH1F*)in1->Get("hpt_gen_1");
    TH1F *hGen1_pt2 = (TH1F*)in1->Get("hpt_gen_2");
    TH1F *hGen1_pt3 = (TH1F*)in1->Get("hpt_gen_3");
    TH1F *hGen1_pt4 = (TH1F*)in1->Get("hpt_gen_4");
    TH1F *hGen1_pt5 = (TH1F*)in1->Get("hpt_gen_5");
    TH1F *hGen1_pt6 = (TH1F*)in1->Get("hpt_gen_6");
    */


    hEff1_pt1->SetMarkerStyle(20);
    hEff1_pt1->SetMarkerColor(kRed+1);
    hEff1_pt1->SetLineColor(kRed+1);
    hEff1_pt2->SetMarkerStyle(30);
    hEff1_pt2->SetMarkerColor(kOrange+1);
    hEff1_pt2->SetLineColor(kOrange+1);
    hEff1_pt3->SetMarkerStyle(26);
    hEff1_pt3->SetMarkerColor(kGreen+1);
    hEff1_pt3->SetLineColor(kGreen+1);
    hEff1_pt4->SetMarkerStyle(27);
    hEff1_pt4->SetMarkerColor(kBlue+1);
    hEff1_pt4->SetLineColor(kBlue+1);
    hEff1_pt5->SetMarkerStyle(28);
    hEff1_pt5->SetMarkerColor(kViolet+1);
    hEff1_pt5->SetLineColor(kViolet+1);
    hEff1_pt6->SetMarkerStyle(31);
    hEff1_pt6->SetMarkerColor(kAzure+1);
    hEff1_pt6->SetLineColor(kAzure+1);

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();

    TLegend *lg1 = new TLegend(0.55,0.20,0.73,0.45);
    lg1->SetFillStyle(0);
    lg1->SetFillColor(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.035);
    //lg1->AddEntry(hEff1_pt1,"|y|: 0-2.4, 10-60 %","lep");
    //lg1->AddEntry(hEff1_pt2,"|y|: 0-0.6, 10-60 %","lep");
    lg1->AddEntry(hEff1_pt3,"|y|: 0.0-1.2, 10-60 %","lep");
    lg1->AddEntry(hEff1_pt4,"|y|: 1.2-1.6, 10-60 %","lep");
    lg1->AddEntry(hEff1_pt5,"|y|: 1.6-1.8, 10-60 %","lep");
    lg1->AddEntry(hEff1_pt6,"|y|: 1.8-2.4, 10-60 %","lep");


    /*
    hEff1_pt1->Divide(hGen1_pt1);
    hEff1_pt2->Divide(hGen1_pt2);
    hEff1_pt3->Divide(hGen1_pt3);
    hEff1_pt4->Divide(hGen1_pt4);
    hEff1_pt5->Divide(hGen1_pt5);
    hEff1_pt6->Divide(hGen1_pt6);
    */


    TCanvas *c1 = new TCanvas("c1","",660,600);
    c1->cd();

    hy1_eff->GetYaxis()->SetRangeUser(0.0,1.2);
    hy1_eff->Draw("E");
    //lt1->DrawLatex(0.20,0.87,"Trigger Matching");
    lt1->DrawLatex(0.20,0.87,"CMS Simulation");
    lt1->DrawLatex(0.20,0.80,"Prompt #psi(2S)");
    lt1->DrawLatex(0.60,0.87,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
    c1->SaveAs("eff_total_y.png");



	c1->cd();
    hEff1_pt3->GetXaxis()->SetRangeUser(0.0,50.0);
    hEff1_pt3->GetYaxis()->SetRangeUser(0.0,1.2);
    hEff1_pt3->GetYaxis()->SetTitle("Reconstruction Efficiency");
    hEff1_pt3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEff1_pt3->GetYaxis()->CenterTitle();
    hEff1_pt3->GetXaxis()->CenterTitle();
    hEff1_pt3->Draw("E");
    //hEff1_pt2->Draw("E same");
    //hEff1_pt3->Draw("E same");
    hEff1_pt4->Draw("E same");
    hEff1_pt5->Draw("E same");
    hEff1_pt6->Draw("E same");

    lg1->Draw("same");
    //lt1->DrawLatex(0.20,0.87,"Trigger Matching");
    lt1->DrawLatex(0.20,0.87,"CMS Simulation");
    lt1->DrawLatex(0.20,0.80,"Prompt #psi(2S)");
    lt1->DrawLatex(0.60,0.87,"PbPb #sqrt{s_{NN}} = 5.02 TeV");

    c1->SaveAs("eff_psi2s_total_pt.png");
    c1->SaveAs("eff_psi2s_total_pt.pdf");

 
    TFile *out = new TFile("eff_psi2s_pbpb_ptwgt1.root","RECREATE");
    out->cd();
    
    hEff1_pt1->SetName("hEff_pt1");
    hEff1_pt2->SetName("hEff_pt2");
    hEff1_pt3->SetName("hEff_pt3");
    hEff1_pt4->SetName("hEff_pt4");
    hEff1_pt5->SetName("hEff_pt5");
    hEff1_pt6->SetName("hEff_pt6");
    hy1_eff->SetName("hy_eff");

    hEff1_pt1->Write();
    hEff1_pt2->Write();
    hEff1_pt3->Write();
    hEff1_pt4->Write();
    hEff1_pt5->Write();
    hEff1_pt6->Write();
    hy1_eff->Write();
    out->Close();


}
