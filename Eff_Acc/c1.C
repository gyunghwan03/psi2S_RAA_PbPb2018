void c1()
{
//=========Macro generated from canvas: c1/
//=========  (Tue Mar  1 14:41:10 2022) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "",164,187,660,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-9.146341,-0.2024096,51.82927,1.243373);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(0);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.03);
   c1->SetTopMargin(0.03);
   c1->SetBottomMargin(0.14);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   Double_t xAxis1[10] = {0, 3, 4.5, 6, 7, 8, 9, 10, 15, 50}; 
   
   TH1D *hEff_pt3__1 = new TH1D("hEff_pt3__1","Eff: Rapidity 1.2-1.6",9, xAxis1);
   hEff_pt3__1->SetBinContent(4,0.1874563);
   hEff_pt3__1->SetBinContent(5,0.2467119);
   hEff_pt3__1->SetBinContent(6,0.3050449);
   hEff_pt3__1->SetBinContent(7,0.3231029);
   hEff_pt3__1->SetBinContent(8,0.3832472);
   hEff_pt3__1->SetBinContent(9,0.4326394);
   hEff_pt3__1->SetBinError(4,0.02071991);
   hEff_pt3__1->SetBinError(5,0.01421113);
   hEff_pt3__1->SetBinError(6,0.01503026);
   hEff_pt3__1->SetBinError(7,0.01456593);
   hEff_pt3__1->SetBinError(8,0.008715345);
   hEff_pt3__1->SetBinError(9,0.008363287);
   hEff_pt3__1->SetMinimum(0);
   hEff_pt3__1->SetMaximum(1.2);
   hEff_pt3__1->SetEntries(28985);
   hEff_pt3__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#00cc00");
   hEff_pt3__1->SetLineColor(ci);

   ci = TColor::GetColor("#00cc00");
   hEff_pt3__1->SetMarkerColor(ci);
   hEff_pt3__1->SetMarkerStyle(26);
   hEff_pt3__1->SetMarkerSize(1.2);
   hEff_pt3__1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hEff_pt3__1->GetXaxis()->SetRange(1,9);
   hEff_pt3__1->GetXaxis()->CenterTitle(true);
   hEff_pt3__1->GetXaxis()->SetLabelFont(42);
   hEff_pt3__1->GetXaxis()->SetTitleSize(0.048);
   hEff_pt3__1->GetXaxis()->SetTitleOffset(1.2);
   hEff_pt3__1->GetXaxis()->SetTitleFont(42);
   hEff_pt3__1->GetYaxis()->SetTitle("Reconstruction Efficiency");
   hEff_pt3__1->GetYaxis()->CenterTitle(true);
   hEff_pt3__1->GetYaxis()->SetLabelFont(42);
   hEff_pt3__1->GetYaxis()->SetTitleSize(0.048);
   hEff_pt3__1->GetYaxis()->SetTitleOffset(1.4);
   hEff_pt3__1->GetYaxis()->SetTitleFont(42);
   hEff_pt3__1->GetZaxis()->SetLabelFont(42);
   hEff_pt3__1->GetZaxis()->SetTitleSize(0.048);
   hEff_pt3__1->GetZaxis()->SetTitleOffset(1);
   hEff_pt3__1->GetZaxis()->SetTitleFont(42);
   hEff_pt3__1->Draw("E");
   Double_t xAxis2[10] = {0, 3, 4.5, 6, 7, 8, 9, 10, 15, 50}; 
   
   TH1D *hEff_pt4__2 = new TH1D("hEff_pt4__2","Eff: Rapidity 1.6-1.8",9, xAxis2);
   hEff_pt4__2->SetBinContent(2,0.1232832);
   hEff_pt4__2->SetBinContent(3,0.2595411);
   hEff_pt4__2->SetBinContent(4,0.280012);
   hEff_pt4__2->SetBinContent(5,0.3862885);
   hEff_pt4__2->SetBinContent(6,0.3847538);
   hEff_pt4__2->SetBinContent(7,0.4106959);
   hEff_pt4__2->SetBinContent(8,0.4770656);
   hEff_pt4__2->SetBinContent(9,0.5565753);
   hEff_pt4__2->SetBinError(2,0.03202246);
   hEff_pt4__2->SetBinError(3,0.02213403);
   hEff_pt4__2->SetBinError(4,0.02079614);
   hEff_pt4__2->SetBinError(5,0.02048252);
   hEff_pt4__2->SetBinError(6,0.02062645);
   hEff_pt4__2->SetBinError(7,0.02175986);
   hEff_pt4__2->SetBinError(8,0.01262612);
   hEff_pt4__2->SetBinError(9,0.01165);
   hEff_pt4__2->SetEntries(16815);
   hEff_pt4__2->SetStats(0);

   ci = TColor::GetColor("#0000cc");
   hEff_pt4__2->SetLineColor(ci);

   ci = TColor::GetColor("#0000cc");
   hEff_pt4__2->SetMarkerColor(ci);
   hEff_pt4__2->SetMarkerStyle(27);
   hEff_pt4__2->SetMarkerSize(1.2);
   hEff_pt4__2->GetXaxis()->SetLabelFont(42);
   hEff_pt4__2->GetXaxis()->SetTitleSize(0.048);
   hEff_pt4__2->GetXaxis()->SetTitleOffset(1.2);
   hEff_pt4__2->GetXaxis()->SetTitleFont(42);
   hEff_pt4__2->GetYaxis()->SetLabelFont(42);
   hEff_pt4__2->GetYaxis()->SetTitleSize(0.048);
   hEff_pt4__2->GetYaxis()->SetTitleOffset(1.4);
   hEff_pt4__2->GetYaxis()->SetTitleFont(42);
   hEff_pt4__2->GetZaxis()->SetLabelFont(42);
   hEff_pt4__2->GetZaxis()->SetTitleSize(0.048);
   hEff_pt4__2->GetZaxis()->SetTitleOffset(1);
   hEff_pt4__2->GetZaxis()->SetTitleFont(42);
   hEff_pt4__2->Draw("E same");
   Double_t xAxis3[10] = {0, 3, 4.5, 6, 7, 8, 9, 10, 15, 50}; 
   
   TH1D *hEff_pt5__3 = new TH1D("hEff_pt5__3","Eff: Rapidity 1.8-2.4",9, xAxis3);
   hEff_pt5__3->SetBinContent(2,0.1105638);
   hEff_pt5__3->SetBinContent(3,0.2227046);
   hEff_pt5__3->SetBinContent(4,0.2801846);
   hEff_pt5__3->SetBinContent(5,0.3668636);
   hEff_pt5__3->SetBinContent(6,0.4032296);
   hEff_pt5__3->SetBinContent(7,0.4192827);
   hEff_pt5__3->SetBinContent(8,0.4323871);
   hEff_pt5__3->SetBinContent(9,0.4786017);
   hEff_pt5__3->SetBinError(2,0.01226563);
   hEff_pt5__3->SetBinError(3,0.01348316);
   hEff_pt5__3->SetBinError(4,0.01643);
   hEff_pt5__3->SetBinError(5,0.01615713);
   hEff_pt5__3->SetBinError(6,0.01589473);
   hEff_pt5__3->SetBinError(7,0.01571094);
   hEff_pt5__3->SetBinError(8,0.008822107);
   hEff_pt5__3->SetBinError(9,0.008084912);
   hEff_pt5__3->SetEntries(33553);
   hEff_pt5__3->SetStats(0);

   ci = TColor::GetColor("#9933ff");
   hEff_pt5__3->SetLineColor(ci);

   ci = TColor::GetColor("#9933ff");
   hEff_pt5__3->SetMarkerColor(ci);
   hEff_pt5__3->SetMarkerStyle(28);
   hEff_pt5__3->SetMarkerSize(1.2);
   hEff_pt5__3->GetXaxis()->SetLabelFont(42);
   hEff_pt5__3->GetXaxis()->SetTitleSize(0.048);
   hEff_pt5__3->GetXaxis()->SetTitleOffset(1.2);
   hEff_pt5__3->GetXaxis()->SetTitleFont(42);
   hEff_pt5__3->GetYaxis()->SetLabelFont(42);
   hEff_pt5__3->GetYaxis()->SetTitleSize(0.048);
   hEff_pt5__3->GetYaxis()->SetTitleOffset(1.4);
   hEff_pt5__3->GetYaxis()->SetTitleFont(42);
   hEff_pt5__3->GetZaxis()->SetLabelFont(42);
   hEff_pt5__3->GetZaxis()->SetTitleSize(0.048);
   hEff_pt5__3->GetZaxis()->SetTitleOffset(1);
   hEff_pt5__3->GetZaxis()->SetTitleFont(42);
   hEff_pt5__3->Draw("E same");
   Double_t xAxis4[10] = {0, 3, 4.5, 6, 7, 8, 9, 10, 15, 50}; 
   
   TH1D *hEff_pt6__4 = new TH1D("hEff_pt6__4","hpt_reco_5",9, xAxis4);
   hEff_pt6__4->SetBinContent(2,9977.409);
   hEff_pt6__4->SetBinContent(3,10945.08);
   hEff_pt6__4->SetBinContent(4,9423.335);
   hEff_pt6__4->SetBinContent(5,11161.67);
   hEff_pt6__4->SetBinContent(6,7486.514);
   hEff_pt6__4->SetBinContent(7,4513.3);
   hEff_pt6__4->SetBinContent(8,5914.456);
   hEff_pt6__4->SetBinContent(9,954.3452);
   hEff_pt6__4->SetBinError(2,1186.738);
   hEff_pt6__4->SetBinError(3,750.044);
   hEff_pt6__4->SetBinError(4,651.4982);
   hEff_pt6__4->SetBinError(5,633.0263);
   hEff_pt6__4->SetBinError(6,387.9437);
   hEff_pt6__4->SetBinError(7,219.6526);
   hEff_pt6__4->SetBinError(8,157.7048);
   hEff_pt6__4->SetBinError(9,21.57804);
   hEff_pt6__4->SetEntries(14204);
   hEff_pt6__4->SetStats(0);

   ci = TColor::GetColor("#3399ff");
   hEff_pt6__4->SetLineColor(ci);

   ci = TColor::GetColor("#3399ff");
   hEff_pt6__4->SetMarkerColor(ci);
   hEff_pt6__4->SetMarkerStyle(31);
   hEff_pt6__4->SetMarkerSize(1.2);
   hEff_pt6__4->GetXaxis()->SetLabelFont(42);
   hEff_pt6__4->GetXaxis()->SetTitleSize(0.048);
   hEff_pt6__4->GetXaxis()->SetTitleOffset(1.2);
   hEff_pt6__4->GetXaxis()->SetTitleFont(42);
   hEff_pt6__4->GetYaxis()->SetLabelFont(42);
   hEff_pt6__4->GetYaxis()->SetTitleSize(0.048);
   hEff_pt6__4->GetYaxis()->SetTitleOffset(1.4);
   hEff_pt6__4->GetYaxis()->SetTitleFont(42);
   hEff_pt6__4->GetZaxis()->SetLabelFont(42);
   hEff_pt6__4->GetZaxis()->SetTitleSize(0.048);
   hEff_pt6__4->GetZaxis()->SetTitleOffset(1);
   hEff_pt6__4->GetZaxis()->SetTitleFont(42);
   hEff_pt6__4->Draw("E same");
   
   TLegend *leg = new TLegend(0.550152,0.1751313,0.731003,0.4238179,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.035);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hEff_pt3","|y|: 0.0-1.2, 0-10 %","lep");

   ci = TColor::GetColor("#00cc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00cc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hEff_pt4","|y|: 1.2-1.6, 0-10 %","lep");

   ci = TColor::GetColor("#0000cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(27);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hEff_pt5","|y|: 1.6-1.8, 0-10 %","lep");

   ci = TColor::GetColor("#9933ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#9933ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(28);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hEff_pt6","|y|: 1.8-2.4, 0-10 %","lep");

   ci = TColor::GetColor("#3399ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3399ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(31);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   leg->Draw();
   TLatex *   tex = new TLatex(0.2,0.87,"CMS Simulation");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8,"Prompt #psi(2S)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6,0.87,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
