void cRAA_Npart()
{
//=========Macro generated from canvas: cRAA/
//=========  (Tue Jul 16 20:49:19 2024) by ROOT version 6.18/04
   TCanvas *cRAA = new TCanvas("cRAA", "",0,53,700,700);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cRAA->Range(-79.20792,-0.2282927,415.8416,1.527805);
   cRAA->SetFillColor(0);
   cRAA->SetBorderMode(0);
   cRAA->SetBorderSize(2);
   cRAA->SetTickx(1);
   cRAA->SetTicky(1);
   cRAA->SetLeftMargin(0.16);
   cRAA->SetRightMargin(0.032);
   cRAA->SetTopMargin(0.05);
   cRAA->SetBottomMargin(0.13);
   cRAA->SetFrameFillStyle(0);
   cRAA->SetFrameBorderMode(0);
   cRAA->SetFrameFillStyle(0);
   cRAA->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[6] = {
   27.12,
   87.19,
   131,
   188.2,
   262.3,
   356.9};
   Double_t Graph0_fy1001[6] = {
   0.4293478,
   0.4964902,
   0.3728716,
   0.4751304,
   0.6624808,
   0.4181731};
   Double_t Graph0_fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[6] = {
   0.1164218,
   0.08055286,
   0.08985873,
   0.1165608,
   0.2366457,
   0.1108037};
   TGraphErrors *gre = new TGraphErrors(6,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("");
   gre->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#000099");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph01001 = new TH1F("Graph_Graph01001","",100,0,400);
   Graph_Graph01001->SetMinimum(0);
   Graph_Graph01001->SetMaximum(1.44);
   Graph_Graph01001->SetDirectory(0);
   Graph_Graph01001->SetStats(0);
   Graph_Graph01001->SetLineStyle(0);
   Graph_Graph01001->SetMarkerStyle(20);
   Graph_Graph01001->GetXaxis()->SetTitle("<N_{Part}>");
   Graph_Graph01001->GetXaxis()->CenterTitle(true);
   Graph_Graph01001->GetXaxis()->SetLabelFont(42);
   Graph_Graph01001->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph01001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph01001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph01001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph01001->GetXaxis()->SetTitleFont(42);
   Graph_Graph01001->GetYaxis()->SetTitle("R_{AA}");
   Graph_Graph01001->GetYaxis()->CenterTitle(true);
   Graph_Graph01001->GetYaxis()->SetLabelFont(42);
   Graph_Graph01001->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph01001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph01001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph01001->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph01001->GetYaxis()->SetTitleFont(42);
   Graph_Graph01001->GetZaxis()->SetLabelFont(42);
   Graph_Graph01001->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph01001->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph01001->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph01001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph01001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph01001);
   
   gre->Draw("ap");
   
   Double_t Graph1_fx1002[6] = {
   27.12,
   87.19,
   131,
   188.2,
   262.3,
   356.9};
   Double_t Graph1_fy1002[6] = {
   0.5175139,
   0.6739952,
   0.3850315,
   0.5122114,
   0.6024247,
   0.2802176};
   Double_t Graph1_fex1002[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1002[6] = {
   0.1455266,
   0.1194875,
   0.09696515,
   0.1311,
   0.2173889,
   0.07986358};
   gre = new TGraphErrors(6,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#660000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#660000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph11002 = new TH1F("Graph_Graph11002","Graph",100,0,389.878);
   Graph_Graph11002->SetMinimum(0.138408);
   Graph_Graph11002->SetMaximum(0.8817596);
   Graph_Graph11002->SetDirectory(0);
   Graph_Graph11002->SetStats(0);
   Graph_Graph11002->SetLineStyle(0);
   Graph_Graph11002->SetMarkerStyle(20);
   Graph_Graph11002->GetXaxis()->SetLabelFont(42);
   Graph_Graph11002->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph11002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph11002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph11002->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph11002->GetXaxis()->SetTitleFont(42);
   Graph_Graph11002->GetYaxis()->SetLabelFont(42);
   Graph_Graph11002->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph11002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph11002->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph11002->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph11002->GetYaxis()->SetTitleFont(42);
   Graph_Graph11002->GetZaxis()->SetLabelFont(42);
   Graph_Graph11002->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph11002->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph11002->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph11002->GetZaxis()->SetTitleOffset(1);
   Graph_Graph11002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph11002);
   
   gre->Draw("p");
   
   Double_t Graph2_fx1003[6] = {
   27.12,
   87.19,
   131,
   188.2,
   262.3,
   356.9};
   Double_t Graph2_fy1003[6] = {
   0.4293478,
   0.4964902,
   0.3728716,
   0.4751304,
   0.6624808,
   0.4181731};
   Double_t Graph2_fex1003[6] = {
   4.3,
   4.3,
   4.3,
   4.3,
   4.3,
   4.3};
   Double_t Graph2_fey1003[6] = {
   0.06700721,
   0.106655,
   0.06215021,
   0.08402691,
   0.1229465,
   0.08527856};
   gre = new TGraphErrors(6,Graph2_fx1003,Graph2_fy1003,Graph2_fex1003,Graph2_fey1003);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");

   ci = 1179;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.4);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#3333ff");
   gre->SetLineColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph21003 = new TH1F("Graph_Graph21003","Graph",100,0,395.038);
   Graph_Graph21003->SetMinimum(0);
   Graph_Graph21003->SetMaximum(1.44);
   Graph_Graph21003->SetDirectory(0);
   Graph_Graph21003->SetStats(0);
   Graph_Graph21003->SetLineStyle(0);
   Graph_Graph21003->SetMarkerStyle(20);
   Graph_Graph21003->GetXaxis()->SetLabelFont(42);
   Graph_Graph21003->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph21003->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph21003->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph21003->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph21003->GetXaxis()->SetTitleFont(42);
   Graph_Graph21003->GetYaxis()->SetLabelFont(42);
   Graph_Graph21003->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph21003->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph21003->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph21003->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph21003->GetYaxis()->SetTitleFont(42);
   Graph_Graph21003->GetZaxis()->SetLabelFont(42);
   Graph_Graph21003->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph21003->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph21003->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph21003->GetZaxis()->SetTitleOffset(1);
   Graph_Graph21003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph21003);
   
   gre->Draw("5");
   
   Double_t Graph3_fx1004[6] = {
   27.12,
   87.19,
   131,
   188.2,
   262.3,
   356.9};
   Double_t Graph3_fy1004[6] = {
   0.5175139,
   0.6739952,
   0.3850315,
   0.5122114,
   0.6024247,
   0.2802176};
   Double_t Graph3_fex1004[6] = {
   4.3,
   4.3,
   4.3,
   4.3,
   4.3,
   4.3};
   Double_t Graph3_fey1004[6] = {
   0.1546556,
   0.2036263,
   0.1172443,
   0.1614459,
   0.1795797,
   0.09753354};
   gre = new TGraphErrors(6,Graph3_fx1004,Graph3_fy1004,Graph3_fex1004,Graph3_fey1004);
   gre->SetName("Graph3");
   gre->SetTitle("Graph");

   ci = 1180;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.4);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#ff3333");
   gre->SetLineColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph31004 = new TH1F("Graph_Graph31004","Graph",100,0,395.038);
   Graph_Graph31004->SetMinimum(0.1131903);
   Graph_Graph31004->SetMaximum(0.9471153);
   Graph_Graph31004->SetDirectory(0);
   Graph_Graph31004->SetStats(0);
   Graph_Graph31004->SetLineStyle(0);
   Graph_Graph31004->SetMarkerStyle(20);
   Graph_Graph31004->GetXaxis()->SetLabelFont(42);
   Graph_Graph31004->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph31004->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph31004->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph31004->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph31004->GetXaxis()->SetTitleFont(42);
   Graph_Graph31004->GetYaxis()->SetLabelFont(42);
   Graph_Graph31004->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph31004->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph31004->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph31004->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph31004->GetYaxis()->SetTitleFont(42);
   Graph_Graph31004->GetZaxis()->SetLabelFont(42);
   Graph_Graph31004->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph31004->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph31004->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph31004->GetZaxis()->SetTitleOffset(1);
   Graph_Graph31004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph31004);
   
   gre->Draw("5");
   
   TLegend *leg = new TLegend(0.68,0.72,0.8,0.82,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.029);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(4000);
   TLegendEntry *entry=leg->AddEntry("Graph0","Prompt #psi(2S)","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Non Prompt #psi(2S)","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#660000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,400,1);
   line->SetLineStyle(7);
   line->Draw();
   TBox *box = new TBox(0,0.9771962,15,1.022804);

   ci = 1181;
   color = new TColor(ci, 0.3, 0.3, 0.3, " ", 0.8);
   box->SetFillColor(ci);
   box->Draw();
   TLatex *   tex = new TLatex(0.2,0.85,"3.5 < p_{T} < 40 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(0.2,0.799,"1.6 < |y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(0.968,0.9528,"PbPb 1.61 nb^{-1}, pp 300.0 pb^{-1} (5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03735);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.92164,0.9095026,"CMS");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.0466875);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.92164,0.8534775,"Preliminary");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0354825);
   tex->SetLineWidth(2);
   tex->Draw();
   cRAA->Modified();
   cRAA->cd();
   cRAA->SetSelected(cRAA);
}
