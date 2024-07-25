void cRAA()
{
//=========Macro generated from canvas: cRAA/
//=========  (Tue Jul 16 20:43:30 2024) by ROOT version 6.18/04
   TCanvas *cRAA = new TCanvas("cRAA", "",2883,362,700,700);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cRAA->Range(-7.920792,-0.2282927,41.58416,1.527805);
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
   7.75,
   10.5,
   13.5,
   17.5,
   22.5,
   32.5};
   Double_t Graph0_fy1001[6] = {
   0.2434115,
   0.1785764,
   0.1934771,
   0.2114064,
   0.2413684,
   0.2736468};
   Double_t Graph0_fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[6] = {
   0.03577576,
   0.01294611,
   0.01636791,
   0.02284889,
   0.04192479,
   0.05337822};
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
   gre->SetMarkerSize(1.4);
   
   TH1F *Graph_Graph01001 = new TH1F("Graph_Graph01001","",100,0,40);
   Graph_Graph01001->SetMinimum(0);
   Graph_Graph01001->SetMaximum(1.44);
   Graph_Graph01001->SetDirectory(0);
   Graph_Graph01001->SetStats(0);
   Graph_Graph01001->SetLineStyle(0);
   Graph_Graph01001->SetMarkerStyle(20);
   Graph_Graph01001->GetXaxis()->SetTitle("p_{T} (GeV/c)");
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
   7.75,
   10.5,
   13.5,
   17.5,
   22.5,
   32.5};
   Double_t Graph1_fy1002[6] = {
   0.2754807,
   0.217654,
   0.2577349,
   0.2550198,
   0.3157402,
   0.3558492};
   Double_t Graph1_fex1002[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1002[6] = {
   0.04177038,
   0.01639407,
   0.02114623,
   0.02562665,
   0.04474357,
   0.05656255};
   gre = new TGraphErrors(6,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#660000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#660000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(1.4);
   
   TH1F *Graph_Graph11002 = new TH1F("Graph_Graph11002","Graph",100,5.275,34.975);
   Graph_Graph11002->SetMinimum(0.1801448);
   Graph_Graph11002->SetMaximum(0.4335269);
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
   7.75,
   10.5,
   13.5,
   17.5,
   22.5,
   32.5};
   Double_t Graph2_fy1003[6] = {
   0.2434115,
   0.1785764,
   0.1934771,
   0.2114064,
   0.2413684,
   0.2736468};
   Double_t Graph2_fex1003[6] = {
   1.25,
   1.5,
   1.5,
   2.5,
   2.5,
   7.5};
   Double_t Graph2_fey1003[6] = {
   0.02599687,
   0.01977726,
   0.03341218,
   0.02190423,
   0.02362113,
   0.02887368};
   gre = new TGraphErrors(6,Graph2_fx1003,Graph2_fy1003,Graph2_fex1003,Graph2_fey1003);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");

   ci = 1183;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.4);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#3333ff");
   gre->SetLineColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph21003 = new TH1F("Graph_Graph21003","Graph",100,0,40);
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
   7.75,
   10.5,
   13.5,
   17.5,
   22.5,
   32.5};
   Double_t Graph3_fy1004[6] = {
   0.2754807,
   0.217654,
   0.2577349,
   0.2550198,
   0.3157402,
   0.3558492};
   Double_t Graph3_fex1004[6] = {
   1.25,
   1.5,
   1.5,
   2.5,
   2.5,
   7.5};
   Double_t Graph3_fey1004[6] = {
   0.04852139,
   0.03220003,
   0.04826296,
   0.02937523,
   0.0272676,
   0.03492241};
   gre = new TGraphErrors(6,Graph3_fx1004,Graph3_fy1004,Graph3_fex1004,Graph3_fey1004);
   gre->SetName("Graph3");
   gre->SetTitle("Graph");

   ci = 1184;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.4);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#ff3333");
   gre->SetLineColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph31004 = new TH1F("Graph_Graph31004","Graph",100,0,40);
   Graph_Graph31004->SetMinimum(0);
   Graph_Graph31004->SetMaximum(1.44);
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
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Non Prompt #psi(2S)","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#660000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,40,1);
   line->SetLineStyle(7);
   line->Draw();
   TBox *box = new TBox(0,0.9683139,2,1.031686);

   ci = 1185;
   color = new TColor(ci, 0.3, 0.3, 0.3, " ", 0.8);
   box->SetFillColor(ci);
   box->Draw();
   TLatex *   tex = new TLatex(0.25,0.85,"Cent. 0-90%");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(0.25,0.779,"|y| < 1.6");
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
