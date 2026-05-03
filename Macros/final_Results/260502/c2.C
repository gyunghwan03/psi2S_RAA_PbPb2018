void c2()
{
//=========Macro generated from canvas: c2/
//=========  (Wed Aug 21 15:05:30 2024) by ROOT version 6.18/04
   TCanvas *c2 = new TCanvas("c2", "",0,72,900,800);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c2->Range(-79.20792,-0.2282927,415.8416,1.527805);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetTickx(1);
   c2->SetTicky(1);
   c2->SetLeftMargin(0.16);
   c2->SetRightMargin(0.032);
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.13);
   c2->SetFrameFillStyle(0);
   c2->SetFrameBorderMode(0);
   c2->SetFrameFillStyle(0);
   c2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[4] = {
   27.12,
   109.1,
   225.2,
   356.9};
   Double_t Graph0_fy1001[4] = {
   0.5620092,
   0.4795678,
   0.2714902,
   0.2366446};
   Double_t Graph0_fex1001[4] = {
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[4] = {
   0.02864934,
   0.01648645,
   0.005321751,
   0.005154834};
   TGraphErrors *gre = new TGraphErrors(4,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
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
   
   Double_t Graph1_fx1002[3] = {
   32.7,
   160.3,
   311.5};
   Double_t Graph1_fy1002[3] = {
   0.719,
   0.558,
   0.372};
   Double_t Graph1_fex1002[3] = {
   0,
   0,
   0};
   Double_t Graph1_fey1002[3] = {
   0.024,
   0.016,
   0.01};
   gre = new TGraphErrors(3,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(15);
   gre->SetMarkerColor(15);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.4);
   
   TH1F *Graph_Graph11002 = new TH1F("Graph_Graph11002","",100,0,400);
   Graph_Graph11002->SetMinimum(0);
   Graph_Graph11002->SetMaximum(1.44);
   Graph_Graph11002->SetDirectory(0);
   Graph_Graph11002->SetStats(0);
   Graph_Graph11002->SetLineStyle(0);
   Graph_Graph11002->SetMarkerStyle(20);
   Graph_Graph11002->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   Graph_Graph11002->GetXaxis()->CenterTitle(true);
   Graph_Graph11002->GetXaxis()->SetLabelFont(42);
   Graph_Graph11002->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph11002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph11002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph11002->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph11002->GetXaxis()->SetTitleFont(42);
   Graph_Graph11002->GetYaxis()->SetTitle("R_{AA}");
   Graph_Graph11002->GetYaxis()->CenterTitle(true);
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
   
   TLegend *leg = new TLegend(0.3351893,0.6025806,0.7048998,0.6825806,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(20);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","#bf{Prompt J/#psi}, Cent. 0-90%","lpf");
   entry->SetFillStyle(1000);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(43);
   entry=leg->AddEntry("Graph1","#bf{Prompt J/#psi}, HIN-16-025, 3 < p_{T} < 30 GeV/c, Cent. 0-100%","lpf");
   entry->SetFillStyle(1000);
   entry->SetLineColor(15);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(15);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(43);
   leg->Draw();
   TLine *line = new TLine(0,1,400,1);
   line->SetLineStyle(7);
   line->Draw();
   TLatex *   tex = new TLatex(0.21,0.87,"3.5 < p_{T} < 40 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->Draw();
      tex = new TLatex(0.21,0.819,"1.6 < |y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->Draw();
      tex = new TLatex(0.968,0.9528,"PbPb 1.7 nb^{-1}, pp 300.2 pb^{-1} (5.02 TeV)");
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
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);
}
