void c1()
{
//=========Macro generated from canvas: c1/
//=========  (Tue Jul 16 16:55:27 2024) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "",1120,1664,1000,702);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.032);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.13);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "",0,0,0.83,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(-76.19047,-0.2282927,400,1.527805);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetLeftMargin(0.16);
   pad1->SetRightMargin(0);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.13);
   pad1->SetFrameFillStyle(0);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameFillStyle(0);
   pad1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[8] = {
   27.12,
   87.19,
   131,
   188.2,
   241,
   283.6,
   331.5,
   382.3};
   Double_t Graph0_fy1001[8] = {
   0.4698968,
   0.3463779,
   0.3309435,
   0.2393069,
   0.1832275,
   0.2244001,
   0.2190444,
   0.1535174};
   Double_t Graph0_fex1001[8] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[8] = {
   0.05198155,
   0.04865679,
   0.03928405,
   0.02804241,
   0.03936525,
   0.03626783,
   0.03261526,
   0.02738867};
   TGraphErrors *gre = new TGraphErrors(8,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
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
   
   Double_t Graph1_fx1002[8] = {
   27.12,
   87.19,
   131,
   188.2,
   241,
   283.6,
   331.5,
   382.3};
   Double_t Graph1_fy1002[8] = {
   0.5974657,
   0.4271403,
   0.3742914,
   0.2825457,
   0.1840907,
   0.1878134,
   0.2761907,
   0.124554};
   Double_t Graph1_fex1002[8] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1002[8] = {
   0.0694333,
   0.06270638,
   0.04664967,
   0.03425445,
   0.040881,
   0.03255164,
   0.04213848,
   0.02389815};
   gre = new TGraphErrors(8,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#660000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#660000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph11002 = new TH1F("Graph_Graph11002","Graph",100,0,417.818);
   Graph_Graph11002->SetMinimum(0);
   Graph_Graph11002->SetMaximum(1.44);
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
   
   TLegend *leg = new TLegend(0.68,0.72,0.8,0.82,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.029);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
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
   TLatex *   tex = new TLatex(0.21,0.87,"6.5 < p_{T} < 40 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(0.21,0.819,"|y| < 1.6");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(1,0.9528,"PbPb 1.61 nb^{-1}, pp 300.0 pb^{-1} (5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03735);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9522,0.9095026,"CMS");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.0466875);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9522,0.8534775,"Preliminary");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0354825);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "",0.83,0,1,1);
   pad2->Draw();
   pad2->cd();
   pad2->Range(0,-0.2282927,1.033058,1.527805);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetTickx(1);
   pad2->SetTicky(1);
   pad2->SetLeftMargin(0);
   pad2->SetRightMargin(0.032);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.13);
   pad2->SetFrameFillStyle(0);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameFillStyle(0);
   pad2->SetFrameBorderMode(0);
   
   TH1D *hmidPR_int__1 = new TH1D("hmidPR_int__1","",3,0,1);
   hmidPR_int__1->SetBinContent(2,0.086);
   hmidPR_int__1->SetBinError(2,0.0182);
   hmidPR_int__1->SetMinimum(0);
   hmidPR_int__1->SetMaximum(1.44);
   hmidPR_int__1->SetEntries(1);

   ci = TColor::GetColor("#000099");
   hmidPR_int__1->SetLineColor(ci);
   hmidPR_int__1->SetLineStyle(0);

   ci = TColor::GetColor("#000099");
   hmidPR_int__1->SetMarkerColor(ci);
   hmidPR_int__1->SetMarkerStyle(20);
   hmidPR_int__1->SetMarkerSize(1.5);
   hmidPR_int__1->GetXaxis()->SetLabelFont(42);
   hmidPR_int__1->GetXaxis()->SetLabelOffset(0);
   hmidPR_int__1->GetXaxis()->SetLabelSize(0);
   hmidPR_int__1->GetXaxis()->SetTitleSize(0.06);
   hmidPR_int__1->GetXaxis()->SetTickLength(0);
   hmidPR_int__1->GetXaxis()->SetTitleOffset(0.9);
   hmidPR_int__1->GetXaxis()->SetTitleFont(42);
   hmidPR_int__1->GetYaxis()->SetLabelFont(42);
   hmidPR_int__1->GetYaxis()->SetLabelOffset(0);
   hmidPR_int__1->GetYaxis()->SetLabelSize(0.12);
   hmidPR_int__1->GetYaxis()->SetTitleSize(0.06);
   hmidPR_int__1->GetYaxis()->SetTitleOffset(1.25);
   hmidPR_int__1->GetYaxis()->SetTitleFont(42);
   hmidPR_int__1->GetZaxis()->SetLabelFont(42);
   hmidPR_int__1->GetZaxis()->SetLabelOffset(0.007);
   hmidPR_int__1->GetZaxis()->SetLabelSize(0.05);
   hmidPR_int__1->GetZaxis()->SetTitleSize(0.06);
   hmidPR_int__1->GetZaxis()->SetTitleOffset(1);
   hmidPR_int__1->GetZaxis()->SetTitleFont(42);
   hmidPR_int__1->Draw("PE");
   
   TH1D *hmidNP_int__2 = new TH1D("hmidNP_int__2","",3,0,1);
   hmidNP_int__2->SetBinContent(2,0.083);
   hmidNP_int__2->SetBinError(2,0.0175);
   hmidNP_int__2->SetMinimum(0);
   hmidNP_int__2->SetMaximum(1.44);
   hmidNP_int__2->SetEntries(1);

   ci = TColor::GetColor("#660000");
   hmidNP_int__2->SetLineColor(ci);
   hmidNP_int__2->SetLineStyle(0);

   ci = TColor::GetColor("#660000");
   hmidNP_int__2->SetMarkerColor(ci);
   hmidNP_int__2->SetMarkerStyle(21);
   hmidNP_int__2->SetMarkerSize(1.5);
   hmidNP_int__2->GetXaxis()->SetLabelFont(42);
   hmidNP_int__2->GetXaxis()->SetLabelOffset(0);
   hmidNP_int__2->GetXaxis()->SetLabelSize(0);
   hmidNP_int__2->GetXaxis()->SetTitleSize(0.06);
   hmidNP_int__2->GetXaxis()->SetTitleOffset(0.9);
   hmidNP_int__2->GetXaxis()->SetTitleFont(42);
   hmidNP_int__2->GetYaxis()->SetLabelFont(42);
   hmidNP_int__2->GetYaxis()->SetLabelOffset(0);
   hmidNP_int__2->GetYaxis()->SetLabelSize(0.12);
   hmidNP_int__2->GetYaxis()->SetTitleSize(0.06);
   hmidNP_int__2->GetYaxis()->SetTitleOffset(1.25);
   hmidNP_int__2->GetYaxis()->SetTitleFont(42);
   hmidNP_int__2->GetZaxis()->SetLabelFont(42);
   hmidNP_int__2->GetZaxis()->SetLabelOffset(0.007);
   hmidNP_int__2->GetZaxis()->SetLabelSize(0.05);
   hmidNP_int__2->GetZaxis()->SetTitleSize(0.06);
   hmidNP_int__2->GetZaxis()->SetTitleOffset(1);
   hmidNP_int__2->GetZaxis()->SetTitleFont(42);
   hmidNP_int__2->Draw("PE SAME");
   line = new TLine(0,1,1,1);
   line->SetLineStyle(7);
   line->Draw();
   pad2->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
