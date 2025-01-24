void c1()
{
//=========Macro generated from canvas: c1/
//=========  (Wed Aug 21 15:21:20 2024) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "",0,72,900,802);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-7.920792,-0.2282927,41.58416,1.527805);
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
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[6] = {
   7.75,
   10.5,
   13.5,
   17.5,
   22.5,
   32.5};
   Double_t Graph0_fy1001[6] = {
   0.3085219,
   0.297992,
   0.3146375,
   0.3308058,
   0.3886052,
   0.429906};
   Double_t Graph0_fex1001[6] = {
   1.25,
   1.5,
   1.5,
   2.5,
   2.5,
   7.5};
   Double_t Graph0_fey1001[6] = {
   0.006563112,
   0.002053768,
   0.003353051,
   0.004840387,
   0.01034174,
   0.0150712};
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
   
   Double_t Graph1D_y1_fx3001[5] = {
   7.75,
   10.5,
   13.5,
   17.5,
   25};
   Double_t Graph1D_y1_fy3001[5] = {
   0.353,
   0.337,
   0.366,
   0.356,
   0.473};
   Double_t Graph1D_y1_felx3001[5] = {
   1.25,
   1.5,
   1.5,
   2.5,
   5};
   Double_t Graph1D_y1_fely3001[5] = {
   0.03590265,
   0.02213594,
   0.02325941,
   0.0254951,
   0.03911521};
   Double_t Graph1D_y1_fehx3001[5] = {
   1.25,
   1.5,
   1.5,
   2.5,
   5};
   Double_t Graph1D_y1_fehy3001[5] = {
   0.03590265,
   0.02213594,
   0.02325941,
   0.0254951,
   0.03911521};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5,Graph1D_y1_fx3001,Graph1D_y1_fy3001,Graph1D_y1_felx3001,Graph1D_y1_fehx3001,Graph1D_y1_fely3001,Graph1D_y1_fehy3001);
   grae->SetName("Graph1D_y1");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->SetLineColor(15);
   grae->SetMarkerColor(15);
   grae->SetMarkerStyle(25);
   grae->SetMarkerSize(1.4);
   
   TH1F *Graph_Graph1D_y13001 = new TH1F("Graph_Graph1D_y13001","",100,0,40);
   Graph_Graph1D_y13001->SetMinimum(0);
   Graph_Graph1D_y13001->SetMaximum(1.44);
   Graph_Graph1D_y13001->SetDirectory(0);
   Graph_Graph1D_y13001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1D_y13001->SetLineColor(ci);
   Graph_Graph1D_y13001->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   Graph_Graph1D_y13001->GetXaxis()->CenterTitle(true);
   Graph_Graph1D_y13001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1D_y13001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1D_y13001->GetXaxis()->SetTitleOffset(1);
   Graph_Graph1D_y13001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1D_y13001->GetYaxis()->SetTitle("R_{AA}");
   Graph_Graph1D_y13001->GetYaxis()->CenterTitle(true);
   Graph_Graph1D_y13001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1D_y13001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1D_y13001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1D_y13001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1D_y13001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1D_y13001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1D_y13001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1D_y13001);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.7026726,0.6103226,0.7628062,0.6903226,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(20);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","#bf{Prompt J/#psi}, New Data","lpf");
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
   entry=leg->AddEntry("Graph1D_y1","#bf{Prompt J/#psi}, HIN-16-025","lpf");
   entry->SetFillStyle(1000);
   entry->SetLineColor(15);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(15);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(43);
   leg->Draw();
   TLine *line = new TLine(0,1,50,1);
   line->SetLineStyle(7);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.87,"6.5 < p_{T} < 40 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
   tex->Draw();
      tex = new TLatex(0.2,0.819,"|y| < 1.6");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(25);
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
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
