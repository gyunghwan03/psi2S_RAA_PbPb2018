#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../Gwak_raa_comparison/jpsi_raa_values.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace
{
constexpr int kMaxBins = 8;

struct CompareSeries
{
  const char *observable;
  const char *rapidity;
  const char *state;
  const char *stateLatex;
  const char *kinematicLatex;
  const char *xTitle;
  const char *outName;
  int n;
  double x[kMaxBins];
  double xBak[kMaxBins];
  double xGwak[kMaxBins];
  double bak[kMaxBins];
  double bakErr[kMaxBins];
  double gwak[kMaxBins];
  double gwakErr[kMaxBins];
  double ratio[kMaxBins];
  double ratioErr[kMaxBins];
  std::string binLabel[kMaxBins];
  double xMin;
  double xMax;
};

void setPlotStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
}

bool readBakHist(const char *fileName, const char *histName, CompareSeries &s)
{
  TFile file(fileName);
  if (file.IsZombie())
  {
    std::cerr << "[compare_Bak_Gwak_Jpsi] Cannot open " << fileName << std::endl;
    return false;
  }

  TH1 *h = dynamic_cast<TH1 *>(file.Get(histName));
  if (!h)
  {
    std::cerr << "[compare_Bak_Gwak_Jpsi] Missing " << histName << " in " << fileName << std::endl;
    return false;
  }
  if (h->GetNbinsX() < s.n)
  {
    std::cerr << "[compare_Bak_Gwak_Jpsi] " << histName << " has only " << h->GetNbinsX()
              << " bins, expected " << s.n << std::endl;
    return false;
  }

  for (int i = 0; i < s.n; ++i)
  {
    s.bak[i] = h->GetBinContent(i + 1);
    s.bakErr[i] = h->GetBinError(i + 1);
  }
  return true;
}

void fillRatio(CompareSeries &s)
{
  for (int i = 0; i < s.n; ++i)
  {
    s.ratio[i] = (s.gwak[i] != 0.0) ? s.bak[i] / s.gwak[i] : 0.0;
    const double bakRel = (s.bak[i] != 0.0) ? s.bakErr[i] / s.bak[i] : 0.0;
    const double gwakRel = (s.gwak[i] != 0.0) ? s.gwakErr[i] / s.gwak[i] : 0.0;
    s.ratioErr[i] = s.ratio[i] * std::sqrt(bakRel * bakRel + gwakRel * gwakRel);
  }
}

TGraphErrors *makeGraph(int n, const double *x, const double *y, const double *yErr)
{
  double xErr[kMaxBins] = {0.0};
  return new TGraphErrors(n, x, y, xErr, yErr);
}

void styleGraph(TGraphErrors *g, Color_t color, Style_t markerStyle, bool openMarker)
{
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(openMarker ? 1.45 : 1.35);
  g->SetLineWidth(2);
}

void drawCmsText()
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.96, 0.945, "PbPb 5.02 TeV");
  latex.SetTextAlign(11);
  latex.SetTextFont(72);
  latex.SetTextSize(0.046);
  latex.DrawLatex(0.15, 0.945, "CMS");
  latex.SetTextFont(42);
  latex.SetTextSize(0.037);
  latex.DrawLatex(0.245, 0.945, "Internal");
}

void drawSeries(const CompareSeries &s)
{
  TCanvas *c = new TCanvas(s.outName, "", 820, 860);
  c->cd();

  TPad *top = new TPad(Form("top_%s", s.outName), "", 0.0, 0.31, 1.0, 1.0);
  TPad *bottom = new TPad(Form("bottom_%s", s.outName), "", 0.0, 0.0, 1.0, 0.31);
  top->SetTopMargin(0.085);
  top->SetBottomMargin(0.02);
  top->SetLeftMargin(0.13);
  top->SetRightMargin(0.04);
  bottom->SetTopMargin(0.03);
  bottom->SetBottomMargin(0.33);
  bottom->SetLeftMargin(0.13);
  bottom->SetRightMargin(0.04);
  top->Draw();
  bottom->Draw();

  TGraphErrors *gBak = makeGraph(s.n, s.xBak, s.bak, s.bakErr);
  TGraphErrors *gGwak = makeGraph(s.n, s.xGwak, s.gwak, s.gwakErr);
  TGraphErrors *gRatio = makeGraph(s.n, s.x, s.ratio, s.ratioErr);
  styleGraph(gBak, kAzure + 2, 20, false);
  styleGraph(gGwak, kOrange + 7, 25, true);
  styleGraph(gRatio, kBlack, 20, false);
  gRatio->SetMarkerSize(0.95);

  top->cd();
  TH1D *frame = new TH1D(Form("frame_%s", s.outName), "", 100, s.xMin, s.xMax);
  frame->SetDirectory(0);
  frame->SetMinimum(0.0);
  frame->SetMaximum(1.08);
  frame->GetYaxis()->SetTitle("R_{AA}");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleSize(0.057);
  frame->GetYaxis()->SetTitleOffset(1.02);
  frame->GetYaxis()->SetLabelSize(0.044);
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->Draw("AXIS");

  TLine *lineOne = new TLine(s.xMin, 1.0, s.xMax, 1.0);
  lineOne->SetLineStyle(7);
  lineOne->SetLineColor(kGray + 2);
  lineOne->Draw();

  gBak->Draw("P SAME");
  gGwak->Draw("P SAME");

  TLegend *leg = new TLegend(0.60, 0.74, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.033);
  leg->AddEntry(gBak, "Bak (260502)", "pe");
  leg->AddEntry(gGwak, "Gwak (TnPL2L3)", "pe");
  leg->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.038);
  latex.DrawLatex(0.17, 0.84, s.stateLatex);
  latex.SetTextSize(0.034);
  latex.DrawLatex(0.17, 0.79, s.kinematicLatex);
  drawCmsText();

  bottom->cd();
  TH1D *ratioFrame = new TH1D(Form("ratioFrame_%s", s.outName), "", 100, s.xMin, s.xMax);
  ratioFrame->SetDirectory(0);
  ratioFrame->SetMinimum(0.92);
  ratioFrame->SetMaximum(1.12);
  ratioFrame->GetXaxis()->SetTitle(s.xTitle);
  ratioFrame->GetYaxis()->SetTitle("Bak/Gwak");
  ratioFrame->GetXaxis()->CenterTitle();
  ratioFrame->GetYaxis()->CenterTitle();
  ratioFrame->GetXaxis()->SetTitleSize(0.12);
  ratioFrame->GetYaxis()->SetTitleSize(0.095);
  ratioFrame->GetXaxis()->SetLabelSize(0.095);
  ratioFrame->GetYaxis()->SetLabelSize(0.082);
  ratioFrame->GetYaxis()->SetTitleOffset(0.56);
  ratioFrame->GetYaxis()->SetNdivisions(505);
  ratioFrame->Draw("AXIS");

  TLine *ratioOne = new TLine(s.xMin, 1.0, s.xMax, 1.0);
  ratioOne->SetLineStyle(7);
  ratioOne->SetLineColor(kGray + 2);
  ratioOne->Draw();
  gRatio->Draw("P SAME");

  c->SaveAs(Form("figs/%s.pdf", s.outName));
  c->SaveAs(Form("figs/%s.png", s.outName));
}

void writeCsvHeader(std::ofstream &csv)
{
  csv << "observable,rapidity,state,bin_label,x,bak,bak_stat,gwak,gwak_stat,bak_over_gwak,bak_over_gwak_stat\n";
}

void writeCsvRows(const CompareSeries &s, std::ofstream &csv)
{
  for (int i = 0; i < s.n; ++i)
  {
    csv << s.observable << ","
        << s.rapidity << ","
        << s.state << ","
        << s.binLabel[i] << ","
        << s.x[i] << ","
        << s.bak[i] << ","
        << s.bakErr[i] << ","
        << s.gwak[i] << ","
        << s.gwakErr[i] << ","
        << s.ratio[i] << ","
        << s.ratioErr[i] << "\n";
  }
}

void printSummary(const CompareSeries &s)
{
  std::cout << "\n[" << s.observable << ", " << s.rapidity << ", " << s.state << "]" << std::endl;
  for (int i = 0; i < s.n; ++i)
  {
    std::cout << "  " << s.binLabel[i]
              << "  Bak=" << s.bak[i] << " +/- " << s.bakErr[i]
              << "  Gwak=" << s.gwak[i] << " +/- " << s.gwakErr[i]
              << "  Bak/Gwak=" << s.ratio[i] << " +/- " << s.ratioErr[i]
              << std::endl;
  }
}

bool fillPtSeries(CompareSeries &s,
                  const char *rootFile,
                  const char *histName,
                  const double *ptEdges,
                  const double *gwakValue,
                  const double *gwakStat,
                  int nBins,
                  const char *rapidity,
                  const char *state,
                  const char *stateLatex,
                  const char *kinematicLatex,
                  const char *outName)
{
  s.observable = "pT";
  s.rapidity = rapidity;
  s.state = state;
  s.stateLatex = stateLatex;
  s.kinematicLatex = kinematicLatex;
  s.xTitle = "p_{T} (GeV/c)";
  s.outName = outName;
  s.n = nBins;
  s.xMin = 0.0;
  s.xMax = 40.0;

  for (int i = 0; i < nBins; ++i)
  {
    const double halfWidth = jpsi_raa::pt_half_width(ptEdges, i);
    const double shift = 0.10 * halfWidth;
    s.x[i] = jpsi_raa::pt_center(ptEdges, i);
    s.xBak[i] = s.x[i] - shift;
    s.xGwak[i] = s.x[i] + shift;
    s.gwak[i] = gwakValue[i];
    s.gwakErr[i] = gwakStat[i];
    s.binLabel[i] = Form("%.1f-%.1f", ptEdges[i], ptEdges[i + 1]);
  }

  if (!readBakHist(rootFile, histName, s))
    return false;
  fillRatio(s);
  return true;
}

bool fillCentSeries(CompareSeries &s,
                    const char *rootFile,
                    const char *histName,
                    const double *npart,
                    const double *gwakValueCentralToPeripheral,
                    const double *gwakStatCentralToPeripheral,
                    int nBins,
                    const char **centLabelsPeripheralToCentral,
                    const char *rapidity,
                    const char *state,
                    const char *stateLatex,
                    const char *kinematicLatex,
                    const char *outName)
{
  s.observable = "Npart";
  s.rapidity = rapidity;
  s.state = state;
  s.stateLatex = stateLatex;
  s.kinematicLatex = kinematicLatex;
  s.xTitle = "<N_{part}>";
  s.outName = outName;
  s.n = nBins;
  s.xMin = 0.0;
  s.xMax = 400.0;

  for (int i = 0; i < nBins; ++i)
  {
    s.x[i] = npart[i];
    s.xBak[i] = npart[i] - 4.0;
    s.xGwak[i] = npart[i] + 4.0;
    s.gwak[i] = gwakValueCentralToPeripheral[nBins - 1 - i];
    s.gwakErr[i] = gwakStatCentralToPeripheral[nBins - 1 - i];
    s.binLabel[i] = centLabelsPeripheralToCentral[i];
  }

  // The 260502 centrality ROOT histograms were written in peripheral-to-central value order.
  if (!readBakHist(rootFile, histName, s))
    return false;
  fillRatio(s);
  return true;
}

} // namespace

void compare_Bak_Gwak_Jpsi()
{
  setPlotStyle();
  gSystem->mkdir("figs", true);

  const char *centMidLabels[jpsi_raa::kNCentMid] = {"50-90", "40-50", "30-40", "20-30", "10-20", "0-10"};
  const char *centFwdLabels[jpsi_raa::kNCentFwd] = {"50-90", "30-50", "10-30", "0-10"};

  CompareSeries series[8];
  int nSeries = 0;

  if (fillPtSeries(series[nSeries],
                   "roots/RAA_JPsi_midRap_pT.root", "hRAA_PR",
                   jpsi_raa::kPtMidBinEdges, jpsi_raa::kTnPL2L3PtMidPr, jpsi_raa::kTnPL2L3PtMidPrStat, jpsi_raa::kNPtMid,
                   "mid", "PR", "Prompt J/#psi", "Cent. 0-90%, |y| < 1.6",
                   "compare_Bak_Gwak_mid_pT_Jpsi_PR"))
    ++nSeries;

  if (fillPtSeries(series[nSeries],
                   "roots/RAA_JPsi_midRap_pT.root", "hRAA_NP",
                   jpsi_raa::kPtMidBinEdges, jpsi_raa::kTnPL2L3PtMidNp, jpsi_raa::kTnPL2L3PtMidNpStat, jpsi_raa::kNPtMid,
                   "mid", "NP", "Nonprompt J/#psi", "Cent. 0-90%, |y| < 1.6",
                   "compare_Bak_Gwak_mid_pT_Jpsi_NP"))
    ++nSeries;

  if (fillPtSeries(series[nSeries],
                   "roots/RAA_JPsi_forRap_pT.root", "hRAA_PR",
                   jpsi_raa::kPtFwdBinEdges, jpsi_raa::kTnPL2L3PtFwdPr, jpsi_raa::kTnPL2L3PtFwdPrStat, jpsi_raa::kNPtFwd,
                   "fwd", "PR", "Prompt J/#psi", "Cent. 0-90%, 1.6 < |y| < 2.4",
                   "compare_Bak_Gwak_fwd_pT_Jpsi_PR"))
    ++nSeries;

  if (fillPtSeries(series[nSeries],
                   "roots/RAA_JPsi_forRap_pT.root", "hRAA_NP",
                   jpsi_raa::kPtFwdBinEdges, jpsi_raa::kTnPL2L3PtFwdNp, jpsi_raa::kTnPL2L3PtFwdNpStat, jpsi_raa::kNPtFwd,
                   "fwd", "NP", "Nonprompt J/#psi", "Cent. 0-90%, 1.6 < |y| < 2.4",
                   "compare_Bak_Gwak_fwd_pT_Jpsi_NP"))
    ++nSeries;

  if (fillCentSeries(series[nSeries],
                     "roots/RAA_JPsi_midRap_Npart.root", "hRAA_PR",
                     jpsi_raa::kMidNpart, jpsi_raa::kTnPL2L3CentMidPr, jpsi_raa::kTnPL2L3CentMidPrStat,
                     jpsi_raa::kNCentMid, centMidLabels,
                     "mid", "PR", "Prompt J/#psi", "6.5 < p_{T} < 40 GeV/c, |y| < 1.6",
                     "compare_Bak_Gwak_mid_Npart_Jpsi_PR"))
    ++nSeries;

  if (fillCentSeries(series[nSeries],
                     "roots/RAA_JPsi_midRap_Npart.root", "hRAA_NP",
                     jpsi_raa::kMidNpart, jpsi_raa::kTnPL2L3CentMidNp, jpsi_raa::kTnPL2L3CentMidNpStat,
                     jpsi_raa::kNCentMid, centMidLabels,
                     "mid", "NP", "Nonprompt J/#psi", "6.5 < p_{T} < 40 GeV/c, |y| < 1.6",
                     "compare_Bak_Gwak_mid_Npart_Jpsi_NP"))
    ++nSeries;

  if (fillCentSeries(series[nSeries],
                     "roots/RAA_JPsi_forRap_Npart_4Bins.root", "hRAA_PR",
                     jpsi_raa::kFwdNpart, jpsi_raa::kTnPL2L3CentFwdPr, jpsi_raa::kTnPL2L3CentFwdPrStat,
                     jpsi_raa::kNCentFwd, centFwdLabels,
                     "fwd", "PR", "Prompt J/#psi", "3.5 < p_{T} < 40 GeV/c, 1.6 < |y| < 2.4",
                     "compare_Bak_Gwak_fwd_Npart_Jpsi_PR"))
    ++nSeries;

  if (fillCentSeries(series[nSeries],
                     "roots/RAA_JPsi_forRap_Npart_4Bins.root", "hRAA_NP",
                     jpsi_raa::kFwdNpart, jpsi_raa::kTnPL2L3CentFwdNp, jpsi_raa::kTnPL2L3CentFwdNpStat,
                     jpsi_raa::kNCentFwd, centFwdLabels,
                     "fwd", "NP", "Nonprompt J/#psi", "3.5 < p_{T} < 40 GeV/c, 1.6 < |y| < 2.4",
                     "compare_Bak_Gwak_fwd_Npart_Jpsi_NP"))
    ++nSeries;

  std::ofstream csv("Bak_Gwak_RAA_comparison_260502.csv");
  writeCsvHeader(csv);

  for (int i = 0; i < nSeries; ++i)
  {
    drawSeries(series[i]);
    writeCsvRows(series[i], csv);
    printSummary(series[i]);
  }

  csv.close();
  std::cout << "\nSaved CSV: Bak_Gwak_RAA_comparison_260502.csv" << std::endl;
}
