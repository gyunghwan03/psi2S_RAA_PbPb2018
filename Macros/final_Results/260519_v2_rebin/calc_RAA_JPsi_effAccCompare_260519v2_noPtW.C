#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
#include "../../../cutsAndBin.h"
#include "JPsiEffAccInputs_260519v2.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace JPsiRAAEffAccCompare260519v2NoPtW
{
const double kNmb = 11968044281.;
const double kLumiPP = 3.002;
const double kLumiPPScale = 1e-9;

TString NoPtWDir() { return "/data/hwan/psi2S_RAA_PbPb2018/Eff_Acc_260519_v2/roots/noPtW_full"; }
TString NoPtWTag() { return "260518full_noPtW"; }

TFile *OpenNoPtW(const char *fileName)
{
  TFile *f = TFile::Open(Form("%s/%s", NoPtWDir().Data(), fileName), "READ");
  if (!f || f->IsZombie())
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] cannot open " << fileName << std::endl;
  return f;
}

struct FitInput {
  double val = 0.0;
  double err = 0.0;
};

struct YieldPair {
  double pr = 0.0;
  double np = 0.0;
  double prErr = 0.0;
  double npErr = 0.0;
};

struct BinDef {
  double ptLow = 0.0;
  double ptHigh = 0.0;
  double yLow = 0.0;
  double yHigh = 0.0;
  double centLow = 0.0;
  double centHigh = 180.0;
  double taa = 6.274;
};

struct AxisDef {
  TString name;
  TString axisTitle;
  TString plotAxisTitle;
  bool isPt = true;
  bool isMid = true;
  std::vector<double> edges;
  std::vector<double> plotX;
  std::vector<BinDef> bins;
};

struct EffAccSet {
  std::unique_ptr<TFile> fEffPbPbPR;
  std::unique_ptr<TFile> fEffPbPbNP;
  std::unique_ptr<TFile> fEffPPPR;
  std::unique_ptr<TFile> fEffPPNP;
  std::unique_ptr<TFile> fAccPbPbPR;
  std::unique_ptr<TFile> fAccPbPbNP;
  std::unique_ptr<TFile> fAccPPPR;
  std::unique_ptr<TFile> fAccPPNP;

  TH1D *hEffPbPbPR = nullptr;
  TH1D *hEffPbPbNP = nullptr;
  TH1D *hEffPPPR = nullptr;
  TH1D *hEffPPNP = nullptr;
  TH1D *hAccPbPbPR = nullptr;
  TH1D *hAccPbPbNP = nullptr;
  TH1D *hAccPPPR = nullptr;
  TH1D *hAccPPNP = nullptr;

  std::unique_ptr<TH1D> hEffPPPRInt;
  std::unique_ptr<TH1D> hEffPPNPInt;
};

TFile *OpenFitFile(const TString &path)
{
  TFile *file = TFile::Open(path, "READ");
  if (!file || file->IsZombie())
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] failed to open " << path << std::endl;
  return file;
}

FitInput ReadFitValue(const TString &path, const char *histName, int bin)
{
  std::unique_ptr<TFile> file(OpenFitFile(path));
  FitInput out;
  if (!file || file->IsZombie())
    return out;

  TH1D *hist = dynamic_cast<TH1D *>(file->Get(histName));
  if (!hist)
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] missing " << histName
              << " in " << path << std::endl;
    return out;
  }

  out.val = hist->GetBinContent(bin);
  out.err = hist->GetBinError(bin);
  return out;
}

TString PbPbLabel(const BinDef &bin)
{
  return getKineLabel(bin.ptLow, bin.ptHigh, bin.yLow, bin.yHigh, 0.0, bin.centLow, bin.centHigh);
}

TString PPLabel(const BinDef &bin)
{
  return getKineLabelpp(bin.ptLow, bin.ptHigh, bin.yLow, bin.yHigh, 0.0);
}

FitInput GetYieldPbPb(const BinDef &bin)
{
  const TString label = PbPbLabel(bin);
  return ReadFitValue(Form("../../Jpsi/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                           label.Data()),
                      "fitResults", 1);
}

FitInput GetYieldPP(const BinDef &bin)
{
  const TString label = PPLabel(bin);
  return ReadFitValue(Form("../../pp_Jpsi/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                           label.Data()),
                      "fitResults", 1);
}

FitInput GetFracPbPb(const BinDef &bin)
{
  const TString label = PbPbLabel(bin);
  return ReadFitValue(Form("../../Jpsi/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                           label.Data()),
                      "2DfitResults", 1);
}

FitInput GetFracPP(const BinDef &bin)
{
  const TString label = PPLabel(bin);
  return ReadFitValue(Form("../../pp_Jpsi/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                           label.Data()),
                      "2DfitResults", 1);
}

double SafeRelErr(double err, double val)
{
  return (val != 0.0) ? err / val : 0.0;
}

YieldPair SplitPromptNonPrompt(const FitInput &yield, const FitInput &frac, double wPR, double wNP)
{
  YieldPair out;
  if (wPR <= 0.0 || wNP <= 0.0)
    return out;

  const double rawPR = yield.val * (1.0 - frac.val);
  const double rawNP = yield.val * frac.val;
  const double rawPRErr = rawPR * std::sqrt(std::pow(SafeRelErr(yield.err, yield.val), 2) +
                                            std::pow(SafeRelErr(frac.err, 1.0 - frac.val), 2));
  const double rawNPErr = rawNP * std::sqrt(std::pow(SafeRelErr(yield.err, yield.val), 2) +
                                            std::pow(SafeRelErr(frac.err, frac.val), 2));

  out.pr = rawPR / wPR;
  out.np = rawNP / wNP;
  out.prErr = rawPRErr / wPR;
  out.npErr = rawNPErr / wNP;
  return out;
}

double HistValue(TH1D *hist, int bin, const char *tag, bool enabled)
{
  if (!enabled)
    return 1.0;
  if (!hist)
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] null histogram for " << tag << std::endl;
    return 1.0;
  }
  const double value = hist->GetBinContent(bin);
  if (value <= 0.0)
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] non-positive " << tag
              << " bin " << bin << ": " << value << std::endl;
    return 1.0;
  }
  return value;
}

std::unique_ptr<TFile> TakeFile(TFile *file)
{
  return std::unique_ptr<TFile>(file);
}

EffAccSet LoadEffAcc(const AxisDef &axis)
{
  EffAccSet out;
  const TString tag = NoPtWTag();
  out.fEffPbPbPR = TakeFile(OpenNoPtW(Form("mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW1_tnp1_%s_noPtW.root", tag.Data())));
  out.fEffPbPbNP = TakeFile(OpenNoPtW(Form("mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW1_tnp1_%s_noPtW.root", tag.Data())));
  out.fEffPPPR = TakeFile(OpenNoPtW(Form("mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW1_tnp1_%s_noPtW.root", tag.Data())));
  out.fEffPPNP = TakeFile(OpenNoPtW(Form("mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtW1_tnp1_%s_noPtW.root", tag.Data())));
  out.fAccPPPR = TakeFile(OpenNoPtW(Form("acceptance_PromptJpsi_GenOnly_wgt0_pp_SysUp0_%s.root", tag.Data())));
  out.fAccPbPbPR = TakeFile(OpenNoPtW(Form("acceptance_PromptJpsi_GenOnly_wgt0_PbPb_SysUp0_%s.root", tag.Data())));
  out.fAccPPNP = TakeFile(OpenNoPtW(Form("acceptance_BtoJpsi_GenOnly_wgt0_pp_SysUp0_%s.root", tag.Data())));
  out.fAccPbPbNP = TakeFile(OpenNoPtW(Form("acceptance_BtoJpsi_GenOnly_wgt0_PbPb_SysUp0_%s.root", tag.Data())));

  const TString rap = axis.isMid ? "mid" : "fwd";
  const TString absY = axis.isMid ? "absy0_1p6" : "absy1p6_2p4";
  const TString accRap = axis.isMid ? "midy" : "Fory";

  if (axis.isPt)
  {
    out.hEffPbPbPR = dynamic_cast<TH1D *>(out.fEffPbPbPR->Get(Form("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_%s", absY.Data())));
    out.hEffPbPbNP = dynamic_cast<TH1D *>(out.fEffPbPbNP->Get(Form("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_%s", absY.Data())));
    out.hEffPPPR = dynamic_cast<TH1D *>(out.fEffPPPR->Get(Form("mc_eff_vs_pt_TnP1_PtW1_%s", absY.Data())));
    out.hEffPPNP = dynamic_cast<TH1D *>(out.fEffPPNP->Get(Form("mc_eff_vs_pt_TnP1_PtW1_%s", absY.Data())));
    out.hAccPbPbPR = dynamic_cast<TH1D *>(out.fAccPbPbPR->Get(Form("hAccPt_2021_%s", accRap.Data())));
    out.hAccPbPbNP = dynamic_cast<TH1D *>(out.fAccPbPbNP->Get(Form("hAccPt_2021_%s", accRap.Data())));
    out.hAccPPPR = dynamic_cast<TH1D *>(out.fAccPPPR->Get(Form("hAccPt_2021_%s", accRap.Data())));
    out.hAccPPNP = dynamic_cast<TH1D *>(out.fAccPPNP->Get(Form("hAccPt_2021_%s", accRap.Data())));
  }
  else
  {
    const TString ptTag = axis.isMid ? "pt_6p5_to_40" : "pt_3_to_40";
    out.hEffPbPbPR = dynamic_cast<TH1D *>(out.fEffPbPbPR->Get(Form("mc_eff_vs_cent_TnP1_PtW1_%s_%s", ptTag.Data(), absY.Data())));
    out.hEffPbPbNP = dynamic_cast<TH1D *>(out.fEffPbPbNP->Get(Form("mc_eff_vs_cent_TnP1_PtW1_%s_%s", ptTag.Data(), absY.Data())));
    out.hEffPPPRInt.reset(JPsiEffAcc260519v2::BuildIntegratedRatioHist(out.fEffPPPR.get(), Form("hist_eff_num_%s", rap.Data()), Form("hist_eff_den_%s", rap.Data()), Form("hEff_ppPR_%s_int", rap.Data())));
    out.hEffPPNPInt.reset(JPsiEffAcc260519v2::BuildIntegratedRatioHist(out.fEffPPNP.get(), Form("hist_eff_num_%s", rap.Data()), Form("hist_eff_den_%s", rap.Data()), Form("hEff_ppNP_%s_int", rap.Data())));
    out.hEffPPPR = out.hEffPPPRInt.get();
    out.hEffPPNP = out.hEffPPNPInt.get();
    out.hAccPbPbPR = dynamic_cast<TH1D *>(out.fAccPbPbPR->Get(Form("hAccPt_2021_%s_Int", accRap.Data())));
    out.hAccPbPbNP = dynamic_cast<TH1D *>(out.fAccPbPbNP->Get(Form("hAccPt_2021_%s_Int", accRap.Data())));
    out.hAccPPPR = dynamic_cast<TH1D *>(out.fAccPPPR->Get(Form("hAccPt_2021_%s_Int", accRap.Data())));
    out.hAccPPNP = dynamic_cast<TH1D *>(out.fAccPPNP->Get(Form("hAccPt_2021_%s_Int", accRap.Data())));
  }

  return out;
}

std::unique_ptr<TH1D> MakeHist(const AxisDef &axis, const TString &name, const TString &title)
{
  std::unique_ptr<TH1D> hist(new TH1D(name, Form("%s;%s;R_{AA}", title.Data(), axis.axisTitle.Data()),
                                      axis.edges.size() - 1, axis.edges.data()));
  hist->SetDirectory(nullptr);
  return hist;
}

void FillRAA(const AxisDef &axis, bool useEff, bool useAcc, TDirectory *outDir)
{
  EffAccSet effAcc = LoadEffAcc(axis);
  const TString tag = Form("%s_%s", useEff ? "eff" : "noEff", useAcc ? "acc" : "noAcc");

  std::unique_ptr<TH1D> hRAAPR = MakeHist(axis, Form("hRAA_PR_%s", tag.Data()), axis.name);
  std::unique_ptr<TH1D> hRAANP = MakeHist(axis, Form("hRAA_NP_%s", tag.Data()), axis.name);
  std::unique_ptr<TH1D> hPbPbPR = MakeHist(axis, Form("hYieldPbPb_PR_%s", tag.Data()), axis.name);
  std::unique_ptr<TH1D> hPbPbNP = MakeHist(axis, Form("hYieldPbPb_NP_%s", tag.Data()), axis.name);
  std::unique_ptr<TH1D> hPPPR = MakeHist(axis, Form("hYieldPP_PR_%s", tag.Data()), axis.name);
  std::unique_ptr<TH1D> hPPNP = MakeHist(axis, Form("hYieldPP_NP_%s", tag.Data()), axis.name);

  std::cout << "\n=== " << axis.name << " useEff=" << useEff << " useAcc=" << useAcc << " ===" << std::endl;
  for (int i = 0; i < static_cast<int>(axis.bins.size()); ++i)
  {
    const BinDef &bin = axis.bins[i];
    const int hbin = i + 1;
    const int eaBin = axis.isPt ? hbin : (axis.isMid ? hbin : hbin);
    const int intBin = axis.isPt ? hbin : 1;

    const double wPbPbPR = HistValue(effAcc.hEffPbPbPR, eaBin, "PbPb PR efficiency", useEff) *
                           HistValue(effAcc.hAccPbPbPR, intBin, "PbPb PR acceptance", useAcc);
    const double wPbPbNP = HistValue(effAcc.hEffPbPbNP, eaBin, "PbPb NP efficiency", useEff) *
                           HistValue(effAcc.hAccPbPbNP, intBin, "PbPb NP acceptance", useAcc);
    const double wPPPR = HistValue(effAcc.hEffPPPR, intBin, "pp PR efficiency", useEff) *
                         HistValue(effAcc.hAccPPPR, intBin, "pp PR acceptance", useAcc);
    const double wPPNP = HistValue(effAcc.hEffPPNP, intBin, "pp NP efficiency", useEff) *
                         HistValue(effAcc.hAccPPNP, intBin, "pp NP acceptance", useAcc);

    const YieldPair pbpb = SplitPromptNonPrompt(GetYieldPbPb(bin), GetFracPbPb(bin), wPbPbPR, wPbPbNP);
    const YieldPair pp = SplitPromptNonPrompt(GetYieldPP(bin), GetFracPP(bin), wPPPR, wPPNP);

    const double ptWidth = bin.ptHigh - bin.ptLow;
    const double yWidth = bin.yHigh - bin.yLow;
    const double centFrac = axis.isPt ? 1.0 : ((bin.centHigh - bin.centLow) / 180.0);

    const double xPPPR = kLumiPPScale * pp.pr / (kLumiPP * 1e2 * ptWidth * 2.0 * yWidth);
    const double xPPNP = kLumiPPScale * pp.np / (kLumiPP * 1e2 * ptWidth * 2.0 * yWidth);
    const double xPbPbPR = pbpb.pr / (kNmb * bin.taa * ptWidth * 2.0 * yWidth * centFrac);
    const double xPbPbNP = pbpb.np / (kNmb * bin.taa * ptWidth * 2.0 * yWidth * centFrac);

    const double raaPR = (xPPPR > 0.0) ? xPbPbPR / xPPPR : 0.0;
    const double raaNP = (xPPNP > 0.0) ? xPbPbNP / xPPNP : 0.0;
    const double raaPRErr = raaPR * std::sqrt(std::pow(SafeRelErr(pbpb.prErr, pbpb.pr), 2) +
                                              std::pow(SafeRelErr(pp.prErr, pp.pr), 2));
    const double raaNPErr = raaNP * std::sqrt(std::pow(SafeRelErr(pbpb.npErr, pbpb.np), 2) +
                                              std::pow(SafeRelErr(pp.npErr, pp.np), 2));

    hRAAPR->SetBinContent(hbin, raaPR);
    hRAAPR->SetBinError(hbin, raaPRErr);
    hRAANP->SetBinContent(hbin, raaNP);
    hRAANP->SetBinError(hbin, raaNPErr);
    hPbPbPR->SetBinContent(hbin, pbpb.pr);
    hPbPbPR->SetBinError(hbin, pbpb.prErr);
    hPbPbNP->SetBinContent(hbin, pbpb.np);
    hPbPbNP->SetBinError(hbin, pbpb.npErr);
    hPPPR->SetBinContent(hbin, pp.pr);
    hPPPR->SetBinError(hbin, pp.prErr);
    hPPNP->SetBinContent(hbin, pp.np);
    hPPNP->SetBinError(hbin, pp.npErr);

    std::cout << "bin " << hbin
              << " RAA_PR=" << raaPR << " +/- " << raaPRErr
              << " RAA_NP=" << raaNP << " +/- " << raaNPErr
              << " weights(PbPbPR,PbPbNP,ppPR,ppNP)="
              << wPbPbPR << "," << wPbPbNP << "," << wPPPR << "," << wPPNP
              << std::endl;
  }

  outDir->cd();
  TDirectory *caseDir = outDir->mkdir(tag);
  caseDir->cd();
  hRAAPR->Write();
  hRAANP->Write();
  hPbPbPR->Write();
  hPbPbNP->Write();
  hPPPR->Write();
  hPPNP->Write();
  outDir->cd();
}

AxisDef MidPt()
{
  AxisDef a;
  a.name = "midRap_pT";
  a.axisTitle = "p_{T} (GeV/c)";
  a.isPt = true;
  a.isMid = true;
  a.edges = {6.5, 9, 12, 15, 20, 25, 40};
  for (int i = 0; i < 6; ++i)
    a.bins.push_back({a.edges[i], a.edges[i + 1], 0.0, 1.6, 0.0, 180.0, 6.274});
  return a;
}

AxisDef FwdPt()
{
  AxisDef a;
  a.name = "fwdRap_pT";
  a.axisTitle = "p_{T} (GeV/c)";
  a.isPt = true;
  a.isMid = false;
  a.edges = {3.5, 6.5, 9, 12, 40};
  for (int i = 0; i < 4; ++i)
    a.bins.push_back({a.edges[i], a.edges[i + 1], 1.6, 2.4, 0.0, 180.0, 6.274});
  return a;
}

AxisDef MidCent()
{
  AxisDef a;
  a.name = "midRap_Npart";
  a.axisTitle = "Centrality (%)";
  a.plotAxisTitle = "<N_{Part}>";
  a.isPt = false;
  a.isMid = true;
  a.edges = {0, 10, 20, 30, 40, 50, 90};
  a.plotX = {356.9, 262.3, 188.2, 131.0, 87.19, 27.12};
  const double taa[6] = {23.05, 14.39, 8.798, 5.124, 2.777, 0.5803};
  const double centHi[7] = {0, 20, 40, 60, 80, 100, 180};
  for (int i = 0; i < 6; ++i)
    a.bins.push_back({6.5, 40.0, 0.0, 1.6, centHi[i], centHi[i + 1], taa[i]});
  return a;
}

AxisDef FwdCent()
{
  AxisDef a;
  a.name = "fwdRap_Npart_4Bins";
  a.axisTitle = "Centrality (%)";
  a.plotAxisTitle = "<N_{Part}>";
  a.isPt = false;
  a.isMid = false;
  a.edges = {0, 10, 30, 50, 90};
  a.plotX = {356.9, 225.2, 109.1, 27.12};
  const double taa[4] = {23.05, 11.60, 3.950, 0.5803};
  const double centHi[5] = {0, 20, 60, 100, 180};
  for (int i = 0; i < 4; ++i)
    a.bins.push_back({3.5, 40.0, 1.6, 2.4, centHi[i], centHi[i + 1], taa[i]});
  return a;
}

void RunAxis(const AxisDef &axis, TFile &outFile)
{
  TDirectory *axisDir = outFile.mkdir(axis.name);
  FillRAA(axis, true, true, axisDir);
  FillRAA(axis, false, false, axisDir);
  FillRAA(axis, true, false, axisDir);
  FillRAA(axis, false, true, axisDir);
}

void StyleRAAHist(TH1D *hist, int color, int marker, int lineStyle)
{
  if (!hist)
    return;
  hist->SetStats(0);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.1);
  hist->SetLineStyle(lineStyle);
  hist->SetLineWidth(2);
  hist->GetYaxis()->SetRangeUser(0.0, 1.5);
}

TString PlotAxisTitle(const AxisDef &axis)
{
  return axis.plotAxisTitle.Length() ? axis.plotAxisTitle : axis.axisTitle;
}

double PlotXMin(const AxisDef &axis)
{
  return axis.isPt ? axis.edges.front() : 0.0;
}

double PlotXMax(const AxisDef &axis)
{
  return axis.isPt ? axis.edges.back() : 400.0;
}

std::unique_ptr<TGraphErrors> MakeRAAGraph(const AxisDef &axis, TH1D *hist, const TString &name)
{
  if (!hist)
    return nullptr;

  const int n = hist->GetNbinsX();
  std::vector<double> x(n), y(n), ex(n), ey(n);
  for (int i = 0; i < n; ++i)
  {
    const int bin = i + 1;
    x[i] = (!axis.plotX.empty() && i < static_cast<int>(axis.plotX.size())) ? axis.plotX[i] : hist->GetBinCenter(bin);
    y[i] = hist->GetBinContent(bin);
    ex[i] = axis.isPt ? 0.5 * hist->GetBinWidth(bin) : 0.0;
    ey[i] = hist->GetBinError(bin);
  }

  std::unique_ptr<TGraphErrors> graph(new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data()));
  graph->SetName(name);
  return graph;
}

void StyleRAAGraph(TGraphErrors *graph, int color, int marker, int lineStyle)
{
  if (!graph)
    return;
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.1);
  graph->SetLineStyle(lineStyle);
  graph->SetLineWidth(2);
}

std::unique_ptr<TGraphErrors> MakeHINPtGraph(const char *fileName, const char *histName,
                                             const char *errName, int nPoints)
{
  std::unique_ptr<TFile> file(TFile::Open(fileName, "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] failed to open HIN file "
              << fileName << std::endl;
    return nullptr;
  }

  TH1 *hist = dynamic_cast<TH1 *>(file->Get(histName));
  TH1 *err = dynamic_cast<TH1 *>(file->Get(errName));
  if (!hist || !err)
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] missing HIN histogram "
              << histName << " or " << errName << " in " << fileName << std::endl;
    return nullptr;
  }

  std::vector<double> x(nPoints), y(nPoints), ex(nPoints), ey(nPoints);
  for (int i = 0; i < nPoints; ++i)
  {
    const int bin = i + 1;
    x[i] = hist->GetXaxis()->GetBinCenter(bin);
    y[i] = hist->GetBinContent(bin);
    ex[i] = 0.5 * hist->GetXaxis()->GetBinWidth(bin);
    ey[i] = err->GetBinContent(bin);
  }
  return std::unique_ptr<TGraphErrors>(new TGraphErrors(nPoints, x.data(), y.data(), ex.data(), ey.data()));
}

std::unique_ptr<TGraphErrors> MakeHINNpartGraph(const char *fileName, const char *histName,
                                                const char *errName, const std::vector<double> &npart)
{
  std::unique_ptr<TFile> file(TFile::Open(fileName, "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] failed to open HIN file "
              << fileName << std::endl;
    return nullptr;
  }

  TH1 *hist = dynamic_cast<TH1 *>(file->Get(histName));
  TH1 *err = dynamic_cast<TH1 *>(file->Get(errName));
  if (!hist || !err)
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] missing HIN histogram "
              << histName << " or " << errName << " in " << fileName << std::endl;
    return nullptr;
  }

  const int nPoints = npart.size();
  std::vector<double> x(nPoints), y(nPoints), ex(nPoints, 0.0), ey(nPoints);
  for (int i = 0; i < nPoints; ++i)
  {
    const int bin = nPoints - i;
    x[i] = npart[i];
    y[i] = hist->GetBinContent(bin);
    ey[i] = err->GetBinContent(bin);
  }
  return std::unique_ptr<TGraphErrors>(new TGraphErrors(nPoints, x.data(), y.data(), ex.data(), ey.data()));
}

void StyleHINGraph(TGraphErrors *graph, int marker)
{
  if (!graph)
    return;
  graph->SetMarkerColor(kGray + 2);
  graph->SetLineColor(kGray + 2);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.25);
}

std::vector<std::unique_ptr<TGraphErrors>> MakeHINGraphs(const AxisDef &axis, bool isPrompt)
{
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  const std::vector<double> npartMidHIN = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  const std::vector<double> npartFwdPRHIN = {32.7, 160.3, 311.5};
  const std::vector<double> npartFwdNPHIN = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};

  if (axis.isPt)
  {
    if (isPrompt && axis.isMid)
      graphs.push_back(MakeHINPtGraph("../roots/RAA_PR_Jpsi_HIN_16_025_mid_pT.root", "Table 18/Hist1D_y1", "Table 18/Hist1D_y1_e1", 5));
    else if (isPrompt)
      graphs.push_back(MakeHINPtGraph("../roots/RAA_PR_Jpsi_HIN_16_025_fwd_pT.root", "Table 22/Hist1D_y1", "Table 22/Hist1D_y1_e1", 3));
    else if (axis.isMid)
      graphs.push_back(MakeHINPtGraph("../roots/RAA_NP_Jpsi_HIN_16_025_mid_pT.root", "Table 28/Hist1D_y1", "Table 28/Hist1D_y1_e1", 6));
    else
      graphs.push_back(MakeHINPtGraph("../roots/RAA_NP_Jpsi_HIN_16_025_fwd_pT.root", "Table 29/Hist1D_y1", "Table 29/Hist1D_y1_e1", 10));
  }
  else
  {
    if (isPrompt && axis.isMid)
      graphs.push_back(MakeHINNpartGraph("../roots/RAA_PR_JPsi_HIN_16_025_mid_Npart.root", "Table 15/Hist1D_y1", "Table 15/Hist1D_y1_e1", npartMidHIN));
    else if (isPrompt)
      graphs.push_back(MakeHINNpartGraph("../roots/RAA_PR_JPsi_HIN_16_025_fwd_Npart.root", "Table 19/Hist1D_y1", "Table 19/Hist1D_y1_e1", npartFwdPRHIN));
    else if (axis.isMid)
      graphs.push_back(MakeHINNpartGraph("../roots/RAA_NP_JPsi_HIN_16_025_Npart.root", "Table 30/Hist1D_y1", "Table 30/Hist1D_y1_e1", npartMidHIN));
    else
    {
      graphs.push_back(MakeHINNpartGraph("../roots/RAA_NP_JPsi_HIN_16_025_Npart.root", "Table 30/Hist1D_y2", "Table 30/Hist1D_y2_e1", npartFwdNPHIN));
      graphs.push_back(MakeHINNpartGraph("../roots/RAA_NP_JPsi_HIN_16_025_fwd_Npart.root", "Table 32/Hist1D_y1", "Table 32/Hist1D_y1_e1", npartFwdNPHIN));
    }
  }

  for (int i = 0; i < static_cast<int>(graphs.size()); ++i)
    StyleHINGraph(graphs[i].get(), (i == 0) ? 25 : 27);
  return graphs;
}

void DrawRAAComponent(const AxisDef &axis, TH1D *hists[4], const char *caseLabels[4],
                      const char *componentTag, const char *componentLabel, bool isPrompt)
{
  if (!hists[0])
  {
    std::cerr << "[calc_RAA_JPsi_effAccCompare_260519v2_noPtW] missing " << componentLabel
              << " RAA histograms for " << axis.name << std::endl;
    return;
  }

  std::unique_ptr<TGraphErrors> graphs[4];
  for (int i = 0; i < 4; ++i)
  {
    graphs[i] = MakeRAAGraph(axis, hists[i], Form("gRAA_%s_%s_%d", axis.name.Data(), componentTag, i));
    if (graphs[i])
      StyleRAAGraph(graphs[i].get(), hists[i]->GetLineColor(), hists[i]->GetMarkerStyle(), hists[i]->GetLineStyle());
  }
  std::vector<std::unique_ptr<TGraphErrors>> hinGraphs = MakeHINGraphs(axis, isPrompt);

  TCanvas *canvas = new TCanvas(Form("cRAA_%s_%s", axis.name.Data(), componentTag), "", 900, 800);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.07);

  graphs[0]->SetTitle(Form("; %s;R_{AA}", PlotAxisTitle(axis).Data()));
  graphs[0]->GetXaxis()->SetLimits(PlotXMin(axis), PlotXMax(axis));
  graphs[0]->SetMinimum(0.0);
  graphs[0]->SetMaximum(1.5);
  graphs[0]->Draw("AP");
  for (int i = 1; i < 4; ++i)
    if (graphs[i])
      graphs[i]->Draw("P same");
  for (auto &hinGraph : hinGraphs)
    if (hinGraph)
      hinGraph->Draw("P same");

  TLine *line = new TLine(PlotXMin(axis), 1.0, PlotXMax(axis), 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->Draw("same");

  TLegend *legCases = new TLegend(0.18, 0.63, axis.isPt ? 0.54 : 0.58, 0.89);
  legCases->SetBorderSize(0);
  legCases->SetFillStyle(0);
  for (int i = 0; i < 4; ++i)
    if (graphs[i])
      legCases->AddEntry(graphs[i].get(), caseLabels[i], "pe");
  if (!hinGraphs.empty() && hinGraphs[0])
    legCases->AddEntry(hinGraphs[0].get(), "HIN-16-025", "pe");
  if (hinGraphs.size() > 1 && hinGraphs[1])
    legCases->AddEntry(hinGraphs[1].get(), "HIN-16-025 low p_{T}", "pe");
  legCases->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.038);
  latex.DrawLatex(0.18, 0.93, Form("J/#psi %s %s", axis.name.Data(), componentLabel));

  canvas->SaveAs(Form("./figs/RAA_JPsi_effAccCompare_260519v2_noPtW_%s_%s.pdf", axis.name.Data(), componentTag));
  canvas->SaveAs(Form("./figs/RAA_JPsi_effAccCompare_260519v2_noPtW_%s_%s.png", axis.name.Data(), componentTag));
}

void DrawAxisRAA(TFile &outFile, const AxisDef &axis)
{
  const char *caseTags[4] = {"eff_acc", "noEff_noAcc", "eff_noAcc", "noEff_acc"};
  const char *caseLabels[4] = {"Eff. #times Acc.", "No eff./acc.", "Eff. only", "Acc. only"};
  const int colors[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
  const int markersPR[4] = {20, 21, 22, 23};
  const int markersNP[4] = {20, 21, 22, 23};
  const int lineStyles[4] = {1, 2, 3, 4};

  TH1D *hPR[4] = {nullptr};
  TH1D *hNP[4] = {nullptr};
  for (int i = 0; i < 4; ++i)
  {
    hPR[i] = dynamic_cast<TH1D *>(outFile.Get(Form("%s/%s/hRAA_PR_%s", axis.name.Data(), caseTags[i], caseTags[i])));
    hNP[i] = dynamic_cast<TH1D *>(outFile.Get(Form("%s/%s/hRAA_NP_%s", axis.name.Data(), caseTags[i], caseTags[i])));
    StyleRAAHist(hPR[i], colors[i], markersPR[i], lineStyles[i]);
    StyleRAAHist(hNP[i], colors[i], markersNP[i], lineStyles[i]);
  }

  DrawRAAComponent(axis, hPR, caseLabels, "prompt", "Prompt", true);
  DrawRAAComponent(axis, hNP, caseLabels, "nonprompt", "Nonprompt", false);
}
}

void calc_RAA_JPsi_effAccCompare_260519v2_noPtW()
{
  using namespace JPsiRAAEffAccCompare260519v2NoPtW;

  gStyle->SetOptStat(0);
  gSystem->mkdir("./figs", true);
  TFile outFile("./roots/RAA_JPsi_effAccCompare_260519v2_noPtW.root", "RECREATE");
  RunAxis(MidPt(), outFile);
  RunAxis(FwdPt(), outFile);
  RunAxis(MidCent(), outFile);
  RunAxis(FwdCent(), outFile);
  DrawAxisRAA(outFile, MidPt());
  DrawAxisRAA(outFile, FwdPt());
  DrawAxisRAA(outFile, MidCent());
  DrawAxisRAA(outFile, FwdCent());
  outFile.Close();

  std::cout << "\nWrote ./roots/RAA_JPsi_effAccCompare_260519v2_noPtW.root" << std::endl;
  std::cout << "Wrote ./figs/RAA_JPsi_effAccCompare_260519v2_noPtW_{midRap_pT,fwdRap_pT,midRap_Npart,fwdRap_Npart_4Bins}_{prompt,nonprompt}.{pdf,png}" << std::endl;
}
