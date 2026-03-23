#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <vector>
#include <cstdio>

using namespace std;

void compare_date_vs_date_fit_channel(TString system = "AA",
                                      TString state = "psi2S",
                                      int PR = 0,
                                      TString dateOld = "251103",
                                      TString dateNew = "260310");
void compare_split_251103_vs_260310(TString system = "AA",
                                    TString state = "psi2S",
                                    int PR = 0,
                                    TString dateSplit = "251103",
                                    TString dateRef = "260310");

namespace {

TString BuildInputFileName(TString system, TString state, int PR, bool is2exp) {
  TString prefix;
  TString tag;

  const bool isPrompt = (PR == 0);

  if (system == "AA") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310"
                        : "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_AA_psi2S_DATA_ctauCut_y0_2p4_260316"
                        : "ratioDataMC_AA_Btopsi2S_DATA_ctauCut_y0_2p4_260316";
    }
  } else if (system == "pp") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310"
                        : "ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_pp_psi2S_DATA_ctauCut_y0_2p4_260310"
                        : "ratioDataMC_pp_Btopsi2S_DATA_ctauCut_y0_2p4_260310";
    }
  }

  if (prefix.IsNull()) return "";

  return is2exp ? Form("%s_2exp.root", prefix.Data()) : Form("%s.root", prefix.Data());
}

TString BuildInputFileNameWithDate(TString system, TString state, int PR, TString dateTag) {
  TString prefix;

  const bool isPrompt = (PR == 0);

  if (system == "AA") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4"
                        : "ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_AA_psi2S_DATA_ctauCut_y0_2p4"
                        : "ratioDataMC_AA_Btopsi2S_DATA_ctauCut_y0_2p4";
    }
  } else if (system == "pp") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4"
                        : "ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_pp_psi2S_DATA_ctauCut_y0_2p4"
                        : "ratioDataMC_pp_Btopsi2S_DATA_ctauCut_y0_2p4";
    }
  }

  if (prefix.IsNull() || dateTag.IsNull()) return "";

  return Form("%s_%s.root", prefix.Data(), dateTag.Data());
}

TString BuildInputFileNameWithRapDate(TString system,
                                      TString state,
                                      int PR,
                                      TString rapTag,
                                      TString dateTag) {
  TString prefix;
  const bool isPrompt = (PR == 0);

  if (system == "AA") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_AA_Jpsi_DATA_ctauCut"
                        : "ratioDataMC_AA_BtoJpsi_DATA_ctauCut";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_AA_psi2S_DATA_ctauCut"
                        : "ratioDataMC_AA_Btopsi2S_DATA_ctauCut";
    }
  } else if (system == "pp") {
    if (state == "Jpsi") {
      prefix = isPrompt ? "ratioDataMC_pp_Jpsi_DATA_ctauCut"
                        : "ratioDataMC_pp_BtoJpsi_DATA_ctauCut";
    } else if (state == "psi2S") {
      prefix = isPrompt ? "ratioDataMC_pp_psi2S_DATA_ctauCut"
                        : "ratioDataMC_pp_Btopsi2S_DATA_ctauCut";
    }
  }

  if (prefix.IsNull() || rapTag.IsNull() || dateTag.IsNull()) return "";
  return Form("%s_%s_%s.root", prefix.Data(), rapTag.Data(), dateTag.Data());
}

TF1* GetRatioFunction(TFile* f, const char* newName) {
  if (!f || f->IsZombie()) return nullptr;

  TF1* func = dynamic_cast<TF1*>(f->Get("dataMC_Ratio1"));
  if (!func) return nullptr;

  TF1* cloned = dynamic_cast<TF1*>(func->Clone(newName));
  return cloned;
}

void SampleMinMax(TF1* f1, TF1* f2, double xMin, double xMax, double& yMin, double& yMax) {
  yMin = 1e9;
  yMax = -1e9;

  const int n = 300;
  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double y1 = f1->Eval(x);
    const double y2 = f2->Eval(x);
    if (TMath::Finite(y1)) {
      if (y1 < yMin) yMin = y1;
      if (y1 > yMax) yMax = y1;
    }
    if (TMath::Finite(y2)) {
      if (y2 < yMin) yMin = y2;
      if (y2 > yMax) yMax = y2;
    }
  }

  if (yMin >= yMax) {
    yMin = 0.0;
    yMax = 2.0;
  } else {
    const double dy = yMax - yMin;
    yMin -= 0.15 * dy;
    yMax += 0.20 * dy;
  }
}

void SampleMinMax3(TF1* f1,
                   TF1* f2,
                   TF1* f3,
                   double xMin,
                   double xMax,
                   double& yMin,
                   double& yMax) {
  yMin = 1e9;
  yMax = -1e9;

  const int n = 300;
  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double y1 = f1->Eval(x);
    const double y2 = f2->Eval(x);
    const double y3 = f3->Eval(x);
    if (TMath::Finite(y1)) { if (y1 < yMin) yMin = y1; if (y1 > yMax) yMax = y1; }
    if (TMath::Finite(y2)) { if (y2 < yMin) yMin = y2; if (y2 > yMax) yMax = y2; }
    if (TMath::Finite(y3)) { if (y3 < yMin) yMin = y3; if (y3 > yMax) yMax = y3; }
  }

  if (yMin >= yMax) {
    yMin = 0.0;
    yMax = 2.0;
  } else {
    const double dy = yMax - yMin;
    yMin -= 0.15 * dy;
    yMax += 0.20 * dy;
  }
}

TGraph* BuildRatioGraph(TF1* fOld, TF1* fNew, double xMin, double xMax) {
  vector<double> vx;
  vector<double> vy;

  const int n = 300;
  vx.reserve(n);
  vy.reserve(n);

  for (int i = 0; i < n; ++i) {
    const double x = xMin + (xMax - xMin) * (i + 0.5) / n;
    const double yOld = fOld->Eval(x);
    const double yNew = fNew->Eval(x);
    if (!TMath::Finite(yOld) || !TMath::Finite(yNew)) continue;
    if (TMath::Abs(yOld) < 1e-12) continue;

    vx.push_back(x);
    vy.push_back(yNew / yOld);
  }

  if (vx.empty()) return nullptr;

  TGraph* g = new TGraph(static_cast<int>(vx.size()), &vx[0], &vy[0]);
  g->SetName("g_new_over_old");
  return g;
}

} // namespace

void compare_old_vs_2exp_fit_files(const char* oldFile,
                                   const char* newFile,
                                   const char* outTag,
                                   const char* plotLabel = "",
                                   const char* oldLegend = "Old fit",
                                   const char* newLegend = "New fit (2exp+const)") {
  gStyle->SetOptStat(0);

  TFile* fOld = TFile::Open(oldFile, "READ");
  TFile* fNew = TFile::Open(newFile, "READ");

  cout << "[INFO] old file: " << oldFile << endl;
  cout << "[INFO] new file: " << newFile << endl;

  if (!fOld || fOld->IsZombie()) {
    cout << "[ERROR] Cannot open old file: " << oldFile << endl;
    return;
  }
  if (!fNew || fNew->IsZombie()) {
    cout << "[ERROR] Cannot open new file: " << newFile << endl;
    fOld->Close();
    return;
  }

  TF1* oldFunc = GetRatioFunction(fOld, "oldFunc");
  TF1* newFunc = GetRatioFunction(fNew, "newFunc");

  if (!oldFunc || !newFunc) {
    cout << "[ERROR] dataMC_Ratio1 not found in one of files" << endl;
    if (fOld) fOld->Close();
    if (fNew) fNew->Close();
    return;
  }

  const double xMin = TMath::Max(oldFunc->GetXmin(), newFunc->GetXmin());
  const double xMax = TMath::Min(oldFunc->GetXmax(), newFunc->GetXmax());

  cout << "[INFO] old formula: " << oldFunc->GetExpFormula("p").Data() << endl;
  cout << "[INFO] new formula: " << newFunc->GetExpFormula("p").Data() << endl;
  cout << "[INFO] old fit range: [" << oldFunc->GetXmin() << ", " << oldFunc->GetXmax() << "]" << endl;
  cout << "[INFO] new fit range: [" << newFunc->GetXmin() << ", " << newFunc->GetXmax() << "]" << endl;
  cout << "[INFO] common range : [" << xMin << ", " << xMax << "]" << endl;
  for (int i = 0; i < oldFunc->GetNpar(); ++i) {
    cout << "[INFO] old p" << i << " = " << oldFunc->GetParameter(i) << endl;
  }
  for (int i = 0; i < newFunc->GetNpar(); ++i) {
    cout << "[INFO] new p" << i << " = " << newFunc->GetParameter(i) << endl;
  }
  if (xMin >= xMax) {
    cout << "[ERROR] Invalid common fit range" << endl;
    fOld->Close();
    fNew->Close();
    return;
  }

  oldFunc->SetLineColor(kBlue + 1);
  oldFunc->SetLineWidth(3);
  oldFunc->SetLineStyle(1);

  newFunc->SetLineColor(kRed + 1);
  newFunc->SetLineWidth(3);
  newFunc->SetLineStyle(2);

  double yMin = 0.0, yMax = 2.0;
  SampleMinMax(oldFunc, newFunc, xMin, xMax, yMin, yMax);

  TCanvas* c = new TCanvas("c_compare_old_new", "old vs 2exp", 850, 800);
  TPad* pTop = new TPad("pTop", "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad("pBot", "", 0.0, 0.0, 1.0, 0.30);

  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.04);
  pTop->SetTicks(1, 1);

  pBot->SetTopMargin(0.03);
  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.04);
  pBot->SetTicks(1, 1);

  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  TH1D* hFrame = new TH1D("hFrame", ";p_{T} (GeV/c);Data/MC fit", 100, xMin, xMax);
  hFrame->SetMinimum(yMin);
  hFrame->SetMaximum(yMax);
  hFrame->GetXaxis()->SetLabelSize(0.0);
  hFrame->Draw("axis");

  oldFunc->Draw("same");
  newFunc->Draw("same");

  TLegend* leg = new TLegend(0.58, 0.73, 0.90, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(oldFunc, oldLegend, "l");
  leg->AddEntry(newFunc, newLegend, "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.035);
  if (plotLabel && TString(plotLabel).Length() > 0) {
    tx.DrawLatex(0.14, 0.92, plotLabel);
  }

  pBot->cd();
  TH1D* hRatioFrame = new TH1D("hRatioFrame", ";p_{T} (GeV/c);new/old", 100, xMin, xMax);
  hRatioFrame->SetMinimum(0.6);
  hRatioFrame->SetMaximum(1.4);
  hRatioFrame->GetXaxis()->SetTitleSize(0.12);
  hRatioFrame->GetXaxis()->SetLabelSize(0.10);
  hRatioFrame->GetYaxis()->SetTitleSize(0.10);
  hRatioFrame->GetYaxis()->SetLabelSize(0.09);
  hRatioFrame->GetYaxis()->SetTitleOffset(0.50);
  hRatioFrame->GetYaxis()->SetNdivisions(505);
  hRatioFrame->Draw("axis");

  TGraph* gRatio = BuildRatioGraph(oldFunc, newFunc, xMin, xMax);
  if (gRatio) {
    gRatio->SetLineColor(kBlack);
    gRatio->SetLineWidth(2);
    gRatio->Draw("L same");
  }

  TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
  unity->SetLineColor(kGray + 2);
  unity->SetLineStyle(2);
  unity->SetLineWidth(2);
  unity->Draw("same");

  gSystem->mkdir("./figs", kTRUE);
  c->SaveAs(Form("./figs/%s.pdf", outTag));
  c->SaveAs(Form("./figs/%s.png", outTag));

  cout << "[INFO] Saved: ./figs/" << outTag << ".pdf/.png" << endl;

  fOld->Close();
  fNew->Close();
}

void compare_old_vs_2exp_fit_channel(TString system = "AA",
                                     TString state = "Jpsi",
                                     int PR = 0) {
  TString oldFile = BuildInputFileName(system, state, PR, false);
  TString newFile = BuildInputFileName(system, state, PR, true);

  if (oldFile.IsNull() || newFile.IsNull()) {
    cout << "[ERROR] Invalid channel input. system=AA/pp, state=Jpsi/psi2S, PR=0/1" << endl;
    return;
  }

  const char* promptLabel = (PR == 0) ? "Prompt" : "NonPrompt";
  TString outTag = Form("compare_old_vs_2exp_%s_%s_%s_y0_2p4",
                        system.Data(), state.Data(), promptLabel);
  TString label = Form("%s %s %s, |y|<2.4", system.Data(), state.Data(), promptLabel);

  compare_old_vs_2exp_fit_files(oldFile.Data(), newFile.Data(), outTag.Data(), label.Data());
}

void compare_old_vs_2exp_fit(TString system = "AA",
                             TString state = "Jpsi",
                             int PR = 0) {
  // Support rq string-dispatch style with escaped or plain quotes.
  if (system.Contains("compare_split_251103_vs_260310") ||
      system.Contains("compare_date_vs_date_fit_channel")) {
    TString callExpr = system;
    callExpr.ReplaceAll("\\\"", "\"");

    if (callExpr.Contains("compare_split_251103_vs_260310")) {
      char sysBuf[32] = {0}, stateBuf[32] = {0}, dateSplitBuf[32] = {0}, dateRefBuf[32] = {0};
      int prBuf = 0;
      const int n = std::sscanf(
        callExpr.Data(),
        "compare_split_251103_vs_260310(\"%31[^\"]\",\"%31[^\"]\",%d,\"%31[^\"]\",\"%31[^\"]\")",
        sysBuf, stateBuf, &prBuf, dateSplitBuf, dateRefBuf);
      if (n == 5) {
        compare_split_251103_vs_260310(sysBuf, stateBuf, prBuf, dateSplitBuf, dateRefBuf);
      } else {
        compare_split_251103_vs_260310();
      }
      return;
    }

    if (callExpr.Contains("compare_date_vs_date_fit_channel")) {
      char sysBuf[32] = {0}, stateBuf[32] = {0}, dateOldBuf[32] = {0}, dateNewBuf[32] = {0};
      int prBuf = 0;
      const int n = std::sscanf(
        callExpr.Data(),
        "compare_date_vs_date_fit_channel(\"%31[^\"]\",\"%31[^\"]\",%d,\"%31[^\"]\",\"%31[^\"]\")",
        sysBuf, stateBuf, &prBuf, dateOldBuf, dateNewBuf);
      if (n == 5) {
        compare_date_vs_date_fit_channel(sysBuf, stateBuf, prBuf, dateOldBuf, dateNewBuf);
      } else {
        compare_date_vs_date_fit_channel();
      }
      return;
    }
  }

  compare_old_vs_2exp_fit_channel(system, state, PR);
}

void compare_old_vs_2exp_fit(TString system,
                             TString state,
                             int PR,
                             TString dateOld,
                             TString dateNew) {
  compare_date_vs_date_fit_channel(system, state, PR, dateOld, dateNew);
}

void compare_old_vs_2exp_fit_all_y0_2p4() {
  compare_old_vs_2exp_fit_channel("AA", "Jpsi", 0);
  compare_old_vs_2exp_fit_channel("AA", "Jpsi", 1);
  compare_old_vs_2exp_fit_channel("AA", "psi2S", 0);
  compare_old_vs_2exp_fit_channel("AA", "psi2S", 1);
  compare_old_vs_2exp_fit_channel("pp", "Jpsi", 0);
  compare_old_vs_2exp_fit_channel("pp", "Jpsi", 1);
  compare_old_vs_2exp_fit_channel("pp", "psi2S", 0);
  compare_old_vs_2exp_fit_channel("pp", "psi2S", 1);
}

void compare_date_vs_date_fit_channel(TString system,
                                      TString state,
                                      int PR,
                                      TString dateOld,
                                      TString dateNew) {
  TString oldFile = BuildInputFileNameWithDate(system, state, PR, dateOld);
  TString newFile = BuildInputFileNameWithDate(system, state, PR, dateNew);

  if (oldFile.IsNull() || newFile.IsNull()) {
    cout << "[ERROR] Invalid channel input. system=AA/pp, state=Jpsi/psi2S, PR=0/1" << endl;
    return;
  }

  const char* promptLabel = (PR == 0) ? "Prompt" : "NonPrompt";
  TString outTag = Form("compare_%s_vs_%s_%s_%s_%s_y0_2p4",
                        dateNew.Data(), dateOld.Data(),
                        system.Data(), state.Data(), promptLabel);
  TString label = Form("%s %s %s, |y|<2.4", system.Data(), state.Data(), promptLabel);
  TString legOld = Form("%s", dateOld.Data());
  TString legNew = Form("%s", dateNew.Data());

  compare_old_vs_2exp_fit_files(oldFile.Data(),
                                newFile.Data(),
                                outTag.Data(),
                                label.Data(),
                                legOld.Data(),
                                legNew.Data());
}

void compare_split_251103_vs_260310(TString system,
                                    TString state,
                                    int PR,
                                    TString dateSplit,
                                    TString dateRef) {
  gStyle->SetOptStat(0);

  TString fileSplitLow  = BuildInputFileNameWithRapDate(system, state, PR, "y0_1p6", dateSplit);
  TString fileSplitHigh = BuildInputFileNameWithRapDate(system, state, PR, "y1p6_2p4", dateSplit);
  TString fileRef       = BuildInputFileNameWithRapDate(system, state, PR, "y0_2p4", dateRef);

  if (fileSplitLow.IsNull() || fileSplitHigh.IsNull() || fileRef.IsNull()) {
    cout << "[ERROR] Invalid channel input. system=AA/pp, state=Jpsi/psi2S, PR=0/1" << endl;
    return;
  }

  TFile* fLow = TFile::Open(fileSplitLow, "READ");
  TFile* fHigh = TFile::Open(fileSplitHigh, "READ");
  TFile* fRef = TFile::Open(fileRef, "READ");

  if (!fLow || fLow->IsZombie()) { cout << "[ERROR] Cannot open " << fileSplitLow << endl; return; }
  if (!fHigh || fHigh->IsZombie()) { cout << "[ERROR] Cannot open " << fileSplitHigh << endl; fLow->Close(); return; }
  if (!fRef || fRef->IsZombie()) { cout << "[ERROR] Cannot open " << fileRef << endl; fLow->Close(); fHigh->Close(); return; }

  TF1* funcLow = GetRatioFunction(fLow, "funcSplitLow");
  TF1* funcHigh = GetRatioFunction(fHigh, "funcSplitHigh");
  TF1* funcRef = GetRatioFunction(fRef, "funcRef");
  if (!funcLow || !funcHigh || !funcRef) {
    cout << "[ERROR] dataMC_Ratio1 not found in one of files" << endl;
    fLow->Close(); fHigh->Close(); fRef->Close();
    return;
  }

  const double xMin = TMath::Max(funcRef->GetXmin(), TMath::Max(funcLow->GetXmin(), funcHigh->GetXmin()));
  const double xMax = TMath::Min(funcRef->GetXmax(), TMath::Min(funcLow->GetXmax(), funcHigh->GetXmax()));
  if (xMin >= xMax) {
    cout << "[ERROR] Invalid common fit range" << endl;
    fLow->Close(); fHigh->Close(); fRef->Close();
    return;
  }

  funcLow->SetLineColor(kBlue + 1);  funcLow->SetLineWidth(3);
  funcHigh->SetLineColor(kGreen + 2); funcHigh->SetLineWidth(3);
  funcRef->SetLineColor(kRed + 1);    funcRef->SetLineWidth(3); funcRef->SetLineStyle(2);

  double yMin = 0.0, yMax = 2.0;
  SampleMinMax3(funcLow, funcHigh, funcRef, xMin, xMax, yMin, yMax);

  TCanvas* c = new TCanvas("c_compare_split_vs_ref", "split vs ref", 850, 800);
  TPad* pTop = new TPad("pTop_split", "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad("pBot_split", "", 0.0, 0.0, 1.0, 0.30);
  pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.12); pTop->SetRightMargin(0.04); pTop->SetTicks(1,1);
  pBot->SetTopMargin(0.03); pBot->SetBottomMargin(0.35); pBot->SetLeftMargin(0.12); pBot->SetRightMargin(0.04); pBot->SetTicks(1,1);
  pTop->Draw(); pBot->Draw();

  pTop->cd();
  TH1D* hFrame = new TH1D("hFrame_split", ";p_{T} (GeV/c);Data/MC fit", 100, xMin, xMax);
  hFrame->SetMinimum(yMin); hFrame->SetMaximum(yMax); hFrame->GetXaxis()->SetLabelSize(0.0); hFrame->Draw("axis");
  funcLow->Draw("same"); funcHigh->Draw("same"); funcRef->Draw("same");

  TLegend* leg = new TLegend(0.47, 0.66, 0.90, 0.89);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(funcLow, Form("%s, |y|<1.6", dateSplit.Data()), "l");
  leg->AddEntry(funcHigh, Form("%s, 1.6<|y|<2.4", dateSplit.Data()), "l");
  leg->AddEntry(funcRef, Form("%s, |y|<2.4", dateRef.Data()), "l");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.DrawLatex(0.14, 0.92, Form("%s %s %s", system.Data(), state.Data(), (PR == 0 ? "Prompt" : "NonPrompt")));

  pBot->cd();
  TH1D* hRatioFrame = new TH1D("hRatioFrame_split", ";p_{T} (GeV/c);split/ref", 100, xMin, xMax);
  hRatioFrame->SetMinimum(0.6); hRatioFrame->SetMaximum(1.4);
  hRatioFrame->GetXaxis()->SetTitleSize(0.12); hRatioFrame->GetXaxis()->SetLabelSize(0.10);
  hRatioFrame->GetYaxis()->SetTitleSize(0.10); hRatioFrame->GetYaxis()->SetLabelSize(0.09);
  hRatioFrame->GetYaxis()->SetTitleOffset(0.50); hRatioFrame->GetYaxis()->SetNdivisions(505);
  hRatioFrame->Draw("axis");

  TGraph* gLowOverRef = BuildRatioGraph(funcRef, funcLow, xMin, xMax);
  if (gLowOverRef) { gLowOverRef->SetLineColor(kBlue + 1); gLowOverRef->SetLineWidth(2); gLowOverRef->Draw("L same"); }
  TGraph* gHighOverRef = BuildRatioGraph(funcRef, funcHigh, xMin, xMax);
  if (gHighOverRef) { gHighOverRef->SetLineColor(kGreen + 2); gHighOverRef->SetLineWidth(2); gHighOverRef->Draw("L same"); }

  TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
  unity->SetLineColor(kGray + 2); unity->SetLineStyle(2); unity->SetLineWidth(2); unity->Draw("same");

  TString outTag = Form("compare_split%s_vs_%s_%s_%s_%s_y0_2p4",
                        dateSplit.Data(), dateRef.Data(), system.Data(), state.Data(),
                        (PR == 0 ? "Prompt" : "NonPrompt"));
  gSystem->mkdir("./figs", kTRUE);
  c->SaveAs(Form("./figs/%s.pdf", outTag.Data()));
  c->SaveAs(Form("./figs/%s.png", outTag.Data()));
  cout << "[INFO] Saved: ./figs/" << outTag << ".pdf/.png" << endl;

  fLow->Close(); fHigh->Close(); fRef->Close();
}
