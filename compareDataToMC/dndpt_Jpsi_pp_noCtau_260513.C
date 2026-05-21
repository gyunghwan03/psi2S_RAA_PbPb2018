#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TNamed.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include "../commonUtility.h"
#include "../cutsAndBin.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace JpsiPpNoCtau260513
{
struct Config {
  bool prompt = true;
  bool forward = false;
  std::vector<double> ptBins;
  double yLow = 0.0;
  double yHigh = 1.6;
  TString rapidityLabel;
  TString componentLabel;
  TString mcPath;
  TString outPath;
  TString plotStem;
};

Config MakeConfig(bool prompt, bool forward)
{
  const TString repo = "/data/hwan/psi2S_RAA_PbPb2018";
  Config cfg;
  cfg.prompt = prompt;
  cfg.forward = forward;
  cfg.componentLabel = prompt ? "Jpsi" : "BtoJpsi";
  cfg.mcPath = prompt
                   ? Form("%s/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root", repo.Data())
                   : Form("%s/skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_miniAOD_251103.root", repo.Data());
  if (forward)
  {
    cfg.ptBins = {3.0, 4.0, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 40.0};
    cfg.yLow = 1.6;
    cfg.yHigh = 2.4;
    cfg.rapidityLabel = "y1p6_2p4";
  }
  else
  {
    cfg.ptBins = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5,
                  20.0, 22.5, 25.0, 27.5, 30.0, 40.0};
    cfg.yLow = 0.0;
    cfg.yHigh = 1.6;
    cfg.rapidityLabel = "y0_1p6";
  }

  cfg.outPath = Form("%s/compareDataToMC/ratioDataMC_pp_%s_DATA_noCtau_%s_260513.root",
                     repo.Data(), cfg.componentLabel.Data(), cfg.rapidityLabel.Data());
  cfg.plotStem = Form("%s/compareDataToMC/dNdpt_plot_%s_Jpsi_pp_noCtau_%s_260513",
                      repo.Data(), prompt ? "Prompt" : "NonPrompt", cfg.rapidityLabel.Data());
  return cfg;
}

valErr GetYield(const Config &cfg, double ptLow, double ptHigh)
{
  const TString repo = "/data/hwan/psi2S_RAA_PbPb2018";
  const TString kineLabel = getKineLabelpp(ptLow, ptHigh, cfg.yLow, cfg.yHigh, 0.0);
  const TString prDir = cfg.prompt ? "PRMC" : "NPMC";
  const TString path = Form("%s/Macros/Jpsi_L_cut/roots_1S_pp/%s/"
                            "Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root",
                            repo.Data(), prDir.Data(), kineLabel.Data());

  valErr ret;
  ret.val = 0.0;
  ret.err = 0.0;
  std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "[dndpt_Jpsi_pp_noCtau_260513] cannot open yield file: " << path << "\n";
    return ret;
  }

  TH1D *fitResults = dynamic_cast<TH1D *>(file->Get("fitResults"));
  if (!fitResults)
  {
    std::cerr << "[dndpt_Jpsi_pp_noCtau_260513] missing fitResults in " << path << "\n";
    return ret;
  }

  ret.val = fitResults->GetBinContent(1);
  ret.err = fitResults->GetBinError(1);
  return ret;
}

bool FillMcPt(const Config &cfg, TH1D *hPt, TH1D *hPtForRatio, TH1D *hFine)
{
  TChain tree("mmepevt");
  tree.Add(cfg.mcPath.Data());
  if (tree.GetEntries() <= 0)
  {
    std::cerr << "[dndpt_Jpsi_pp_noCtau_260513] no MC entries from " << cfg.mcPath << "\n";
    return false;
  }

  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu];
  Int_t event = 0;
  Int_t nDimu = 0;
  float vz = 0.0;
  int recoQQsign[nMaxDimu];

  tree.SetBranchAddress("event", &event);
  tree.SetBranchAddress("nDimu", &nDimu);
  tree.SetBranchAddress("vz", &vz);
  tree.SetBranchAddress("recoQQsign", recoQQsign);
  tree.SetBranchAddress("mass", mass);
  tree.SetBranchAddress("y", y);
  tree.SetBranchAddress("pt", pt);

  const double massLow = 2.6;
  const double massHigh = 3.5;
  const double ptMin = cfg.ptBins.front();
  const double ptMax = cfg.ptBins.back();
  const Long64_t nEvt = tree.GetEntries();
  std::cout << "[dndpt_Jpsi_pp_noCtau_260513] " << cfg.componentLabel << " "
            << cfg.rapidityLabel << " MC entries: " << nEvt << "\n";

  for (Long64_t i = 0; i < nEvt; ++i)
  {
    tree.GetEntry(i);
    for (int j = 0; j < nDimu; ++j)
    {
      const double absY = std::fabs(y[j]);
      if (!(mass[j] > massLow && mass[j] < massHigh &&
            pt[j] > ptMin && pt[j] < ptMax &&
            absY > cfg.yLow && absY < cfg.yHigh))
        continue;

      hPt->Fill(pt[j]);
      hPtForRatio->Fill(pt[j]);
      hFine->Fill(pt[j]);
    }
  }

  return hPt->Integral() > 0.0 && hPtForRatio->Integral() > 0.0;
}

TF1 *MakeRatioFit(const Config &cfg)
{
  const double xmin = cfg.forward ? 1.0 : 3.0;
  TF1 *fit = new TF1("fitRatio1",
                     "( [0] + [1]*x + [2]*x*x + [4]*x*x*x ) / "
                     "( (x-[3])*(x-[3])*(x-[3]) )",
                     xmin, 40.0);
  fit->SetNpx(1000);
  return fit;
}

TCanvas *BuildCanvas(const Config &cfg, TH1D *hData, TH1D *hMc, TH1D *ratio, TF1 *fit)
{
  TLegend *leg = new TLegend(0.65, 0.75, 0.85, 0.85);
  leg->AddEntry(hData, "Data", "p");
  leg->AddEntry(hMc, "MC", "l");
  leg->SetLineColor(kWhite);

  TCanvas *c_A = new TCanvas("canvas_A", "pp dN/dpT and Data/MC", 4, 4, 550, 520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
  pad_A_1->SetTicks(1, 1);
  pad_A_1->SetFillColor(0);
  pad_A_1->SetBorderMode(0);
  pad_A_1->SetBorderSize(2);
  pad_A_1->SetTopMargin(0.05646528);
  pad_A_1->Draw();
  pad_A_1->cd();

  hData->Draw();
  hMc->Draw("same hist");
  leg->Draw("same");
  hData->SetAxisRange(0.0, std::max(hData->GetMaximum(), hMc->GetMaximum()) + 0.05, "Y");
  hData->GetXaxis()->SetLabelSize(0);
  hData->GetYaxis()->SetTitleSize(0.04);
  hData->GetYaxis()->SetTitleOffset(1.00);
  hData->GetYaxis()->SetTitle("dN/dp_{T}");

  c_A->cd();
  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
  pad_A_2->SetFillColor(0);
  pad_A_2->SetBorderMode(0);
  pad_A_2->SetBorderSize(2);
  pad_A_2->SetTicks(1, 1);
  pad_A_2->SetBottomMargin(0.4361001);
  pad_A_2->Draw();
  pad_A_2->cd();

  ratio->SetName("WeightFactor");
  ratio->Draw();
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetXaxis()->SetTitleSize(0.15);
  ratio->GetXaxis()->CenterTitle();
  ratio->GetXaxis()->SetLabelOffset(0.04);
  ratio->GetXaxis()->SetLabelSize(0.15);
  ratio->GetXaxis()->SetTickSize(0.03);
  ratio->GetYaxis()->SetTickSize(0.04);
  ratio->GetYaxis()->SetNdivisions(404);
  ratio->GetYaxis()->SetTitle("Data/MC");
  ratio->GetYaxis()->SetTitleOffset(0.25);
  ratio->GetYaxis()->SetTitleSize(0.15);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetLabelSize(0.15);
  ratio->SetAxisRange(0.0, 4.0, "Y");

  if (fit)
  {
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    fit->Draw("same");
  }
  TLine *unity = new TLine(cfg.ptBins.front(), 1.0, cfg.ptBins.back(), 1.0);
  unity->SetLineStyle(2);
  unity->Draw("same");

  c_A->cd();
  c_A->SetTitle(Form("%s pp no-ctau %s", cfg.componentLabel.Data(), cfg.rapidityLabel.Data()));
  c_A->Modified();
  c_A->Update();
  return c_A;
}

void RunOne(bool prompt, bool forward, bool write)
{
  Config cfg = MakeConfig(prompt, forward);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  TH1::SetDefaultSumw2();

  const int nPtBins = static_cast<int>(cfg.ptBins.size()) - 1;
  TH1D *hData = new TH1D("hptData", ";p_{T}(GeV/c);", nPtBins, cfg.ptBins.data());
  TH1D *hDataRatio = new TH1D("hptData1", ";p_{T}(GeV/c);", nPtBins, cfg.ptBins.data());
  TH1D *hMc = new TH1D("hptMC", ";p_{T}(GeV/c);", nPtBins, cfg.ptBins.data());
  TH1D *hMcRatio = new TH1D("hptMC1", ";p_{T}(GeV/c);", nPtBins, cfg.ptBins.data());
  TH1D *hMcFine = new TH1D("hptMC2", ";p_{T}(GeV/c);", 100, 0.0, 50.0);

  if (!FillMcPt(cfg, hMc, hMcRatio, hMcFine))
    return;

  for (int ib = 1; ib <= nPtBins; ++ib)
  {
    const double ptLow = cfg.ptBins[ib - 1];
    const double ptHigh = cfg.ptBins[ib];
    const valErr yield = GetYield(cfg, ptLow, ptHigh);
    std::cout << "[yield] " << cfg.componentLabel << " " << cfg.rapidityLabel << " "
              << ptLow << "-" << ptHigh << ": " << yield.val << " +/- " << yield.err << "\n";
    hData->SetBinContent(ib, yield.val);
    hData->SetBinError(ib, yield.err);
    hDataRatio->SetBinContent(ib, yield.val);
    hDataRatio->SetBinError(ib, yield.err);
  }

  if (hData->Integral() <= 0.0 || hDataRatio->Integral() <= 0.0)
  {
    std::cerr << "[dndpt_Jpsi_pp_noCtau_260513] zero data integral for "
              << cfg.componentLabel << " " << cfg.rapidityLabel << "\n";
    return;
  }

  hMc->Scale(1.0 / hMc->Integral());
  hMcRatio->Scale(1.0 / hMcRatio->Integral());
  hMcFine->Scale(1.0 / hMcFine->Integral());
  hData->Scale(1.0 / hData->Integral());
  hDataRatio->Scale(1.0 / hDataRatio->Integral());

  TH1ScaleByWidth(hMc);
  TH1ScaleByWidth(hData);
  handsomeTH1(hMc, 1);
  handsomeTH1(hData, 1);
  handsomeTH1(hDataRatio, 1);

  hDataRatio->Divide(hMcRatio);
  std::unique_ptr<TF1> fitRatio(MakeRatioFit(cfg));
  hDataRatio->Fit(fitRatio.get(), "IE", "", 3.0, 12.0);
  TFitResultPtr fitResult = hDataRatio->Fit(fitRatio.get(), "S", "", 3.0, 40.0);
  if (fitResult.Get())
    fitResult.Get()->Print("V");

  std::unique_ptr<TCanvas> c_A(BuildCanvas(cfg, hData, hMc, hDataRatio, fitRatio.get()));

  if (!write)
    return;

  std::unique_ptr<TFile> out(TFile::Open(cfg.outPath, "RECREATE"));
  if (!out || out->IsZombie())
  {
    std::cerr << "[dndpt_Jpsi_pp_noCtau_260513] cannot write " << cfg.outPath << "\n";
    return;
  }

  hDataRatio->Write("WeightFactor");
  fitRatio->SetName("dataMC_Ratio1");
  fitRatio->Write();
  fitRatio->SetName("fitRatio1");
  fitRatio->Write();
  hData->Write("hptData_norm_width");
  hMc->Write("hptMC_norm_width");
  hMcFine->Write();
  if (c_A)
    c_A->Write();
  TNamed note("note", "pp no-ctau dN/dpT source: yield from pp mass fits, no 1/ctauEff correction, same binning as 251103 ctauCut dndpt macros");
  note.Write();
  out->Close();

  if (c_A)
  {
    c_A->SaveAs(Form("%s.pdf", cfg.plotStem.Data()));
    c_A->SaveAs(Form("%s.png", cfg.plotStem.Data()));
  }
  std::cout << "[wrote] " << cfg.outPath << "\n";
}

void RunAll()
{
  RunOne(true, false, true);
  RunOne(true, true, true);
  RunOne(false, false, true);
  RunOne(false, true, true);
}
}

void dndpt_Jpsi_pp_noCtau_260513()
{
  JpsiPpNoCtau260513::RunAll();
}
