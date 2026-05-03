#include "/data/hwan/psi2S_RAA_PbPb2018/Macros/final_Results/Gwak_raa_comparison/jpsi_raa_values.h"

#include "TFile.h"
#include "TH1.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kRepo = "/data/hwan/psi2S_RAA_PbPb2018";
const std::string kGwakDir = "/data/users/pjgwak/work/daily_code_tracker/2026/04/11_run2_eff_study/acc_eff_skim/skim_roots";
const std::string kOutCsv = kRepo + "/Macros/final_Results/Gwak_raa_comparison/raa_acc_eff_compare_260501.csv";
const std::string kOutMd = kRepo + "/Macros/final_Results/Gwak_raa_comparison/raa_acc_eff_compare_260501.md";

double NaN()
{
  return std::numeric_limits<double>::quiet_NaN();
}

std::string RefAccFile(bool isPbPb, bool isPrompt)
{
  return kGwakDir + "/acc_" + (isPbPb ? "PbPb2018" : "pp2018") + "_ppInput_isMC1_" +
         (isPrompt ? "PR" : "NP") + "_ncollW0_genW1_ptW1.root";
}

std::string NewAccFile(bool isPbPb, bool isPrompt)
{
  return kRepo + "/Eff_Acc_260429/roots/acceptance_" + (isPrompt ? "PromptJpsi" : "BtoJpsi") +
         "_GenOnly_wgt1_" + (isPbPb ? "PbPb" : "pp") + "_SysUp0_260429_1d.root";
}

std::string RefEffFile(bool isPbPb, bool isPrompt)
{
  if (isPbPb)
    return kGwakDir + "/eff_PbPb2018_isMC1_" + (isPrompt ? std::string("PR") : std::string("NP")) +
           "_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root";
  return kGwakDir + "/eff_pp5p02TeV_isMC1_" + (isPrompt ? std::string("PR") : std::string("NP")) +
         "_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root";
}

std::string NewEffFile(bool isPbPb, bool isPrompt)
{
  if (isPbPb)
    return kRepo + "/Eff_Acc_260429/roots/mc_eff_vs_pt_cent_0_to_180_rap_" +
           (isPrompt ? std::string("prompt") : std::string("nprompt")) +
           "_pbpb_JPsi_PtWnomi_tnp1_260430_1d.root";
  return kRepo + "/Eff_Acc_260429/roots/mc_eff_vs_pt_rap_" +
         (isPrompt ? std::string("prompt") : std::string("nprompt")) +
         "_pp_Jpsi_PtWnomi_tnp1_260429_1d.root";
}

TFile *OpenFile(const std::string &path)
{
  TFile *file = TFile::Open(path.c_str(), "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "[ERROR] cannot open " << path << "\n";
    return nullptr;
  }
  return file;
}

TH1 *GetHist(TFile *file, const char *key)
{
  if (!file)
    return nullptr;
  TH1 *hist = dynamic_cast<TH1 *>(file->Get(key));
  if (!hist)
    std::cerr << "[ERROR] missing " << key << " in " << file->GetName() << "\n";
  return hist;
}

double BinValue(TH1 *hist, int bin)
{
  if (!hist || bin < 1 || bin > hist->GetNbinsX())
    return NaN();
  return hist->GetBinContent(bin);
}

double SumRatio(TFile *file, const char *numKey, const char *denKey)
{
  TH1 *num = GetHist(file, numKey);
  TH1 *den = GetHist(file, denKey);
  if (!num || !den)
    return NaN();

  double sumNum = 0.0;
  double sumDen = 0.0;
  for (int ib = 1; ib <= num->GetNbinsX(); ++ib)
    sumNum += num->GetBinContent(ib);
  for (int ib = 1; ib <= den->GetNbinsX(); ++ib)
    sumDen += den->GetBinContent(ib);

  if (sumDen == 0.0)
    return NaN();
  return sumNum / sumDen;
}

double RatioOrBin(double num, double den)
{
  if (!std::isfinite(num) || !std::isfinite(den) || den == 0.0)
    return NaN();
  return num / den;
}

std::string F(double value, int precision = 9)
{
  if (!std::isfinite(value))
    return "nan";
  std::ostringstream os;
  os << std::fixed << std::setprecision(precision) << value;
  return os.str();
}

std::string BinLabel(const double *edges, int idx)
{
  std::ostringstream os;
  os << F(edges[idx], 1) << "-" << F(edges[idx + 1], 1);
  return os.str();
}

struct Row
{
  std::string observable;
  std::string state;
  std::string bin;
  double bakRaa;
  double gwakRaa;
  double raaRatio;
  double raaDiffPct;
  double accPpGwak;
  double accPbPbGwak;
  double accPpBak;
  double accPbPbBak;
  double accOnlyRatio;
  double effPpGwak;
  double effPbPbGwak;
  double effPpBak;
  double effPbPbBak;
  double accEffRatio;
  double residualAfterAccPct;
  double residualAfterAccEffPct;
};

void PushPtRows(std::vector<Row> &rows,
                const char *observable,
                const char *raaFileName,
                const double *edges,
                const double *refPr,
                const double *refNp,
                int nBins,
                bool isMid)
{
  TFile *raaFile = OpenFile(kRepo + "/Macros/final_Results/260430/roots/" + raaFileName);
  TH1 *raaPr = GetHist(raaFile, "hRAA_PR");
  TH1 *raaNp = GetHist(raaFile, "hRAA_NP");

  const char *accKey = isMid ? "hAccPt_2021_midy" : "hAccPt_2021_Fory";
  const char *refAccKey = isMid ? "hist_acc_mid" : "hist_acc_fwd";
  const char *effKeyPp = isMid ? "mc_eff_vs_pt_TnP1_PtW1_absy0_1p6" : "mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4";
  const char *effKeyPbPb = isMid ? "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6"
                                 : "mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4";
  const char *refEffKey = isMid ? "hist_eff_mid" : "hist_eff_fwd";

  for (int stateIdx = 0; stateIdx < 2; ++stateIdx)
  {
    const bool isPrompt = stateIdx == 0;
    const char *state = isPrompt ? "PR" : "NP";
    TH1 *raa = isPrompt ? raaPr : raaNp;
    const double *ref = isPrompt ? refPr : refNp;

    TFile *refAccPpFile = OpenFile(RefAccFile(false, isPrompt));
    TFile *refAccPbPbFile = OpenFile(RefAccFile(true, isPrompt));
    TFile *newAccPpFile = OpenFile(NewAccFile(false, isPrompt));
    TFile *newAccPbPbFile = OpenFile(NewAccFile(true, isPrompt));
    TFile *refEffPpFile = OpenFile(RefEffFile(false, isPrompt));
    TFile *refEffPbPbFile = OpenFile(RefEffFile(true, isPrompt));
    TFile *newEffPpFile = OpenFile(NewEffFile(false, isPrompt));
    TFile *newEffPbPbFile = OpenFile(NewEffFile(true, isPrompt));

    TH1 *refAccPp = GetHist(refAccPpFile, refAccKey);
    TH1 *refAccPbPb = GetHist(refAccPbPbFile, refAccKey);
    TH1 *newAccPp = GetHist(newAccPpFile, accKey);
    TH1 *newAccPbPb = GetHist(newAccPbPbFile, accKey);
    TH1 *refEffPp = GetHist(refEffPpFile, refEffKey);
    TH1 *refEffPbPb = GetHist(refEffPbPbFile, refEffKey);
    TH1 *newEffPp = GetHist(newEffPpFile, effKeyPp);
    TH1 *newEffPbPb = GetHist(newEffPbPbFile, effKeyPbPb);

    for (int ib = 1; ib <= nBins; ++ib)
    {
      const double gwakRaa = ref[ib - 1];
      const double bakRaa = BinValue(raa, ib);
      const double raaRatio = RatioOrBin(bakRaa, gwakRaa);

      const double accPpGwak = BinValue(refAccPp, ib);
      const double accPbPbGwak = BinValue(refAccPbPb, ib);
      const double accPpBak = BinValue(newAccPp, ib);
      const double accPbPbBak = BinValue(newAccPbPb, ib);
      const double accOnlyRatio = RatioOrBin(RatioOrBin(accPpBak, accPbPbBak), RatioOrBin(accPpGwak, accPbPbGwak));

      const double effPpGwak = BinValue(refEffPp, ib);
      const double effPbPbGwak = BinValue(refEffPbPb, ib);
      const double effPpBak = BinValue(newEffPp, ib);
      const double effPbPbBak = BinValue(newEffPbPb, ib);
      const double accEffRatio = RatioOrBin(RatioOrBin(accPpBak * effPpBak, accPbPbBak * effPbPbBak),
                                           RatioOrBin(accPpGwak * effPpGwak, accPbPbGwak * effPbPbGwak));

      rows.push_back({observable,
                      state,
                      BinLabel(edges, ib - 1),
                      bakRaa,
                      gwakRaa,
                      raaRatio,
                      100.0 * (raaRatio - 1.0),
                      accPpGwak,
                      accPbPbGwak,
                      accPpBak,
                      accPbPbBak,
                      accOnlyRatio,
                      effPpGwak,
                      effPbPbGwak,
                      effPpBak,
                      effPbPbBak,
                      accEffRatio,
                      100.0 * (RatioOrBin(raaRatio, accOnlyRatio) - 1.0),
                      100.0 * (RatioOrBin(raaRatio, accEffRatio) - 1.0)});
    }

    delete refAccPpFile;
    delete refAccPbPbFile;
    delete newAccPpFile;
    delete newAccPbPbFile;
    delete refEffPpFile;
    delete refEffPbPbFile;
    delete newEffPpFile;
    delete newEffPbPbFile;
  }

  delete raaFile;
}

void PushCentRows(std::vector<Row> &rows,
                  const char *observable,
                  const char *raaFileName,
                  const double *npart,
                  const double *refPr,
                  const double *refNp,
                  int nBins,
                  bool isMid)
{
  TFile *raaFile = OpenFile(kRepo + "/Macros/final_Results/260430/roots/" + raaFileName);
  TH1 *raaPr = GetHist(raaFile, "hRAA_PR");
  TH1 *raaNp = GetHist(raaFile, "hRAA_NP");

  const char *accIntKey = isMid ? "hAccPt_2021_midy_Int" : "hAccPt_2021_Fory_Int";
  const char *refAccNumKey = isMid ? "hist_acc_num_mid" : "hist_acc_num_fwd";
  const char *refAccDenKey = isMid ? "hist_acc_den_mid" : "hist_acc_den_fwd";
  const char *refEffNumKey = isMid ? "hist_eff_num_mid" : "hist_eff_num_fwd";
  const char *refEffDenKey = isMid ? "hist_eff_den_mid" : "hist_eff_den_fwd";
  const char *effCentKey = isMid ? "mc_eff_vs_cent_TnP1_PtW1_pt_6p5_to_40_absy0_1p6"
                                 : "mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4";
  const char *refEffCentKey = isMid ? "hist_eff_cent_mid" : "hist_eff_cent_fwd";

  for (int stateIdx = 0; stateIdx < 2; ++stateIdx)
  {
    const bool isPrompt = stateIdx == 0;
    const char *state = isPrompt ? "PR" : "NP";
    TH1 *raa = isPrompt ? raaPr : raaNp;
    const double *ref = isPrompt ? refPr : refNp;

    TFile *refAccPpFile = OpenFile(RefAccFile(false, isPrompt));
    TFile *refAccPbPbFile = OpenFile(RefAccFile(true, isPrompt));
    TFile *newAccPpFile = OpenFile(NewAccFile(false, isPrompt));
    TFile *newAccPbPbFile = OpenFile(NewAccFile(true, isPrompt));
    TFile *refEffPpFile = OpenFile(RefEffFile(false, isPrompt));
    TFile *refEffPbPbFile = OpenFile(RefEffFile(true, isPrompt));
    TFile *newEffPpFile = OpenFile(NewEffFile(false, isPrompt));
    TFile *newEffPbPbFile = OpenFile(NewEffFile(true, isPrompt));

    TH1 *newAccPp = GetHist(newAccPpFile, accIntKey);
    TH1 *newAccPbPb = GetHist(newAccPbPbFile, accIntKey);
    TH1 *refEffPbPb = GetHist(refEffPbPbFile, refEffCentKey);
    TH1 *newEffPbPb = GetHist(newEffPbPbFile, effCentKey);

    const double accPpGwak = SumRatio(refAccPpFile, refAccNumKey, refAccDenKey);
    const double accPbPbGwak = SumRatio(refAccPbPbFile, refAccNumKey, refAccDenKey);
    const double accPpBak = BinValue(newAccPp, 1);
    const double accPbPbBak = BinValue(newAccPbPb, 1);
    const double accOnlyRatio = RatioOrBin(RatioOrBin(accPpBak, accPbPbBak), RatioOrBin(accPpGwak, accPbPbGwak));
    const double effPpGwak = SumRatio(refEffPpFile, refEffNumKey, refEffDenKey);
    const double effPpBak = SumRatio(newEffPpFile, refEffNumKey, refEffDenKey);

    for (int ib = 1; ib <= nBins; ++ib)
    {
      // RAA output bins and Npart are stored in peripheral-to-central order.
      // The header centrality arrays are drawn through makeCentGraph(), which
      // reverses y[nBins - 1 - i], so use the same convention here.
      const int effBin = nBins - ib + 1;
      const double gwakRaa = ref[nBins - ib];
      const double bakRaa = BinValue(raa, ib);
      const double raaRatio = RatioOrBin(bakRaa, gwakRaa);
      const double effPbPbGwak = BinValue(refEffPbPb, effBin);
      const double effPbPbBak = BinValue(newEffPbPb, effBin);
      const double accEffRatio = RatioOrBin(RatioOrBin(accPpBak * effPpBak, accPbPbBak * effPbPbBak),
                                           RatioOrBin(accPpGwak * effPpGwak, accPbPbGwak * effPbPbGwak));

      rows.push_back({observable,
                      state,
                      "Npart=" + F(npart[ib - 1], 2),
                      bakRaa,
                      gwakRaa,
                      raaRatio,
                      100.0 * (raaRatio - 1.0),
                      accPpGwak,
                      accPbPbGwak,
                      accPpBak,
                      accPbPbBak,
                      accOnlyRatio,
                      effPpGwak,
                      effPbPbGwak,
                      effPpBak,
                      effPbPbBak,
                      accEffRatio,
                      100.0 * (RatioOrBin(raaRatio, accOnlyRatio) - 1.0),
                      100.0 * (RatioOrBin(raaRatio, accEffRatio) - 1.0)});
    }

    delete refAccPpFile;
    delete refAccPbPbFile;
    delete newAccPpFile;
    delete newAccPbPbFile;
    delete refEffPpFile;
    delete refEffPbPbFile;
    delete newEffPpFile;
    delete newEffPbPbFile;
  }

  delete raaFile;
}

void WriteCsv(const std::vector<Row> &rows)
{
  std::ofstream out(kOutCsv);
  out << "observable,state,bin,bak_raa,gwak_header_raa,bak_over_gwak_raa,raa_diff_pct,"
         "gwak_acc_pp,gwak_acc_pbpb,bak_acc_pp,bak_acc_pbpb,acc_only_ratio,"
         "gwak_eff_pp,gwak_eff_pbpb,bak_eff_pp,bak_eff_pbpb,acc_eff_ratio,"
         "residual_after_acc_pct,residual_after_acc_eff_pct\n";
  for (const Row &r : rows)
  {
    out << r.observable << "," << r.state << "," << r.bin << ","
        << F(r.bakRaa) << "," << F(r.gwakRaa) << "," << F(r.raaRatio) << "," << F(r.raaDiffPct) << ","
        << F(r.accPpGwak) << "," << F(r.accPbPbGwak) << "," << F(r.accPpBak) << "," << F(r.accPbPbBak) << ","
        << F(r.accOnlyRatio) << "," << F(r.effPpGwak) << "," << F(r.effPbPbGwak) << ","
        << F(r.effPpBak) << "," << F(r.effPbPbBak) << "," << F(r.accEffRatio) << ","
        << F(r.residualAfterAccPct) << "," << F(r.residualAfterAccEffPct) << "\n";
  }
}

void WriteMarkdown(const std::vector<Row> &rows)
{
  std::ofstream out(kOutMd);
  out << "# J/psi RAA vs Gwak Header and Acc/Eff Inputs\n\n";
  out << "Generated by `compare_raa_acc_eff_260501.C`.\n\n";
  out << "Formula checked: `RAA ~ raw_ratio * (A_pp * eff_pp) / (A_PbPb * eff_PbPb)`.\n";
  out << "`acc_only_ratio` is `(A_pp/A_PbPb)_Bak / (A_pp/A_PbPb)_Gwak`.\n";
  out << "`acc_eff_ratio` is the same correction ratio including efficiency.\n\n";
  out << "| obs | state | bin | Bak RAA | Gwak header | Bak/Gwak RAA | RAA diff % | acc-only ratio | acc+eff ratio | residual after acc % | residual after acc+eff % |\n";
  out << "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";
  for (const Row &r : rows)
  {
    out << "| " << r.observable << " | " << r.state << " | " << r.bin << " | "
        << F(r.bakRaa, 6) << " | " << F(r.gwakRaa, 6) << " | " << F(r.raaRatio, 6) << " | "
        << F(r.raaDiffPct, 3) << " | " << F(r.accOnlyRatio, 6) << " | " << F(r.accEffRatio, 6) << " | "
        << F(r.residualAfterAccPct, 3) << " | " << F(r.residualAfterAccEffPct, 3) << " |\n";
  }
}
} // namespace

void compare_raa_acc_eff_260501()
{
  std::vector<Row> rows;
  PushPtRows(rows, "mid_pt", "RAA_JPsi_midRap_pT.root", jpsi_raa::kPtMidBinEdges,
             jpsi_raa::kTnPL2L3PtMidPr, jpsi_raa::kTnPL2L3PtMidNp, jpsi_raa::kNPtMid, true);
  PushPtRows(rows, "fwd_pt", "RAA_JPsi_forRap_pT.root", jpsi_raa::kPtFwdBinEdges,
             jpsi_raa::kTnPL2L3PtFwdPr, jpsi_raa::kTnPL2L3PtFwdNp, jpsi_raa::kNPtFwd, false);
  PushCentRows(rows, "mid_cent", "RAA_JPsi_midRap_Npart.root", jpsi_raa::kMidNpart,
               jpsi_raa::kTnPL2L3CentMidPr, jpsi_raa::kTnPL2L3CentMidNp, jpsi_raa::kNCentMid, true);
  PushCentRows(rows, "fwd_cent", "RAA_JPsi_forRap_Npart_4Bins.root", jpsi_raa::kFwdNpart,
               jpsi_raa::kTnPL2L3CentFwdPr, jpsi_raa::kTnPL2L3CentFwdNp, jpsi_raa::kNCentFwd, false);

  WriteCsv(rows);
  WriteMarkdown(rows);

  std::cout << "wrote " << kOutCsv << "\n";
  std::cout << "wrote " << kOutMd << "\n\n";
  std::cout << "Key columns:\n";
  std::cout << "  Bak/Gwak RAA: direct comparison between 260430 ROOT RAA and jpsi_raa_values.h\n";
  std::cout << "  acc-only ratio: predicted Bak/Gwak shift from acceptance only\n";
  std::cout << "  residual after acc %: remaining percentage after dividing out acc-only ratio\n";
}
