#include "make_ptw_aggregated_260519.C"

void run_make_ptw_aggregated_260519_pp()
{
  using namespace PtWAgg260519;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  McSample mcPpPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_pp_Jpsi_isMC1_241011.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcPpNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_pp_BtoJpsi_isMC1_miniAOD_251103.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};

  TString outDir = Form("%s/compareDataToMC/%s/%s", kRepo, kOutBase, kCandTag);
  if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir)) {
    std::cerr << "[err] cannot mkdir " << outDir << "\n";
    return;
  }

  std::cout << "\n##### G_aggregateHighPt_2exp pp only #####\n";
  ProcessComp(false, true, false, mcPpPr, outDir);
  ProcessComp(false, true, true, mcPpPr, outDir);
  ProcessComp(false, false, false, mcPpNp, outDir);
  ProcessComp(false, false, true, mcPpNp, outDir);
}
