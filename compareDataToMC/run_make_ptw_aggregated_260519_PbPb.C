#include "make_ptw_aggregated_260519.C"

void run_make_ptw_aggregated_260519_PbPb()
{
  using namespace PtWAgg260519;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  McSample mcAaPr{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_Prompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};
  McSample mcAaNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_NonPrompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};

  TString outDir = Form("%s/compareDataToMC/%s/%s", kRepo, kOutBase, kCandTag);
  if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir)) {
    std::cerr << "[err] cannot mkdir " << outDir << "\n";
    return;
  }

  std::cout << "\n##### G_aggregateHighPt_2exp PbPb only #####\n";
  ProcessComp(true, true, false, mcAaPr, outDir);
  ProcessComp(true, true, true, mcAaPr, outDir);
  ProcessComp(true, false, false, mcAaNp, outDir);
  ProcessComp(true, false, true, mcAaNp, outDir);
}
