#include "make_ptw_NP_rational_260520.C"

void run_make_ptw_NP_rational_260520_PbPb()
{
  using namespace PtWNPRational260520;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  McSample mcAaNp{TString(kRepo) + "/skimmedFiles/OniaFlowSkim_JpsiTrig_NonPrompt_miniAOD_JPsi_isMC1_HFNom_240530.root",
                  2.6, 3.5, 0.0, 2.4, {}, false};

  TString outDir = Form("%s/compareDataToMC/%s/%s", kRepo, kOutBase, kCandTag);
  if (gSystem->mkdir(outDir, true) != 0 && gSystem->AccessPathName(outDir)) {
    std::cerr << "[err] cannot mkdir " << outDir << "\n";
    return;
  }

  std::cout << "\n##### H_NP_rational PbPb only #####\n";
  ProcessNP(true, false, mcAaNp, outDir);
  ProcessNP(true, true, mcAaNp, outDir);
}
