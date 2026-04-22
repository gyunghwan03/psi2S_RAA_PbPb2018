#include <iostream>
#include <vector>
#include <set>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

using namespace std;

void run_all_eff_ctauCut_v5(bool runDuplicateJpsiPP = false)
{
  vector<TString> macroFiles = {
    "get_Eff_JPsi_pbpb_ctauCut_v5.C",
    "get_Eff_Jpsi_pp_ctauCut_v4.C",
    "get_Eff_psi_pbpb_ctauCut_v4.C"
  };

  if (runDuplicateJpsiPP) {
    macroFiles.push_back("get_Eff_Jpsi_pp_ctauCut_v4.C");
  }

  set<TString> loaded;
  for (const auto &macro : macroFiles) {
    if (loaded.find(macro) != loaded.end()) continue;
    gROOT->ProcessLine(Form(".L %s", macro.Data()));
    loaded.insert(macro);
  }

  for (int state : {1, 2}) {
    cout << "[run_all_eff_ctauCut_v5] Running state = " << state << endl;

    gROOT->ProcessLine(Form("get_Eff_JPsi_pbpb_ctauCut_v5(%d);", state));
    gROOT->ProcessLine(Form("get_Eff_Jpsi_pp_ctauCut_v4(%d);", state));
    gROOT->ProcessLine(Form("get_Eff_psi_pbpb_ctauCut_v4(%d);", state));

    if (runDuplicateJpsiPP) {
      gROOT->ProcessLine(Form("get_Eff_Jpsi_pp_ctauCut_v4(%d);", state));
    }
  }

  cout << "[run_all_eff_ctauCut_v5] Done." << endl;
}

void run_all_eff_ctauCut_v5_parallel(bool runDuplicateJpsiPP = false)
{
  TString workDir = gSystem->pwd();

  vector<TString> jobs = {
    "get_Eff_JPsi_pbpb_ctauCut_v5.C",
    "get_Eff_Jpsi_pp_ctauCut_v4.C",
    "get_Eff_psi_pbpb_ctauCut_v4.C",
    "get_Eff_psi_pp_ctauCut_v4.C"
  };

  if (runDuplicateJpsiPP) {
    jobs.push_back("get_Eff_Jpsi_pp_ctauCut_v4.C");
  }

  TString shellCmd = Form("cd %s && ( ", workDir.Data());
  for (const auto &macro : jobs) {
    for (int state : {1, 2}) {
      TString macroBase = macro;
      macroBase.ReplaceAll(".C", "");
      TString logName = Form("log_%s_state%d.txt", macroBase.Data(), state);

      shellCmd += Form("root -l -b -q '%s(%d)' > %s 2>&1 & ",
                       macro.Data(), state, logName.Data());
      cout << "[run_all_eff_ctauCut_v5_parallel] Queued: " << macro << "(" << state << ") -> " << logName << endl;
    }
  }
  shellCmd += "wait )";

  cout << "[run_all_eff_ctauCut_v5_parallel] Running all jobs in parallel..." << endl;
  int exitCode = gSystem->Exec(shellCmd.Data());

  if (exitCode == 0) {
    cout << "[run_all_eff_ctauCut_v5_parallel] All jobs finished." << endl;
  } else {
    cout << "[run_all_eff_ctauCut_v5_parallel] Finished with non-zero exit code: " << exitCode << endl;
  }
}
