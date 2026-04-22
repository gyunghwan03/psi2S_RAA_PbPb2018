#include "plot_compare_acc_Bak_vs_Gwak.C"

bool compare_acceptance_jpsi(bool isPbPb = true,
                             bool isPrompt = true,
                             const TString &studyRemark = "",
                             int acceptanceGenWeight = 1,
                             int acceptancePtWeight = 1,
                             int studyWeightOption = 1,
                             int studyPtWeightSystematic = 0,
                             const TString &outputTag = "",
                             const TString &acceptanceFilePath = "",
                             const TString &studyFilePath = "",
                             const TString &acceptanceLabel = "acceptance_1d",
                             const TString &studyLabel = "JpsiaccStudy")
{
  if (!outputTag.IsNull() ||
      !acceptanceLabel.EqualTo("acceptance_1d") ||
      !studyLabel.EqualTo("JpsiaccStudy"))
  {
    std::cout << "[INFO] compatibility wrapper ignores outputTag and custom legend labels.\n";
  }

  return plot_compare_acc_Bak_vs_Gwak(isPbPb,
                                      isPrompt,
                                      studyRemark.Data(),
                                      acceptanceGenWeight,
                                      acceptancePtWeight,
                                      studyWeightOption,
                                      studyPtWeightSystematic,
                                      studyFilePath.Data(),
                                      acceptanceFilePath.Data(),
                                      "");
}

void run_compare_acc_jpsi_pbpb_prompt(const TString &studyRemark = "20260420")
{
  plot_compare_acc_pbpb_prompt_Bak_vs_Gwak(studyRemark.Data());
}

void run_compare_acc_jpsi_pp_prompt(const TString &studyRemark = "20260419")
{
  plot_compare_acc_pp_prompt_Bak_vs_Gwak(studyRemark.Data());
}

void run_compare_acc_jpsi_pbpb_nonprompt(const TString &studyRemark,
                                         int acceptanceGenWeight = 1,
                                         int acceptancePtWeight = 1,
                                         int studyWeightOption = 1,
                                         int studyPtWeightSystematic = 0)
{
  plot_compare_acc_pbpb_nonprompt_Bak_vs_Gwak(studyRemark.Data(),
                                              acceptanceGenWeight,
                                              acceptancePtWeight,
                                              studyWeightOption,
                                              studyPtWeightSystematic);
}

void run_compare_acc_jpsi_pp_nonprompt(const TString &studyRemark,
                                       int acceptanceGenWeight = 1,
                                       int acceptancePtWeight = 1,
                                       int studyWeightOption = 1,
                                       int studyPtWeightSystematic = 0)
{
  plot_compare_acc_pp_nonprompt_Bak_vs_Gwak(studyRemark.Data(),
                                            acceptanceGenWeight,
                                            acceptancePtWeight,
                                            studyWeightOption,
                                            studyPtWeightSystematic);
}
