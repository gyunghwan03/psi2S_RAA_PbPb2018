#ifndef RAA_COMPARISON_JPSI_RAA_VALUES_SAME_ACC_H
#define RAA_COMPARISON_JPSI_RAA_VALUES_SAME_ACC_H

namespace jpsi_raa_same_acc
{
// Active same-acceptance Run-2 comparison points used in the plots.
// The plotting macros overlay these points on top of the default Run-2 set.

constexpr int kNPtMid = 6;
constexpr int kNPtFwd = 4;
constexpr int kNCentMid = 6;
constexpr int kNCentFwd = 4;

// pT binning used by the current points.
constexpr double kPtMidBinEdges[kNPtMid + 1] = {6.5, 9.0, 12.0, 15.0, 20.0, 25.0, 40.0};
constexpr double kPtFwdBinEdges[kNPtFwd + 1] = {3.5, 6.5, 9.0, 12.0, 40.0};

// <Npart> points used for the current centrality plots.
constexpr double kMidNpart[kNCentMid] = {27.12, 87.19, 131.0, 188.2, 262.3, 356.9};
constexpr double kFwdNpart[kNCentFwd] = {27.12, 109.1, 225.2, 356.9};
constexpr double kMidNpartXErr[kNCentMid] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
constexpr double kFwdNpartXErr[kNCentFwd] = {4.3, 4.3, 4.3, 4.3};

// TnPL2L3 prompt/nonprompt RAA vs pT.
constexpr double kTnPL2L3PtMidPr[kNPtMid] = {
    0.317, 0.304, 0.326,
    0.348, 0.410, 0.448};
constexpr double kTnPL2L3PtMidPrStat[kNPtMid] = {
    0.0024, 0.0021, 0.0034,
    0.0050, 0.0108, 0.0159};

constexpr double kTnPL2L3PtMidNp[kNPtMid] = {
    0.406, 0.366, 0.363,
    0.384, 0.402, 0.437};
constexpr double kTnPL2L3PtMidNpStat[kNPtMid] = {
    0.0048, 0.0034, 0.0046,
    0.0057, 0.0101, 0.0132};

constexpr double kTnPL2L3PtFwdPr[kNPtFwd] = {
    0.372, 0.289, 0.297, 0.329};
constexpr double kTnPL2L3PtFwdPrStat[kNPtFwd] = {
    0.0033, 0.0023, 0.0032, 0.0048};

constexpr double kTnPL2L3PtFwdNp[kNPtFwd] = {
    0.462, 0.367, 0.340, 0.351};
constexpr double kTnPL2L3PtFwdNpStat[kNPtFwd] = {
    0.0077, 0.0046, 0.0052, 0.0060};

// Inclusive same-acceptance 0-90% points, kept here for completeness.
constexpr double kTnPL2L3PtMidInclusivePr = 0.337;
constexpr double kTnPL2L3PtMidInclusivePrStat = 0.0018;
constexpr double kTnPL2L3PtMidInclusiveNp = 0.459;
constexpr double kTnPL2L3PtMidInclusiveNpStat = 0.0025;
constexpr double kTnPL2L3PtFwdInclusivePr = 0.362;
constexpr double kTnPL2L3PtFwdInclusivePrStat = 0.0015;
constexpr double kTnPL2L3PtFwdInclusiveNp = 0.488;
constexpr double kTnPL2L3PtFwdInclusiveNpStat = 0.0042;

// TnPL2L3 prompt/nonprompt RAA vs centrality.
constexpr double kTnPL2L3CentMidPr[kNCentMid] = {
    0.218, 0.284, 0.347,
    0.416, 0.477, 0.549};
constexpr double kTnPL2L3CentMidPrStat[kNCentMid] = {
    0.0029, 0.0035, 0.0043,
    0.0057, 0.0081, 0.0087};

constexpr double kTnPL2L3CentMidNp[kNCentMid] = {
    0.344, 0.407, 0.451,
    0.518, 0.572, 0.612};
constexpr double kTnPL2L3CentMidNpStat[kNCentMid] = {
    0.0042, 0.0049, 0.0058,
    0.0076, 0.0102, 0.0109};

constexpr double kTnPL2L3CentFwdPr[kNCentFwd] = {
    0.278, 0.350, 0.451, 0.560};
constexpr double kTnPL2L3CentFwdPrStat[kNCentFwd] = {
    0.0040, 0.0028, 0.0038, 0.0059};

constexpr double kTnPL2L3CentFwdNp[kNCentFwd] = {
    0.377, 0.437, 0.525, 0.580};
constexpr double kTnPL2L3CentFwdNpStat[kNCentFwd] = {
    0.0103, 0.0069, 0.0092, 0.0146};

// Inclusive same-acceptance 0-90% points, kept here for completeness.
constexpr double kTnPL2L3CentMidInclusivePr = 0.337;
constexpr double kTnPL2L3CentMidInclusivePrStat = 0.0018;
constexpr double kTnPL2L3CentMidInclusiveNp = 0.459;
constexpr double kTnPL2L3CentMidInclusiveNpStat = 0.0025;
constexpr double kTnPL2L3CentFwdInclusivePr = 0.362;
constexpr double kTnPL2L3CentFwdInclusivePrStat = 0.0015;
constexpr double kTnPL2L3CentFwdInclusiveNp = 0.488;
constexpr double kTnPL2L3CentFwdInclusiveNpStat = 0.0042;

inline double pt_center(const double *edges, int idx)
{
  return 0.5 * (edges[idx] + edges[idx + 1]);
}

inline double pt_half_width(const double *edges, int idx)
{
  return 0.5 * (edges[idx + 1] - edges[idx]);
}
} // namespace jpsi_raa_same_acc

#endif
