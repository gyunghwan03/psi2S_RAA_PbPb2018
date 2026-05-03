#ifndef RAA_COMPARISON_JPSI_RAA_VALUES_H
#define RAA_COMPARISON_JPSI_RAA_VALUES_H

namespace jpsi_raa
{
// Active Run-2 comparison points used in the plots.
// The current plotting macros draw only the TnPL2L3 set against HIN-16-025.

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
    0.340, 0.309, 0.328,
    0.350, 0.411, 0.450};
constexpr double kTnPL2L3PtMidPrStat[kNPtMid] = {
    0.0024, 0.0021, 0.0034,
    0.0050, 0.0108, 0.0159};

constexpr double kTnPL2L3PtMidNp[kNPtMid] = {
    0.435, 0.373, 0.365,
    0.386, 0.403, 0.439};
constexpr double kTnPL2L3PtMidNpStat[kNPtMid] = {
    0.0048, 0.0034, 0.0046,
    0.0057, 0.0101, 0.0132};

constexpr double kTnPL2L3PtFwdPr[kNPtFwd] = {
    0.428, 0.294, 0.299, 0.334};
constexpr double kTnPL2L3PtFwdPrStat[kNPtFwd] = {
    0.0033, 0.0023, 0.0032, 0.0048};

constexpr double kTnPL2L3PtFwdNp[kNPtFwd] = {
    0.527, 0.374, 0.343, 0.358};
constexpr double kTnPL2L3PtFwdNpStat[kNPtFwd] = {
    0.0077, 0.0046, 0.0052, 0.0060};

// TnPL2L3 prompt/nonprompt RAA vs centrality.
constexpr double kTnPL2L3CentMidPr[kNCentMid] = {
    0.264, 0.343, 0.420,
    0.503, 0.577, 0.664};
constexpr double kTnPL2L3CentMidPrStat[kNCentMid] = {
    0.0029, 0.0035, 0.0043,
    0.0057, 0.0081, 0.0087};

constexpr double kTnPL2L3CentMidNp[kNCentMid] = {
    0.430, 0.508, 0.562,
    0.646, 0.714, 0.763};
constexpr double kTnPL2L3CentMidNpStat[kNCentMid] = {
    0.0042, 0.0049, 0.0058,
    0.0076, 0.0102, 0.0109};

constexpr double kTnPL2L3CentFwdPr[kNCentFwd] = {
    0.368, 0.463, 0.597, 0.742};
constexpr double kTnPL2L3CentFwdPrStat[kNCentFwd] = {
    0.0040, 0.0028, 0.0038, 0.0059};

constexpr double kTnPL2L3CentFwdNp[kNCentFwd] = {
    0.508, 0.590, 0.709, 0.783};
constexpr double kTnPL2L3CentFwdNpStat[kNCentFwd] = {
    0.0103, 0.0069, 0.0092, 0.0146};

inline double pt_center(const double *edges, int idx)
{
  return 0.5 * (edges[idx] + edges[idx + 1]);
}

inline double pt_half_width(const double *edges, int idx)
{
  return 0.5 * (edges[idx + 1] - edges[idx]);
}
} // namespace jpsi_raa

#endif
