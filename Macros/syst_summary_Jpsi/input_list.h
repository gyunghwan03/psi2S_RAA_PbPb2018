#pragma once
#include <iostream>
#include <string>
using std::vector;
using std::string;

// Make input list
vector<string> pp_mid_pt = {
    "2DFitResult_pt6.5-9.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt9.0-12.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt12.0-15.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt15.0-20.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt20.0-25.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt25.0-40.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    //"2DFitResult_pt30.0-50.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root"
};
vector<string> pb_mid_pt = {
    "2DFitResult_pt6.5-9.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt9.0-12.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt12.0-15.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt15.0-20.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt20.0-25.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt25.0-40.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    //"2DFitResult_pt30.0-50.0_y0.0-1.6_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root"
};

vector<string> pp_mid_cent = {
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root"
};
vector<string> pb_mid_cent = {
    //"2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality0-10_PRw_Effw0_Accw0_PtW0_TnP0.root",
    //"2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality10-20_PRw_Effw0_Accw0_PtW0_TnP0.root",
    //"2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality20-30_PRw_Effw0_Accw0_PtW0_TnP0.root",
    //"2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality30-40_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality0-20_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality20-40_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality40-60_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality60-80_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality80-100_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-40.0_y0.0-1.6_muPt0.0_centrality100-180_PRw_Effw0_Accw0_PtW0_TnP0.root"
};


vector<string> pp_fwd_pt = {
    "2DFitResult_pt3.5-6.5_y1.6-2.4_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-9.0_y1.6-2.4_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt9.0-12.0_y1.6-2.4_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt12.0-40.0_y1.6-2.4_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root"
};
vector<string> pb_fwd_pt = {
    "2DFitResult_pt3.5-6.5_y1.6-2.4_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt6.5-9.0_y1.6-2.4_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt9.0-12.0_y1.6-2.4_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt12.0-40.0_y1.6-2.4_muPt0.0_centrality0-180_PRw_Effw0_Accw0_PtW0_TnP0.root"
};


vector<string> pp_fwd_cent = {
    "2DFitResult_pt3.5-40.0_y1.6-2.4_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root"
};
vector<string> pb_fwd_cent = {
    "2DFitResult_pt3.5-40.0_y1.6-2.4_muPt0.0_centrality0-20_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt3.5-40.0_y1.6-2.4_muPt0.0_centrality20-60_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt3.5-40.0_y1.6-2.4_muPt0.0_centrality60-100_PRw_Effw0_Accw0_PtW0_TnP0.root",
    "2DFitResult_pt3.5-40.0_y1.6-2.4_muPt0.0_centrality100-180_PRw_Effw0_Accw0_PtW0_TnP0.root"
};
