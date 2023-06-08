#ifndef TreeSetting_h
#define TreeSetting_h

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"

using namespace std;

enum mode
{
  TAG,
  STATUP,
  STATDOWN,
  SYSUP,
  SYSDOWN
};
std::map<int, std::string> mode_str =
{
  {TAG, "TagChange"},
  {STATUP, "statUp"},
  {STATDOWN, "statDo"},
  {SYSUP, "sysUp"},
  {SYSDOWN, "sysDo"}
};

#if !defined(TreeSetting_h)
extern const long int maxBranchSize = 100000;
#endif

// import the tree to the RooDataSet
int		 cBin;
int		 nTrkHits1;
int		 nTrkHits2;
int		 normChi2_global1;
int		 normChi2_global2;
int		 nMuValHits1;
int		 nMuValHits2;
int		 StationsMatched1;
int		 StationsMatched2;
int		 nPixWMea1;
int		 nPixWMea2;
int		 nTrkWMea1;
int		 nTrkWMea2;
int		 ctau;
int		 ctau3D;
int		 highPurity1;
int		 highPurity2;
double		 mass;
double		 pt;
double		 y;
double		 pt1;
double		 pt2;
double		 eta1;
double		 eta2;
double		 eta;
double		 dxy1;
double		 dxy2;
double		 dz1;
double		 dz2;
double		 dxyErr1;
double		 dxyErr2;
double		 dzErr1;
double		 dzErr2;
double		 QQVtxProb;
double		 QQMassErr;
double		 normChi2_inner1;
double		 normChi2_inner2;
double		 ptErr_inner1;
double		 ptErr_inner2;
double		 weight;
double		 tnp_weight;
double		 pt_weight;

//  muon id 
/////////////////////////////////////////
////// Gen QQ 
/////////////////////////////////////////
// NO Gen Info
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*Conversion~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//

class SetTree
{
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* tree)
{
  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("cBin", &cBin);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("pt1", &pt1);
  tree->SetBranchAddress("pt2", &pt2);
  tree->SetBranchAddress("eta1", &eta1);
  tree->SetBranchAddress("eta2", &eta2);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("tnp_weight", &tnp_weight);
  tree->SetBranchAddress("pt_weight", &pt_weight);
  tree->SetBranchAddress("QQVtxProb", &QQVtxProb);
  tree->SetBranchAddress("nTrkWMea1", &nTrkWMea1);
  tree->SetBranchAddress("nTrkWMea2", &nTrkWMea2);
  tree->SetBranchAddress("nPixWMea2", &nPixWMea2);
  tree->SetBranchAddress("nPixWMea1", &nPixWMea1);
  tree->SetBranchAddress("dxy1", &dxy1);
  tree->SetBranchAddress("dxy2", &dxy2);
  tree->SetBranchAddress("dz1", &dz1);
  tree->SetBranchAddress("dz2", &dz2);

};

class SetTree_SYSREF
{
  public:
    SetTree_SYSREF(){};

    virtual ~SetTree_SYSREF();
    virtual void TreeSetting(TTree* tree);
};

SetTree_SYSREF::~SetTree_SYSREF()
{
}


void SetTree_SYSREF::TreeSetting(TTree* tree)
{
  tree->SetBranchStatus("weight", 0);
  tree->SetBranchStatus("tnp_weight", 0);
  tree->SetBranchStatus("pt_weight", 0);
  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("cBin", &cBin);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("pt1", &pt1);
  tree->SetBranchAddress("pt2", &pt2);
  tree->SetBranchAddress("eta1", &eta1);
  tree->SetBranchAddress("eta2", &eta2);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("QQVtxProb", &QQVtxProb);
  tree->SetBranchAddress("nTrkWMea1", &nTrkWMea1);
  tree->SetBranchAddress("nTrkWMea2", &nTrkWMea2);
  tree->SetBranchAddress("nPixWMea2", &nPixWMea2);
  tree->SetBranchAddress("nPixWMea1", &nPixWMea1);
  tree->SetBranchAddress("dxy1", &dxy1);
  tree->SetBranchAddress("dxy2", &dxy2);
  tree->SetBranchAddress("dz1", &dz1);
  tree->SetBranchAddress("dz2", &dz2);
};

class SetTreeSYS
{
  public:
    SetTreeSYS(){};

    virtual ~SetTreeSYS();
    virtual void TreeSetting(TTree* tree);
    double i_weight, i_tnp_weight, i_pt_weight;
};

SetTreeSYS::~SetTreeSYS()
{
}


void SetTreeSYS::TreeSetting(TTree* tree)
{
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("weight", 1);
  tree->SetBranchStatus("tnp_weight", 1);
  tree->SetBranchStatus("pt_weight", 1);
  tree->SetBranchAddress("weight", &i_weight);
  tree->SetBranchAddress("tnp_weight", &i_tnp_weight);
  tree->SetBranchAddress("pt_weight", &i_pt_weight);

};


#endif 
