#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
// #include "Style_jaebeom.h"
// #include "tnp_weight_lowptPbPb.h"
#include "../tnp_weight_lowptPbPb_num_den_new.h"
#include <TAttMarker.h>

using namespace std;
// v3 : change systematics up and down
// v4 : change TnP part to makeMuonSkimTree_Jpsi.C

void decay_length_OniaTree_v2_1S(
		int state = 2,
		int type = 1, // 1 : pt dependent ctau cut, 2 : cent dependent ctau cut, 3 : whole pt and cent dependent ctau cut
		bool isTnP = true, int isPtWeight = 0, // 0 : nominal , 1 : up, -1 : down
		bool isMC = true,
		// bool isTnP = true, bool isPtWeight = false,
		float ptLow = 3.0, float ptHigh = 40.0,
		float yLow = 0.0, float yHigh = 2.4,
		float cLow = 0, float cHigh = 180
		// bool isTnP = false, bool isPtWeight = false, int state=1
		// bool isTnP = true, bool isPtWeight = false, int state=1
		// bool isTnP = false, bool isPtWeight = true, int state=1
		)
{

	// const int nCores = 4;
	// ROOT::EnableImplicitMT();
	// ROOT::TProcessExecutor mpe(nCores);

	gStyle->SetOptStat(0);
	int kTrigSel = 12; // jpsi=12,Upsilon=13

	float muPtCut = 0; // 3.5, 1.8

	TString ptSys;
	if (isPtWeight == 0)
		ptSys = "nomi";
	else if (isPtWeight == 1)
		ptSys = "up";
	else if (isPtWeight == -1)
		ptSys = "down";
	// jpsi
	//float massLow = 0.0;
	//float massHigh = 10.0;
	float massLow = 2.6;
	float massHigh = 3.5;

	double min = 0;
	double max = ptHigh;
	const int numBins = 9; // 50;//(max-min)/binwidth;  //31//


	TString kineLabel = Form("state%d_pt%.1f-%.1f_y%.1f-%.1f_cent%.0f-%.0f", 
			state, ptLow, ptHigh, yLow, yHigh, cLow, cHigh);

	float pos_x_mass = 0.55;
	float pos_y = 0.65;
	float pos_y_diff = 0.071;
	int text_color = 1;
	float text_size = 16;
	TString perc = "%";

	// input files
	// PbPb
	TString fPRMC_path = "/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
	TString fNPMC_path = "/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root"; // PbPb_non prompt

	// Separate TChains for PRMC and NPMC (for ctau3D cut calculation)
	TChain *treePRMC = new TChain("hionia/myTree");
	TChain *treeNPMC = new TChain("hionia/myTree");
	treePRMC->Add(fPRMC_path.Data());
	treeNPMC->Add(fNPMC_path.Data());

	// Combined tree for main analysis
	TChain *mytree = new TChain("hionia/myTree");
	mytree->Add(fPRMC_path.Data());
	if(state==2) mytree->Add(fNPMC_path.Data());
	// TFile *inf = new TFile("../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root","READ");
	// TTree *mytree = (TTree*)inf->Get("myTree");
	// if(state==1) mytree->Add(inputMC2.Data());
	cout << "dmoon chk entries : " << mytree->GetEntries() << endl;
	cout << "PRMC entries : " << treePRMC->GetEntries() << endl;
	cout << "NPMC entries : " << treeNPMC->GetEntries() << endl;
	// TString outFileName = "mc_eff_vs_pt_cent_prompt_pbpb_Jpsi.root";
	// if(state==2) outFileName = "mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi.root";

	// pT reweighting function
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y1p6_2p4_241016.root", "read");
	// TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_JPsi_DATA_y0_1p6_241008.root", "read");
	// ratioDataMC_AA_Jpsi_DATA_y0.0-2.4_210915.root
	// TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","read");
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_250116.root","read"); //v2
	// before
	// TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_241126.root", "read"); // v4

	// affter 250121
	// TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_All_y_211218.root", "read");
	// TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_y1p6_2p4_241126.root", "read"); // v4
	// Jpsi v2 weighting function
	//TFile *fPtW1 = new TFile("./ratioDataMC_AA_Jpsi_DATA_All_y_211218.root", "read");
	//TFile *fPtW2 = new TFile("./ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root", "read");
	// TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Psi2S_DATA_y0_1p6_230521.root","read");
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_PbPb_Psi2S_DATA_y1p6_2p4_230520.root","read");
	// TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y0_1p6_230321.root", "read");
	// TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_Psi2S_DATA_y1p6_2p4_230621.root", "read");
	// pp Jpsi weight Fucntion
	TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
	TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root", "read");

	if( state == 2 ) {
		TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
		TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root", "read");
	}
	// if (state == 2) {
	//  TFile *fPtW1 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoPsi2S_DATA_y0_1p6_240515.root", "read");
	//  TFile *fPtW2 = new TFile("../compareDataToMC/ratioDataMC_AA_BtoPsi2S_DATA_y1p6_2p4_240515.root", "read");
	//}
	TF1 *fptw1 = (TF1 *)fPtW1->Get("dataMC_Ratio1");
	TF1 *fptw2 = (TF1 *)fPtW2->Get("dataMC_Ratio1");

	// double ptBin_all[16] = {0,6.5,7.5,8.5,9.5,11,13,15,17.5,20,22.5,25,27.5,30,40,50};
	// double ptBin_all[17] = {0,3.,4.5,5.5,6.5,7.5,8.5,9.5,11,13,15,17.5,20,22.5,25,27.5,50};
	double ptBin_for[6] = {0, 3.5, 6.5, 9, 12, 40};
	double ptBin_mid[8] = {0, 6.5, 9, 12, 15, 20, 25, 40};
	double centBin_for[5] = {0, 20, 60, 100, 180};
	double centBin_mid[7] = {0, 20, 40, 60, 80, 100, 180};

	float nbins = 4000;
	float xmin = -1; float xmax = 3;

	TH1D *hpt_ctau_1 = new TH1D("hpt_ctau_1", "hpt_ctau_1", 5, ptBin_for);
	TH1D *hpt_ctau_2 = new TH1D("hpt_ctau_2", "hpt_ctau_2", 7, ptBin_mid);


	TH1D *hcent_ctau_1 = new TH1D("hcent_ctau_1", "hcent_ctau_1", 4, centBin_for);
	TH1D *hcent_ctau_2 = new TH1D("hcent_ctau_2", "hcent_ctau_2", 6, centBin_mid);

	TH1D *hInt_ctau_1 = new TH1D("hInt_ctau_1", "", 1, 0, 50);
	TH1D *hInt_ctau_2 = new TH1D("hInt_ctau_2", "", 1, 0, 50);

	TH1D* h_decayPRMC = new TH1D("h_decayPRMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
	TH1D* h_decayNPMC = new TH1D("h_decayNPMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
	TH1D* h_deffPRMC = new TH1D("h_deffPRMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);
	TH1D* h_deffNPMC = new TH1D("h_deffNPMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);

	// For lcutv calculation
	TLine *lcutv = nullptr;
	TLine *lcuth = nullptr;
	TLine *lresi = nullptr;
	double ctau3D_cut = 0.0;
	double pr_eff = 0.0;
	double np_res = 0.0;

	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// for systematics  ////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	// parameters for all
	double p10 = 0.0, p11 = 0.0, p12 = 0.0, p13 = 0.0;
	// parameters for forward
	double p20 = 0.0, p21 = 0.0, p22 = 0.0, p23 = 0.0;
	// errors for all
	double e10 = 0.0, e11 = 0.0, e12 = 0.0, e13 = 0.0;
	// errors for forward
	double e20 = 0.0, e21 = 0.0, e22 = 0.0, e23 = 0.0;

	p10 = fptw1->GetParameter(0);
	p11 = fptw1->GetParameter(1);
	p12 = fptw1->GetParameter(2);
	p13 = fptw1->GetParameter(3);

	e10 = fptw1->GetParError(0);
	e11 = fptw1->GetParError(1);
	e12 = fptw1->GetParError(2);
	e13 = fptw1->GetParError(3);

	p20 = fptw2->GetParameter(0);
	p21 = fptw2->GetParameter(1);
	p22 = fptw2->GetParameter(2);
	p23 = fptw2->GetParameter(3);

	e20 = fptw2->GetParError(0);
	e21 = fptw2->GetParError(1);
	e22 = fptw2->GetParError(2);
	e23 = fptw2->GetParError(3);

	TF1 *fptw1sys = new TF1("fptw1sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 6.5, 50.0);
	TF1 *fptw2sys = new TF1("fptw2sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 3.0, 50.0);

	cout << "p10 : " << p10 << ", e10 : " << e10 << ", p20 : " << p20 << ", e20 : " << e20 << endl;

	// sys for nominal
	if (isPtWeight == 0)
	{
		p10 = p10;
		p11 = p11;
		p12 = p12;
		p13 = p13;

		p20 = p20;
		p21 = p21;
		p22 = p22;
		p23 = p23;
	}

	// sys for up
	if (isPtWeight == 1)
	{
		p10 = p10 + e10;
		p11 = p11 + e11;
		p12 = p12 + e12;
		p13 = p13 + e13;

		p20 = p20 + e20;
		p21 = p21 + e21;
		p22 = p22 + e22;
		p23 = p23 + e23;
	}

	// sys for down
	if (isPtWeight == -1)
	{
		p10 = p10 - e10;
		p11 = p11 - e11;
		p12 = p12 - e12;
		p13 = p13 - e13;

		p20 = p20 - e20;
		p21 = p21 - e21;
		p22 = p22 - e22;
		p23 = p23 - e23;
	}

	// cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

	fptw1sys->SetParameter(0, p10);
	fptw1sys->SetParameter(1, p11);
	fptw1sys->SetParameter(2, p12);
	fptw1sys->SetParameter(3, p13);

	fptw2sys->SetParameter(0, p20);
	fptw2sys->SetParameter(1, p21);
	fptw2sys->SetParameter(2, p22);
	fptw2sys->SetParameter(3, p23);

	hpt_ctau_1->Sumw2();
	hpt_ctau_2->Sumw2();

	hcent_ctau_1->Sumw2(); 
	hcent_ctau_2->Sumw2();

	const int maxBranchSize = 1000;
	Int_t Centrality;
	ULong64_t HLTriggers;
	TBranch *b_Centrality;                 //!
	TBranch *b_HLTriggers;                 //!

	Short_t Reco_QQ_size;
	Short_t Reco_mu_size;
	TClonesArray *Reco_QQ_4mom;
	TClonesArray *Reco_mu_4mom;
	ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
	ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_QQ_size]
	Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_QQ_size;                //!
	TBranch *b_Reco_mu_size;                //!
	TBranch *b_Reco_QQ_4mom;                //!
	TBranch *b_Reco_mu_4mom;                //!
	TBranch *b_Reco_QQ_trig;                //!
	TBranch *b_Reco_mu_trig;                //!
	TBranch *b_Reco_QQ_VtxProb;             //!

	Reco_QQ_4mom = 0;
	Reco_mu_4mom = 0;
	mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
	mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
	mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
	mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

	//  muon id
	Short_t Reco_QQ_mupl_idx[maxBranchSize];
	Short_t Reco_QQ_mumi_idx[maxBranchSize];
	Short_t Reco_QQ_whichGen[maxBranchSize];   //[Reco_QQ_size]
	TBranch *b_Reco_QQ_mupl_idx;
	TBranch *b_Reco_QQ_mumi_idx;
	TBranch *b_Reco_QQ_whichGen;   //!
	mytree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
	mytree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

	Float_t Reco_mu_dxy[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_dxy;             //!
	mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
	Float_t Reco_mu_dz[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_dz;             //!
	mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
	Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nTrkWMea;           //!
	mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
	Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nPixWMea;           //!
	mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
	Short_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_QQ_sign;             //!
	mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

	Int_t Reco_mu_SelectionType[maxBranchSize];
	TBranch *b_Reco_mu_SelectionType;
	mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

	Float_t Reco_QQ_ctau3D[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_QQ_ctau3D;              //!
	mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);

	Short_t Reco_mu_whichGen[maxBranchSize];
	TBranch *b_Reco_mu_whichGen;
	Float_t Gen_weight;
	TBranch *b_Gen_weight;
	mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
	mytree->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);

	TLorentzVector *JP_Reco = new TLorentzVector;
	TLorentzVector *mupl_Reco = new TLorentzVector;
	TLorentzVector *mumi_Reco = new TLorentzVector;

	double weight = 1;
	double tnp_weight = 1;
	double tnp_trig_weight = 1;
	double tnp_trig_weight_mupl = -1;
	double tnp_trig_weight_mumi = -1;
	double pt_weight = 1;

	double pt, y, eta, pt1, pt2, eta1, eta2, phi1, phi2;
	double dcosMuTheta;
	double ptimb;
	//  TH2D* hpt_tnp_trig = new TH2D("hpt_tnp_trig","hpt_tnp_trig",numBins,ptBin,40,0,2);

	int kL2filter = 16; // jpsi=16,Upsilon=38
	int kL3filter = 17; // jpsi=17,Upsilon=39

	int count = 0;
	int counttnp = 0;

	//============================================================================
	// Calculate ctau3D cut (lcutv) like psuedo_proper_decay_length_1S.C
	//============================================================================
	cout << "============================================" << endl;
	cout << "Calculating ctau3D cut from PRMC and NPMC..." << endl;
	cout << "============================================" << endl;

	// Set up branches for PRMC tree
	const int maxBranchSize_ctau = 1000;
	Int_t Centrality_pr;
	ULong64_t HLTriggers_pr;
	Short_t Reco_QQ_size_pr;
	TClonesArray *Reco_QQ_4mom_pr = 0;
	TClonesArray *Reco_mu_4mom_pr = 0;
	ULong64_t Reco_QQ_trig_pr[maxBranchSize_ctau];
	Float_t Reco_QQ_VtxProb_pr[maxBranchSize_ctau];
	Short_t Reco_QQ_mupl_idx_pr[maxBranchSize_ctau];
	Short_t Reco_QQ_mumi_idx_pr[maxBranchSize_ctau];
	Short_t Reco_QQ_sign_pr[maxBranchSize_ctau];
	Float_t Reco_QQ_ctau3D_pr[maxBranchSize_ctau];
	Int_t Reco_mu_SelectionType_pr[maxBranchSize_ctau];
	Int_t Reco_mu_nTrkWMea_pr[maxBranchSize_ctau];
	Int_t Reco_mu_nPixWMea_pr[maxBranchSize_ctau];
	Float_t Reco_mu_dxy_pr[maxBranchSize_ctau];
	Float_t Reco_mu_dz_pr[maxBranchSize_ctau];
	Short_t Reco_mu_whichGen_pr[maxBranchSize_ctau];
	Short_t Reco_QQ_whichGen_pr[maxBranchSize_ctau];
	ULong64_t Reco_mu_trig_pr[maxBranchSize_ctau];
	Float_t Gen_weight_pr;

	treePRMC->SetBranchAddress("Centrality", &Centrality_pr);
	treePRMC->SetBranchAddress("HLTriggers", &HLTriggers_pr);
	treePRMC->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size_pr);
	treePRMC->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom_pr);
	treePRMC->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom_pr);
	treePRMC->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig_pr);
	treePRMC->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb_pr);
	treePRMC->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx_pr);
	treePRMC->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx_pr);
	treePRMC->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign_pr);
	treePRMC->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D_pr);
	treePRMC->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType_pr);
	treePRMC->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea_pr);
	treePRMC->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea_pr);
	treePRMC->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy_pr);
	treePRMC->SetBranchAddress("Reco_mu_dz", Reco_mu_dz_pr);
	treePRMC->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen_pr);
	treePRMC->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen_pr);
	treePRMC->SetBranchAddress("Reco_mu_trig", Reco_mu_trig_pr);
	treePRMC->SetBranchAddress("Gen_weight", &Gen_weight_pr);

	// Same for NPMC tree
	Int_t Centrality_np;
	ULong64_t HLTriggers_np;
	Short_t Reco_QQ_size_np;
	TClonesArray *Reco_QQ_4mom_np = 0;
	TClonesArray *Reco_mu_4mom_np = 0;
	ULong64_t Reco_QQ_trig_np[maxBranchSize_ctau];
	Float_t Reco_QQ_VtxProb_np[maxBranchSize_ctau];
	Short_t Reco_QQ_mupl_idx_np[maxBranchSize_ctau];
	Short_t Reco_QQ_mumi_idx_np[maxBranchSize_ctau];
	Short_t Reco_QQ_sign_np[maxBranchSize_ctau];
	Float_t Reco_QQ_ctau3D_np[maxBranchSize_ctau];
	Int_t Reco_mu_SelectionType_np[maxBranchSize_ctau];
	Int_t Reco_mu_nTrkWMea_np[maxBranchSize_ctau];
	Int_t Reco_mu_nPixWMea_np[maxBranchSize_ctau];
	Float_t Reco_mu_dxy_np[maxBranchSize_ctau];
	Float_t Reco_mu_dz_np[maxBranchSize_ctau];
	Short_t Reco_mu_whichGen_np[maxBranchSize_ctau];
	Short_t Reco_QQ_whichGen_np[maxBranchSize_ctau];
	ULong64_t Reco_mu_trig_np[maxBranchSize_ctau];
	Float_t Gen_weight_np;

	treeNPMC->SetBranchAddress("Centrality", &Centrality_np);
	treeNPMC->SetBranchAddress("HLTriggers", &HLTriggers_np);
	treeNPMC->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size_np);
	treeNPMC->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom_np);
	treeNPMC->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom_np);
	treeNPMC->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig_np);
	treeNPMC->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb_np);
	treeNPMC->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx_np);
	treeNPMC->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx_np);
	treeNPMC->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign_np);
	treeNPMC->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D_np);
	treeNPMC->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType_np);
	treeNPMC->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea_np);
	treeNPMC->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea_np);
	treeNPMC->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy_np);
	treeNPMC->SetBranchAddress("Reco_mu_dz", Reco_mu_dz_np);
	treeNPMC->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen_np);
	treeNPMC->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen_np);
	treeNPMC->SetBranchAddress("Reco_mu_trig", Reco_mu_trig_np);
	treeNPMC->SetBranchAddress("Gen_weight", &Gen_weight_np);

	if (type == 3) {
	// Fill PRMC ctau3D histogram
	cout << "Filling PRMC ctau3D histogram..." << endl;
	for (int iev = 0; iev < treePRMC->GetEntries(); ++iev) {
		if (iev % 10000000 == 0) cout << "PRMC: " << iev << " / " << treePRMC->GetEntries() << endl;
		treePRMC->GetEntry(iev);

		bool HLTPass_pr = ((HLTriggers_pr & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

		for (Int_t irqq = 0; irqq < Reco_QQ_size_pr; irqq++) {
			TLorentzVector *JP_pr = (TLorentzVector*)Reco_QQ_4mom_pr->At(irqq);
			TLorentzVector *mupl_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mupl_idx_pr[irqq]);
			TLorentzVector *mumi_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mumi_idx_pr[irqq]);

			bool HLTFilterPass_pr = ((Reco_QQ_trig_pr[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
			if (!HLTFilterPass_pr) continue;
			if (Reco_mu_whichGen_pr[Reco_QQ_mupl_idx_pr[irqq]] == -1) continue;
			if (Reco_mu_whichGen_pr[Reco_QQ_mumi_idx_pr[irqq]] == -1) continue;
			if (Reco_QQ_whichGen_pr[irqq] == -1) continue;
			if (Reco_QQ_sign_pr[irqq] != 0) continue;
			//if (Reco_QQ_VtxProb_pr[irqq] < 0.005) continue;
			if (fabs(JP_pr->Rapidity()) > yHigh || fabs(JP_pr->Rapidity()) < yLow) continue;
			if (JP_pr->Pt() < ptLow || JP_pr->Pt() > ptHigh) continue;
			if (JP_pr->M() < massLow || JP_pr->M() > massHigh) continue;
			if (!(IsAcceptanceQQ(mupl_pr->Pt(), mupl_pr->Eta()) && IsAcceptanceQQ(mumi_pr->Pt(), mumi_pr->Eta()))) continue;
			if (Centrality_pr < cLow || Centrality_pr > cHigh) continue;

			// Soft muon cuts
			bool passMuonTypePl_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,1))) && 
				(Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,3)));
			bool passMuonTypeMi_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,1))) && 
				(Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,3)));
			if (!passMuonTypePl_pr || !passMuonTypeMi_pr) continue;

			bool muplSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 5) &&
				(Reco_mu_nPixWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 0) &&
				(fabs(Reco_mu_dxy_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 20.) && passMuonTypePl_pr;
			bool mumiSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 5) &&
				(Reco_mu_nPixWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 0) &&
				(fabs(Reco_mu_dxy_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 20.) && passMuonTypeMi_pr;
			if (!(muplSoft_pr && mumiSoft_pr)) continue;

			// L2 filter cut (same as line 713 in main loop)
			if (!((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
						(Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

			// isTnP condition: L2/L3 filter SelDone check (same as lines 715-752 in main loop)
			if (isTnP) {
				bool mupl_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
				bool mupl_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
				bool mumi_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
				bool mumi_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));

				bool mupl_isL2_pr = (mupl_L2Filter_pr && !mupl_L3Filter_pr);
				bool mupl_isL3_pr = (mupl_L2Filter_pr && mupl_L3Filter_pr);
				bool mumi_isL2_pr = (mumi_L2Filter_pr && !mumi_L3Filter_pr);
				bool mumi_isL3_pr = (mumi_L2Filter_pr && mumi_L3Filter_pr);
				bool SelDone_pr = false;

				if ((mupl_isL2_pr && mumi_isL3_pr) || (mupl_isL3_pr && mumi_isL2_pr) || (mupl_isL3_pr && mumi_isL3_pr)) {
					SelDone_pr = true;
				}
				if (!SelDone_pr) continue;
			}

			if (HLTPass_pr && HLTFilterPass_pr) {
				// Calculate weight (Ncoll * Gen_weight * TnP * pT weight)
				double weight_pr = findNcoll(Centrality_pr) * Gen_weight_pr;

				// pT weight
				double pt_weight_pr = 1.0;
				if (fabs(JP_pr->Rapidity()) < 1.6 && JP_pr->Pt() >= 6.5)
					pt_weight_pr = fptw1sys->Eval(JP_pr->Pt());
				if (fabs(JP_pr->Rapidity()) >= 1.6 && fabs(JP_pr->Rapidity()) < 2.4)
					pt_weight_pr = fptw2sys->Eval(JP_pr->Pt());

				// TnP weight
				double tnp_weight_pr = 1.0;
				if (isTnP) {
					double pt1_pr = mupl_pr->Pt();
					double pt2_pr = mumi_pr->Pt();
					double eta1_pr = mupl_pr->Eta();
					double eta2_pr = mumi_pr->Eta();
					tnp_weight_pr = tnp_weight_muid_pbpb(pt1_pr, eta1_pr, 0) * tnp_weight_muid_pbpb(pt2_pr, eta2_pr, 0);
					tnp_weight_pr *= tnp_weight_trk_pbpb(eta1_pr, 0) * tnp_weight_trk_pbpb(eta2_pr, 0);

					// Trigger TnP weight
					bool mupl_L2Filter_pr_tnp = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mupl_L3Filter_pr_tnp = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mumi_L2Filter_pr_tnp = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mumi_L3Filter_pr_tnp = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mupl_isL2_pr_tnp = (mupl_L2Filter_pr_tnp && !mupl_L3Filter_pr_tnp);
					bool mupl_isL3_pr_tnp = (mupl_L2Filter_pr_tnp && mupl_L3Filter_pr_tnp);
					bool mumi_isL2_pr_tnp = (mumi_L2Filter_pr_tnp && !mumi_L3Filter_pr_tnp);
					bool mumi_isL3_pr_tnp = (mumi_L2Filter_pr_tnp && mumi_L3Filter_pr_tnp);

					double tnp_trig_weight_mupl_pr = 1.0;
					double tnp_trig_weight_mumi_pr = 1.0;
					if (mupl_isL2_pr_tnp && mumi_isL3_pr_tnp) {
						tnp_trig_weight_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
						tnp_trig_weight_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
					} else if (mupl_isL3_pr_tnp && mumi_isL2_pr_tnp) {
						tnp_trig_weight_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
						tnp_trig_weight_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
					} else if (mupl_isL3_pr_tnp && mumi_isL3_pr_tnp) {
						int t[2] = {-1, 1};
						int l = rand() % 2;
						if (t[l] == -1) {
							tnp_trig_weight_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
							tnp_trig_weight_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
						} else {
							tnp_trig_weight_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
							tnp_trig_weight_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
						}
					}
					tnp_weight_pr *= tnp_trig_weight_mupl_pr * tnp_trig_weight_mumi_pr;
				}

				double total_weight_pr = weight_pr * pt_weight_pr * tnp_weight_pr;
				h_decayPRMC->Fill(Reco_QQ_ctau3D_pr[irqq], total_weight_pr);
			}
		}
	}

	// Fill NPMC ctau3D histogram
	cout << "Filling NPMC ctau3D histogram..." << endl;
	for (int iev = 0; iev < treeNPMC->GetEntries(); ++iev) {
		if (iev % 10000000 == 0) cout << "NPMC: " << iev << " / " << treeNPMC->GetEntries() << endl;
		treeNPMC->GetEntry(iev);

		bool HLTPass_np = ((HLTriggers_np & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

		for (Int_t irqq = 0; irqq < Reco_QQ_size_np; irqq++) {
			TLorentzVector *JP_np = (TLorentzVector*)Reco_QQ_4mom_np->At(irqq);
			TLorentzVector *mupl_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mupl_idx_np[irqq]);
			TLorentzVector *mumi_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mumi_idx_np[irqq]);

			bool HLTFilterPass_np = ((Reco_QQ_trig_np[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
			if (!HLTFilterPass_np) continue;
			if (Reco_mu_whichGen_np[Reco_QQ_mupl_idx_np[irqq]] == -1) continue;
			if (Reco_mu_whichGen_np[Reco_QQ_mumi_idx_np[irqq]] == -1) continue;
			if (Reco_QQ_whichGen_np[irqq] == -1) continue;
			if (Reco_QQ_sign_np[irqq] != 0) continue;
			//if (Reco_QQ_VtxProb_np[irqq] < 0.005) continue;
			if (fabs(JP_np->Rapidity()) > yHigh || fabs(JP_np->Rapidity()) < yLow) continue;
			if (JP_np->Pt() < ptLow || JP_np->Pt() > ptHigh) continue;
			if (JP_np->M() < massLow || JP_np->M() > massHigh) continue;
			if (!(IsAcceptanceQQ(mupl_np->Pt(), mupl_np->Eta()) && IsAcceptanceQQ(mumi_np->Pt(), mumi_np->Eta()))) continue;
			if (Centrality_np < cLow || Centrality_np > cHigh) continue;

			// Soft muon cuts
			bool passMuonTypePl_np = (Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,1))) && 
				(Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,3)));
			bool passMuonTypeMi_np = (Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,1))) && 
				(Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,3)));
			if (!passMuonTypePl_np || !passMuonTypeMi_np) continue;

			bool muplSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mupl_idx_np[irqq]] > 5) &&
				(Reco_mu_nPixWMea_np[Reco_QQ_mupl_idx_np[irqq]] > 0) &&
				(fabs(Reco_mu_dxy_np[Reco_QQ_mupl_idx_np[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz_np[Reco_QQ_mupl_idx_np[irqq]]) < 20.) && passMuonTypePl_np;
			bool mumiSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 5) &&
				(Reco_mu_nPixWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 0) &&
				(fabs(Reco_mu_dxy_np[Reco_QQ_mumi_idx_np[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz_np[Reco_QQ_mumi_idx_np[irqq]]) < 20.) && passMuonTypeMi_np;
			if (!(muplSoft_np && mumiSoft_np)) continue;

			// L2 filter cut (same as line 713 in main loop)
			if (!((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
						(Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

			// isTnP condition: L2/L3 filter SelDone check (same as lines 715-752 in main loop)
			if (isTnP) {
				bool mupl_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
				bool mupl_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
				bool mumi_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
				bool mumi_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));

				bool mupl_isL2_np = (mupl_L2Filter_np && !mupl_L3Filter_np);
				bool mupl_isL3_np = (mupl_L2Filter_np && mupl_L3Filter_np);
				bool mumi_isL2_np = (mumi_L2Filter_np && !mumi_L3Filter_np);
				bool mumi_isL3_np = (mumi_L2Filter_np && mumi_L3Filter_np);
				bool SelDone_np = false;

				if ((mupl_isL2_np && mumi_isL3_np) || (mupl_isL3_np && mumi_isL2_np) || (mupl_isL3_np && mumi_isL3_np)) {
					SelDone_np = true;
				}
				if (!SelDone_np) continue;
			}

			if (HLTPass_np && HLTFilterPass_np) {
				// Calculate weight (Ncoll * Gen_weight * TnP * pT weight)
				double weight_np = findNcoll(Centrality_np) * Gen_weight_np;

				// pT weight
				double pt_weight_np = 1.0;
				if (fabs(JP_np->Rapidity()) < 1.6 && JP_np->Pt() >= 6.5)
					pt_weight_np = fptw1sys->Eval(JP_np->Pt());
				if (fabs(JP_np->Rapidity()) >= 1.6 && fabs(JP_np->Rapidity()) < 2.4)
					pt_weight_np = fptw2sys->Eval(JP_np->Pt());

				// TnP weight
				double tnp_weight_np = 1.0;
				if (isTnP) {
					double pt1_np = mupl_np->Pt();
					double pt2_np = mumi_np->Pt();
					double eta1_np = mupl_np->Eta();
					double eta2_np = mumi_np->Eta();
					tnp_weight_np = tnp_weight_muid_pbpb(pt1_np, eta1_np, 0) * tnp_weight_muid_pbpb(pt2_np, eta2_np, 0);
					tnp_weight_np *= tnp_weight_trk_pbpb(eta1_np, 0) * tnp_weight_trk_pbpb(eta2_np, 0);

					// Trigger TnP weight
					bool mupl_L2Filter_np_tnp = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mupl_L3Filter_np_tnp = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mumi_L2Filter_np_tnp = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mumi_L3Filter_np_tnp = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mupl_isL2_np_tnp = (mupl_L2Filter_np_tnp && !mupl_L3Filter_np_tnp);
					bool mupl_isL3_np_tnp = (mupl_L2Filter_np_tnp && mupl_L3Filter_np_tnp);
					bool mumi_isL2_np_tnp = (mumi_L2Filter_np_tnp && !mumi_L3Filter_np_tnp);
					bool mumi_isL3_np_tnp = (mumi_L2Filter_np_tnp && mumi_L3Filter_np_tnp);

					double tnp_trig_weight_mupl_np = 1.0;
					double tnp_trig_weight_mumi_np = 1.0;
					if (mupl_isL2_np_tnp && mumi_isL3_np_tnp) {
						tnp_trig_weight_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
						tnp_trig_weight_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
					} else if (mupl_isL3_np_tnp && mumi_isL2_np_tnp) {
						tnp_trig_weight_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
						tnp_trig_weight_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
					} else if (mupl_isL3_np_tnp && mumi_isL3_np_tnp) {
						int t[2] = {-1, 1};
						int l = rand() % 2;
						if (t[l] == -1) {
							tnp_trig_weight_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
							tnp_trig_weight_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
						} else {
							tnp_trig_weight_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
							tnp_trig_weight_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
						}
					}
					tnp_weight_np *= tnp_trig_weight_mupl_np * tnp_trig_weight_mumi_np;
				}

				double total_weight_np = weight_np * pt_weight_np * tnp_weight_np;
				h_decayNPMC->Fill(Reco_QQ_ctau3D_np[irqq], total_weight_np);
			}
		}
	}

	cout << "PRMC entries in ctau histogram: " << h_decayPRMC->GetEntries() << endl;
	cout << "NPMC entries in ctau histogram: " << h_decayNPMC->GetEntries() << endl;
	cout << "PRMC weighted integral: " << h_decayPRMC->Integral() << endl;
	cout << "NPMC weighted integral: " << h_decayNPMC->Integral() << endl;

	// Calculate lcutv (ctau3D cut) like psuedo_proper_decay_length_1S.C
	// Use Integral() for weighted histograms instead of GetEntries()
	double totalPRMC = h_decayPRMC->Integral();
	double totalNPMC = h_decayNPMC->Integral();

	for (int bins = 0; bins < h_decayPRMC->GetNbinsX(); bins++) {
		if (state == 1) { // Prompt selection: find cut where PR efficiency ~90%
			h_deffPRMC->SetBinContent(bins, h_decayPRMC->Integral(0, bins) / totalPRMC);
			h_deffNPMC->SetBinContent(bins, 1 - (h_decayNPMC->Integral(0, bins) / totalNPMC));

			if (h_decayPRMC->Integral(0, bins) < (0.899 * totalPRMC) || h_decayPRMC->Integral(0, bins) > (totalPRMC * 0.9099)) continue;
			cout << bins << "th bin, ctau3D: " << h_decayPRMC->GetBinCenter(bins) 
				<< ", PR eff: " << h_decayPRMC->Integral(0, bins) / totalPRMC * 100 << "%" << endl;
			lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
			lcuth = new TLine(xmin, h_decayPRMC->Integral(0, bins) / totalPRMC, 
					h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0, bins) / totalPRMC);
			lresi = new TLine(xmin, 1 - (h_decayNPMC->Integral(0, bins) / totalNPMC), 
					h_decayNPMC->GetBinCenter(bins), 1 - (h_decayNPMC->Integral(0, bins) / totalNPMC));
		}
		else if (state == 2) { // Non-prompt selection: find cut where PR rejection ~97-98%
			h_deffNPMC->SetBinContent(bins, h_decayNPMC->Integral(0, bins) / totalNPMC);
			h_deffPRMC->SetBinContent(bins, 1 - (h_decayPRMC->Integral(0, bins) / totalPRMC));

			if (h_decayPRMC->Integral(0, bins) < (0.97 * totalPRMC) || h_decayPRMC->Integral(0, bins) > (totalPRMC * 0.98)) continue;
			cout << bins << "th bin, ctau3D: " << h_decayPRMC->GetBinCenter(bins) 
				<< ", NP eff: " << (1 - h_decayNPMC->Integral(0, bins) / totalNPMC) * 100 << "%" << endl;
			lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
			lcuth = new TLine(xmin, h_decayNPMC->Integral(0, bins) / totalNPMC, 
					h_decayNPMC->GetBinCenter(bins), h_decayNPMC->Integral(0, bins) / totalNPMC);
			lresi = new TLine(xmin, h_decayPRMC->Integral(0, bins) / totalPRMC, 
					h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0, bins) / totalPRMC);
		}
	}

	if (lcutv) {
		ctau3D_cut = lcutv->GetX1();
		if (state == 1) {
			pr_eff = lcuth->GetY1();
			np_res = lresi->GetY1();
			cout << "============================================" << endl;
			cout << "ctau3D cut (for Prompt): " << ctau3D_cut << endl;
			cout << "Prompt J/psi efficiency: " << pr_eff * 100 << "%" << endl;
			cout << "Non-prompt J/psi residual: " << np_res * 100 << "%" << endl;
			cout << "Apply cut: ctau3D <= " << ctau3D_cut << endl;
			cout << "============================================" << endl;
		} else if (state == 2) {
			np_res = 1 - lcuth->GetY1();  // NP efficiency
			pr_eff = 1 - lresi->GetY1();  // PR residual
			cout << "============================================" << endl;
			cout << "ctau3D cut (for Non-Prompt): " << ctau3D_cut << endl;
			cout << "Non-prompt J/psi efficiency: " << np_res * 100 << "%" << endl;
			cout << "Prompt J/psi residual: " << pr_eff * 100 << "%" << endl;
			cout << "Apply cut: ctau3D >= " << ctau3D_cut << endl;
			cout << "============================================" << endl;
		}
	} else {
		cout << "ERROR: Could not find ctau3D cut!" << endl;
	}

	//============================================================================
	// Draw ctau3D cut plots (like psuedo_proper_decay_length_1S.C)
	//============================================================================
	gSystem->mkdir("./figs_ctau3D_v2/", 1);
	gSystem->mkdir("./roots_ctau3D_v2/", 1);


	// Draw ctau3D efficiency plot
	TCanvas* c_decayL = new TCanvas("c_decayL", "", 600, 600);
	h_deffPRMC->SetLineColor(kBlue);
	h_deffPRMC->SetLineWidth(2);
	h_deffPRMC->GetYaxis()->SetRangeUser(0, 1.1);
	h_deffPRMC->GetXaxis()->SetTitle("ctau3D (mm)");
	h_deffPRMC->GetYaxis()->SetTitle("Efficiency");
	h_deffPRMC->Draw("l");

	h_deffNPMC->SetLineColor(kRed+2);
	h_deffNPMC->SetLineWidth(2);
	h_deffNPMC->Draw("same l");

	// Draw reference line at y=1
	TLine *lref = new TLine(xmin, 1, xmax, 1);
	lref->SetLineStyle(2);
	lref->Draw("same");

	// Draw cut line
	if (lcutv) {
		lcutv->SetLineColor(kBlack);
		lcutv->SetLineWidth(2);
		lcutv->SetLineStyle(2);
		lcutv->Draw("same");
	}

	// Draw text labels
	if (ptLow != 0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", (float)ptLow, ptHigh), pos_x_mass, pos_y, text_color, text_size);
	if (yLow == 0) drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), pos_x_mass, pos_y - pos_y_diff * 0.5, text_color, text_size);
	else drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), pos_x_mass, pos_y - pos_y_diff, text_color, text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass, pos_y - pos_y_diff * 1.5, text_color, text_size);
	drawText(Form("Centrality %.0f-%.0f%s", cLow / 2, cHigh / 2, perc.Data()), pos_x_mass, pos_y - pos_y_diff * 2, text_color, text_size);

	if (lcutv) {
		if (state == 1) {
			drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
			drawText(Form("PR J/#psi eff: %.3f", pr_eff), pos_x_mass, pos_y - pos_y_diff * 3.5, kBlue, text_size);
			drawText(Form("NP J/#psi res: %.3f", np_res), pos_x_mass, pos_y - pos_y_diff * 4, kRed + 2, text_size);
		} else if (state == 2) {
			drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
			drawText(Form("NP J/#psi eff: %.3f", np_res), pos_x_mass, pos_y - pos_y_diff * 3.5, kRed + 2, text_size);
			drawText(Form("PR J/#psi res: %.3f", pr_eff), pos_x_mass, pos_y - pos_y_diff * 4, kBlue, text_size);
		}
	}

	// Legend
	TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	if (state == 1) {
		leg->AddEntry(h_deffPRMC, "Prompt J/#psi (CDF)", "l");
		leg->AddEntry(h_deffNPMC, "Non-Prompt J/#psi (1-CDF)", "l");
	} else {
		leg->AddEntry(h_deffNPMC, "Non-Prompt J/#psi (CDF)", "l");
		leg->AddEntry(h_deffPRMC, "Prompt J/#psi (1-CDF)", "l");
	}
	leg->Draw();

	c_decayL->Update();
	if (state == 1) {
		c_decayL->SaveAs(Form("./figs_ctau3D_v2/ctau3D_cut_PRMC_%s.pdf", kineLabel.Data()));
		c_decayL->SaveAs(Form("./figs_ctau3D_v2/ctau3D_cut_PRMC_%s.png", kineLabel.Data()));
	} else {
		c_decayL->SaveAs(Form("./figs_ctau3D_v2/ctau3D_cut_NPMC_%s.pdf", kineLabel.Data()));
		c_decayL->SaveAs(Form("./figs_ctau3D_v2/ctau3D_cut_NPMC_%s.png", kineLabel.Data()));
	}

	// Draw ctau3D distribution plot
	TCanvas* c_ctau_dist = new TCanvas("c_ctau_dist", "", 600, 600);
	c_ctau_dist->SetLogy();

	h_decayPRMC->SetLineColor(kBlue);
	h_decayPRMC->SetLineWidth(2);
	h_decayPRMC->GetXaxis()->SetTitle("ctau3D (mm)");
	h_decayPRMC->GetYaxis()->SetTitle("Counts");
	h_decayPRMC->GetXaxis()->SetRangeUser(-0.5, 1.5);
	h_decayPRMC->DrawNormalized("hist");

	h_decayNPMC->SetLineColor(kRed+2);
	h_decayNPMC->SetLineWidth(2);
	h_decayNPMC->DrawNormalized("hist same");

	// Draw cut line on distribution
	if (lcutv) {
		TLine *lcutv_dist = new TLine(ctau3D_cut, 0, ctau3D_cut, h_decayPRMC->GetMaximum() / h_decayPRMC->Integral());
		lcutv_dist->SetLineColor(kBlack);
		lcutv_dist->SetLineWidth(2);
		lcutv_dist->SetLineStyle(2);
		lcutv_dist->Draw("same");
	}

	TLegend *leg2 = new TLegend(0.55, 0.75, 0.88, 0.88);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->AddEntry(h_decayPRMC, "Prompt J/#psi MC", "l");
	leg2->AddEntry(h_decayNPMC, "Non-Prompt J/#psi MC", "l");
	leg2->Draw();

	if (ptLow != 0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", (float)ptLow, ptHigh), pos_x_mass, pos_y, text_color, text_size);
	drawText(Form("Centrality %.0f-%.0f%s", cLow / 2, cHigh / 2, perc.Data()), pos_x_mass, pos_y - pos_y_diff, text_color, text_size);

	c_ctau_dist->Update();
	if (state == 1) {
		c_ctau_dist->SaveAs(Form("./figs_ctau3D_v2/ctau3D_dist_PRMC_%s.pdf", kineLabel.Data()));
		c_ctau_dist->SaveAs(Form("./figs_ctau3D_v2/ctau3D_dist_PRMC_%s.png", kineLabel.Data()));
	} else {
		c_ctau_dist->SaveAs(Form("./figs_ctau3D_v2/ctau3D_dist_NPMC_%s.pdf", kineLabel.Data()));
		c_ctau_dist->SaveAs(Form("./figs_ctau3D_v2/ctau3D_dist_NPMC_%s.png", kineLabel.Data()));
	}

	// Save histograms and cut value to ROOT file
	TString rootFileName;
	if (state == 1) {
		rootFileName = Form("./roots_ctau3D_v2/ctau3D_cut_PRMC_%s.root", kineLabel.Data());
	} else {
		rootFileName = Form("./roots_ctau3D_v2/ctau3D_cut_NPMC_%s.root", kineLabel.Data());
	}
	TFile *wf = new TFile(rootFileName, "RECREATE");
	wf->cd();
	h_deffPRMC->Write();
	h_deffNPMC->Write();
	h_decayPRMC->Write();
	h_decayNPMC->Write();

	// Save the ctau3D cut value
	TH1D *h_ctau3D_cut = new TH1D("h_ctau3D_cut", "ctau3D cut value", 1, 0, 1);
	h_ctau3D_cut->SetBinContent(1, ctau3D_cut);
	h_ctau3D_cut->Write();

	TH1D *h_eff = new TH1D("h_eff", "efficiency", 1, 0, 1);
	if (state == 1) h_eff->SetBinContent(1, pr_eff);
	else h_eff->SetBinContent(1, np_res);
	h_eff->Write();

	TH1D *h_res = new TH1D("h_res", "residual", 1, 0, 1);
	if (state == 1) h_res->SetBinContent(1, np_res);
	else h_res->SetBinContent(1, pr_eff);
	h_res->Write();

	wf->Close();
	cout << "Saved ROOT file: " << rootFileName << endl;
}

	// Define pt bins based on rapidity
	int nPtBins_for = 5;  // for forward (1.6-2.4)
	int nPtBins_mid = 7;  // for mid (0-1.6)
	const int nCentBins_for = 4;  // centBin_for has 5 edges = 4 bins
	const int nCentBins_mid = 6;  // centBin_mid has 7 edges = 6 bins

	if(type == 1 || type == 3)
	{
		//============================================================================
		// Calculate ctau3D cut for each pt bin and draw pt vs ctau3D cut
		// OPTIMIZED: Read tree only once and fill all pt bin histograms simultaneously
		//============================================================================
		cout << endl;
		cout << "============================================" << endl;
		cout << "Calculating ctau3D cut for each pt bin (OPTIMIZED)..." << endl;
		cout << "============================================" << endl;


		// Arrays to store ctau3D cut values for each pt bin
		double ctau3D_cut_ptBin_for[5] = {0};
		double ctau3D_cut_ptBin_mid[7] = {0};
		double pr_eff_ptBin_for[5] = {0};
		double pr_eff_ptBin_mid[7] = {0};
		double np_res_ptBin_for[5] = {0};
		double np_res_ptBin_mid[7] = {0};

		// Create all histograms at once
		TH1D* h_decay_pr_pt_for[5];
		TH1D* h_decay_np_pt_for[5];
		TH1D* h_decay_pr_pt_mid[7];
		TH1D* h_decay_np_pt_mid[7];

		// Arrays to store ctau3D cut values for each cent bin
		double ctau3D_cut_centBin_for[nCentBins_for] = {0};
		double ctau3D_cut_centBin_mid[nCentBins_mid] = {0};
		double pr_eff_centBin_for[nCentBins_for] = {0};
		double pr_eff_centBin_mid[nCentBins_mid] = {0};
		double np_res_centBin_for[nCentBins_for] = {0};
		double np_res_centBin_mid[nCentBins_mid] = {0};

		// Integrated values for forward and mid rapidity
		double ctau3D_cut_int_for = 0, ctau3D_cut_int_mid = 0;
		double pr_eff_int_for = 0, pr_eff_int_mid = 0;
		double np_res_int_for = 0, np_res_int_mid = 0;

		for (int ipt = 0; ipt < nPtBins_for; ipt++) {
			h_decay_pr_pt_for[ipt] = new TH1D(Form("h_decay_pr_pt_for_%d", ipt), "", nbins, xmin, xmax);
			h_decay_np_pt_for[ipt] = new TH1D(Form("h_decay_np_pt_for_%d", ipt), "", nbins, xmin, xmax);
		}
		for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
			h_decay_pr_pt_mid[ipt] = new TH1D(Form("h_decay_pr_pt_mid_%d", ipt), "", nbins, xmin, xmax);
			h_decay_np_pt_mid[ipt] = new TH1D(Form("h_decay_np_pt_mid_%d", ipt), "", nbins, xmin, xmax);
		}

		// Fill PRMC histograms for ALL pt bins in ONE loop
		cout << "Filling PRMC histograms for all pt bins..." << endl;
		for (int iev = 0; iev < treePRMC->GetEntries(); ++iev) {
			if (iev % 10000000 == 0) cout << "PRMC pt-binned: " << iev << " / " << treePRMC->GetEntries() << endl;
			treePRMC->GetEntry(iev);
			bool HLTPass_pr = ((HLTriggers_pr & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

			for (Int_t irqq = 0; irqq < Reco_QQ_size_pr; irqq++) {
				TLorentzVector *JP_pr = (TLorentzVector*)Reco_QQ_4mom_pr->At(irqq);
				TLorentzVector *mupl_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mupl_idx_pr[irqq]);
				TLorentzVector *mumi_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mumi_idx_pr[irqq]);

				bool HLTFilterPass_pr = ((Reco_QQ_trig_pr[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
				if (!HLTFilterPass_pr) continue;
				if (Reco_mu_whichGen_pr[Reco_QQ_mupl_idx_pr[irqq]] == -1) continue;
				if (Reco_mu_whichGen_pr[Reco_QQ_mumi_idx_pr[irqq]] == -1) continue;
				if (Reco_QQ_whichGen_pr[irqq] == -1) continue;
				if (Reco_QQ_sign_pr[irqq] != 0) continue;
				//if (Reco_QQ_VtxProb_pr[irqq] < 0.005) continue;
				if (JP_pr->M() < massLow || JP_pr->M() > massHigh) continue;
				if (!(IsAcceptanceQQ(mupl_pr->Pt(), mupl_pr->Eta()) && IsAcceptanceQQ(mumi_pr->Pt(), mumi_pr->Eta()))) continue;
				if (Centrality_pr < cLow || Centrality_pr > cHigh) continue;

				bool passMuonTypePl_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,1))) && 
					(Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,3)));
				bool passMuonTypeMi_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,1))) && 
					(Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,3)));
				if (!passMuonTypePl_pr || !passMuonTypeMi_pr) continue;

				bool muplSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 5) &&
					(Reco_mu_nPixWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 0) &&
					(fabs(Reco_mu_dxy_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 0.3) &&
					(fabs(Reco_mu_dz_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 20.) && passMuonTypePl_pr;
				bool mumiSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 5) &&
					(Reco_mu_nPixWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 0) &&
					(fabs(Reco_mu_dxy_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 0.3) &&
					(fabs(Reco_mu_dz_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 20.) && passMuonTypeMi_pr;
				if (!(muplSoft_pr && mumiSoft_pr)) continue;

				if (!((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
							(Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

				if (isTnP) {
					bool mupl_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mupl_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mumi_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mumi_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mupl_isL2_pr = (mupl_L2Filter_pr && !mupl_L3Filter_pr);
					bool mupl_isL3_pr = (mupl_L2Filter_pr && mupl_L3Filter_pr);
					bool mumi_isL2_pr = (mumi_L2Filter_pr && !mumi_L3Filter_pr);
					bool mumi_isL3_pr = (mumi_L2Filter_pr && mumi_L3Filter_pr);
					bool SelDone_pr = ((mupl_isL2_pr && mumi_isL3_pr) || (mupl_isL3_pr && mumi_isL2_pr) || (mupl_isL3_pr && mumi_isL3_pr));
					if (!SelDone_pr) continue;
				}

				if (HLTPass_pr && HLTFilterPass_pr) {
					double rap = fabs(JP_pr->Rapidity());
					double pt = JP_pr->Pt();
					double ctau = Reco_QQ_ctau3D_pr[irqq];

					// Calculate weight (Ncoll * Gen_weight * TnP * pT weight)
					double weight_pr_pt = findNcoll(Centrality_pr) * Gen_weight_pr;

					// pT weight
					double pt_weight_pr_pt = 1.0;
					if (rap < 1.6 && pt >= 6.5)
						pt_weight_pr_pt = fptw1sys->Eval(pt);
					if (rap >= 1.6 && rap < 2.4)
						pt_weight_pr_pt = fptw2sys->Eval(pt);

					// TnP weight
					double tnp_weight_pr_pt = 1.0;
					if (isTnP) {
						double pt1_pr = mupl_pr->Pt();
						double pt2_pr = mumi_pr->Pt();
						double eta1_pr = mupl_pr->Eta();
						double eta2_pr = mumi_pr->Eta();
						tnp_weight_pr_pt = tnp_weight_muid_pbpb(pt1_pr, eta1_pr, 0) * tnp_weight_muid_pbpb(pt2_pr, eta2_pr, 0);
						tnp_weight_pr_pt *= tnp_weight_trk_pbpb(eta1_pr, 0) * tnp_weight_trk_pbpb(eta2_pr, 0);

						// Trigger TnP weight
						bool mupl_L2Filter_pr_t = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3Filter_pr_t = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2Filter_pr_t = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3Filter_pr_t = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_pr_t = (mupl_L2Filter_pr_t && !mupl_L3Filter_pr_t);
						bool mupl_isL3_pr_t = (mupl_L2Filter_pr_t && mupl_L3Filter_pr_t);
						bool mumi_isL2_pr_t = (mumi_L2Filter_pr_t && !mumi_L3Filter_pr_t);
						bool mumi_isL3_pr_t = (mumi_L2Filter_pr_t && mumi_L3Filter_pr_t);

						double tnp_trig_mupl_pr = 1.0, tnp_trig_mumi_pr = 1.0;
						if (mupl_isL2_pr_t && mumi_isL3_pr_t) {
							tnp_trig_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
							tnp_trig_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
						} else if (mupl_isL3_pr_t && mumi_isL2_pr_t) {
							tnp_trig_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
							tnp_trig_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
						} else if (mupl_isL3_pr_t && mumi_isL3_pr_t) {
							int t_pr[2] = {-1, 1}; int l_pr = rand() % 2;
							if (t_pr[l_pr] == -1) {
								tnp_trig_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
								tnp_trig_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
							} else {
								tnp_trig_mupl_pr = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
								tnp_trig_mumi_pr = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
							}
						}
						tnp_weight_pr_pt *= tnp_trig_mupl_pr * tnp_trig_mumi_pr;
					}

					double total_weight_pr_pt = weight_pr_pt * pt_weight_pr_pt * tnp_weight_pr_pt;

					// Fill forward rapidity histograms
					if (rap >= 1.6 && rap < 2.4) {
						for (int ipt = 0; ipt < nPtBins_for; ipt++) {
							if (pt >= ptBin_for[ipt] && pt < ptBin_for[ipt+1]) {
								h_decay_pr_pt_for[ipt]->Fill(ctau, total_weight_pr_pt);
								break;
							}
						}
					}
					// Fill mid rapidity histograms
					if (rap < 1.6) {
						for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
							if (pt >= ptBin_mid[ipt] && pt < ptBin_mid[ipt+1]) {
								h_decay_pr_pt_mid[ipt]->Fill(ctau, total_weight_pr_pt);
								break;
							}
						}
					}
				}
			}
		}

		// Fill NPMC histograms for ALL pt bins in ONE loop
		cout << "Filling NPMC histograms for all pt bins..." << endl;
		for (int iev = 0; iev < treeNPMC->GetEntries(); ++iev) {
			if (iev % 10000000 == 0) cout << "NPMC pt-binned: " << iev << " / " << treeNPMC->GetEntries() << endl;
			treeNPMC->GetEntry(iev);
			bool HLTPass_np = ((HLTriggers_np & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

			for (Int_t irqq = 0; irqq < Reco_QQ_size_np; irqq++) {
				TLorentzVector *JP_np = (TLorentzVector*)Reco_QQ_4mom_np->At(irqq);
				TLorentzVector *mupl_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mupl_idx_np[irqq]);
				TLorentzVector *mumi_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mumi_idx_np[irqq]);

				bool HLTFilterPass_np = ((Reco_QQ_trig_np[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
				if (!HLTFilterPass_np) continue;
				if (Reco_mu_whichGen_np[Reco_QQ_mupl_idx_np[irqq]] == -1) continue;
				if (Reco_mu_whichGen_np[Reco_QQ_mumi_idx_np[irqq]] == -1) continue;
				if (Reco_QQ_whichGen_np[irqq] == -1) continue;
				if (Reco_QQ_sign_np[irqq] != 0) continue;
				//if (Reco_QQ_VtxProb_np[irqq] < 0.005) continue;
				if (JP_np->M() < massLow || JP_np->M() > massHigh) continue;
				if (!(IsAcceptanceQQ(mupl_np->Pt(), mupl_np->Eta()) && IsAcceptanceQQ(mumi_np->Pt(), mumi_np->Eta()))) continue;
				if (Centrality_np < cLow || Centrality_np > cHigh) continue;

				bool passMuonTypePl_np = (Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,1))) && 
					(Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,3)));
				bool passMuonTypeMi_np = (Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,1))) && 
					(Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,3)));
				if (!passMuonTypePl_np || !passMuonTypeMi_np) continue;

				bool muplSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mupl_idx_np[irqq]] > 5) &&
					(Reco_mu_nPixWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 0) &&
					(fabs(Reco_mu_dxy_np[Reco_QQ_mupl_idx_np[irqq]]) < 0.3) &&
					(fabs(Reco_mu_dz_np[Reco_QQ_mupl_idx_np[irqq]]) < 20.) && passMuonTypePl_np;
				bool mumiSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 5) &&
					(Reco_mu_nPixWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 0) &&
					(fabs(Reco_mu_dxy_np[Reco_QQ_mumi_idx_np[irqq]]) < 0.3) &&
					(fabs(Reco_mu_dz_np[Reco_QQ_mumi_idx_np[irqq]]) < 20.) && passMuonTypeMi_np;
				if (!(muplSoft_np && mumiSoft_np)) continue;

				if (!((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
							(Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

				if (isTnP) {
					bool mupl_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mupl_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mumi_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
					bool mumi_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
					bool mupl_isL2_np = (mupl_L2Filter_np && !mupl_L3Filter_np);
					bool mupl_isL3_np = (mupl_L2Filter_np && mupl_L3Filter_np);
					bool mumi_isL2_np = (mumi_L2Filter_np && !mumi_L3Filter_np);
					bool mumi_isL3_np = (mumi_L2Filter_np && mumi_L3Filter_np);
					bool SelDone_np = ((mupl_isL2_np && mumi_isL3_np) || (mupl_isL3_np && mumi_isL2_np) || (mupl_isL3_np && mumi_isL3_np));
					if (!SelDone_np) continue;
				}

				if (HLTPass_np && HLTFilterPass_np) {
					double rap = fabs(JP_np->Rapidity());
					double pt = JP_np->Pt();
					double ctau = Reco_QQ_ctau3D_np[irqq];

					// Calculate weight (Ncoll * Gen_weight * TnP * pT weight)
					double weight_np_pt = findNcoll(Centrality_np) * Gen_weight_np;

					// pT weight
					double pt_weight_np_pt = 1.0;
					if (rap < 1.6 && pt >= 6.5)
						pt_weight_np_pt = fptw1sys->Eval(pt);
					if (rap >= 1.6 && rap < 2.4)
						pt_weight_np_pt = fptw2sys->Eval(pt);

					// TnP weight
					double tnp_weight_np_pt = 1.0;
					if (isTnP) {
						double pt1_np = mupl_np->Pt();
						double pt2_np = mumi_np->Pt();
						double eta1_np = mupl_np->Eta();
						double eta2_np = mumi_np->Eta();
						tnp_weight_np_pt = tnp_weight_muid_pbpb(pt1_np, eta1_np, 0) * tnp_weight_muid_pbpb(pt2_np, eta2_np, 0);
						tnp_weight_np_pt *= tnp_weight_trk_pbpb(eta1_np, 0) * tnp_weight_trk_pbpb(eta2_np, 0);

						// Trigger TnP weight
						bool mupl_L2Filter_np_t = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3Filter_np_t = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2Filter_np_t = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3Filter_np_t = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_np_t = (mupl_L2Filter_np_t && !mupl_L3Filter_np_t);
						bool mupl_isL3_np_t = (mupl_L2Filter_np_t && mupl_L3Filter_np_t);
						bool mumi_isL2_np_t = (mumi_L2Filter_np_t && !mumi_L3Filter_np_t);
						bool mumi_isL3_np_t = (mumi_L2Filter_np_t && mumi_L3Filter_np_t);

						double tnp_trig_mupl_np = 1.0, tnp_trig_mumi_np = 1.0;
						if (mupl_isL2_np_t && mumi_isL3_np_t) {
							tnp_trig_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
							tnp_trig_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
						} else if (mupl_isL3_np_t && mumi_isL2_np_t) {
							tnp_trig_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
							tnp_trig_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
						} else if (mupl_isL3_np_t && mumi_isL3_np_t) {
							int t_np[2] = {-1, 1}; int l_np = rand() % 2;
							if (t_np[l_np] == -1) {
								tnp_trig_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
								tnp_trig_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
							} else {
								tnp_trig_mupl_np = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
								tnp_trig_mumi_np = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
							}
						}
						tnp_weight_np_pt *= tnp_trig_mupl_np * tnp_trig_mumi_np;
					}

					double total_weight_np_pt = weight_np_pt * pt_weight_np_pt * tnp_weight_np_pt;

					// Fill forward rapidity histograms
					if (rap >= 1.6 && rap < 2.4) {
						for (int ipt = 0; ipt < nPtBins_for; ipt++) {
							if (pt >= ptBin_for[ipt] && pt < ptBin_for[ipt+1]) {
								h_decay_np_pt_for[ipt]->Fill(ctau, total_weight_np_pt);
								break;
							}
						}
					}
					// Fill mid rapidity histograms
					if (rap < 1.6) {
						for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
							if (pt >= ptBin_mid[ipt] && pt < ptBin_mid[ipt+1]) {
								h_decay_np_pt_mid[ipt]->Fill(ctau, total_weight_np_pt);
								break;
							}
						}
					}
				}
			}
		}

		if(type == 1 || type == 3)
		{
			// Calculate ctau3D cut for each pt bin (forward rapidity)
			cout << "\n--- Forward rapidity (1.6-2.4) pt bins ---" << endl;
			for (int ipt = 0; ipt < nPtBins_for; ipt++) {
				// Use Integral() for weighted histograms
				double totalPR_pt = h_decay_pr_pt_for[ipt]->Integral();
				double totalNP_pt = h_decay_np_pt_for[ipt]->Integral();

				if (totalPR_pt > 0 && totalNP_pt > 0) {
					for (int bins = 0; bins < h_decay_pr_pt_for[ipt]->GetNbinsX(); bins++) {
						if (state == 1) {
							if (h_decay_pr_pt_for[ipt]->Integral(0, bins) < (0.899 * totalPR_pt) || h_decay_pr_pt_for[ipt]->Integral(0, bins) > (totalPR_pt * 0.9099)) continue;
							ctau3D_cut_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->GetBinCenter(bins);
							pr_eff_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR_pt;
							np_res_ptBin_for[ipt] = h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP_pt;
						} else {
							if (h_decay_pr_pt_for[ipt]->Integral(0, bins) < (0.97 * totalPR_pt) || h_decay_pr_pt_for[ipt]->Integral(0, bins) > (totalPR_pt * 0.98)) continue;
							ctau3D_cut_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->GetBinCenter(bins);
							pr_eff_ptBin_for[ipt] = 1 - h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR_pt;
							np_res_ptBin_for[ipt] = 1 - h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP_pt;
						}
					}
					cout << Form("pt %.1f-%.1f: ctau3D cut = %.4f, eff = %.3f, res = %.3f, integral PR=%.1f NP=%.1f", 
							ptBin_for[ipt], ptBin_for[ipt+1], ctau3D_cut_ptBin_for[ipt], pr_eff_ptBin_for[ipt], np_res_ptBin_for[ipt], totalPR_pt, totalNP_pt) << endl;

					// Draw efficiency plot for this pt bin
					TH1D* h_eff_pr_pt_for = (TH1D*)h_decay_pr_pt_for[ipt]->Clone(Form("h_eff_pr_pt_for_%d", ipt));
					TH1D* h_eff_np_pt_for = (TH1D*)h_decay_np_pt_for[ipt]->Clone(Form("h_eff_np_pt_for_%d", ipt));

					// Calculate CDF for NPMC and 1-CDF for PRMC
					for (int ib = 1; ib <= h_eff_pr_pt_for->GetNbinsX(); ib++) {
						h_eff_pr_pt_for->SetBinContent(ib, 1.0 - h_decay_pr_pt_for[ipt]->Integral(0, ib) / totalPR_pt);  // 1-CDF
						h_eff_np_pt_for->SetBinContent(ib, h_decay_np_pt_for[ipt]->Integral(0, ib) / totalNP_pt);  // CDF
					}

					TCanvas* c_eff_pt_for = new TCanvas(Form("c_eff_pt_for_%d", ipt), "", 800, 600);
					c_eff_pt_for->SetGrid();

					h_eff_np_pt_for->SetLineColor(kRed);
					h_eff_np_pt_for->SetLineWidth(2);
					h_eff_np_pt_for->GetXaxis()->SetTitle("ctau3D (mm)");
					h_eff_np_pt_for->GetYaxis()->SetTitle("Efficiency");
					h_eff_np_pt_for->GetYaxis()->SetRangeUser(0, 1.1);
					h_eff_np_pt_for->GetXaxis()->SetRangeUser(-1, 3);
					h_eff_np_pt_for->Draw("hist");

					h_eff_pr_pt_for->SetLineColor(kBlue);
					h_eff_pr_pt_for->SetLineWidth(2);
					h_eff_pr_pt_for->Draw("hist same");

					// Draw vertical line at ctau3D cut
					TLine* line_pt_for = new TLine(ctau3D_cut_ptBin_for[ipt], 0, ctau3D_cut_ptBin_for[ipt], 1);
					line_pt_for->SetLineStyle(2);
					line_pt_for->SetLineWidth(2);
					line_pt_for->Draw();

					TLegend* leg_pt_for = new TLegend(0.15, 0.75, 0.45, 0.88);
					leg_pt_for->SetBorderSize(0);
					leg_pt_for->SetFillStyle(0);
					leg_pt_for->AddEntry(h_eff_np_pt_for, "Non-Prompt J/#psi (CDF)", "l");
					leg_pt_for->AddEntry(h_eff_pr_pt_for, "Prompt J/#psi (1-CDF)", "l");
					leg_pt_for->Draw();

					drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptBin_for[ipt], ptBin_for[ipt+1]), 0.50, 0.85, 1, 18);
					drawText("1.6 < |y^{#mu#mu}| < 2.4", 0.50, 0.80, 1, 16);
					drawText(Form("Centrality %.0f-%.0f%%", cLow/2, cHigh/2), 0.50, 0.75, 1, 16);
					drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_ptBin_for[ipt]), 0.50, 0.65, 1, 16);
					if (state == 1) {
						drawText(Form("NP J/#psi eff: %.3f", np_res_ptBin_for[ipt]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_ptBin_for[ipt]), 0.50, 0.55, 2, 16);
					} else {
						drawText(Form("NP J/#psi eff: %.3f", 1-np_res_ptBin_for[ipt]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", pr_eff_ptBin_for[ipt]), 0.50, 0.55, 2, 16);
					}

					TString stateStr = (state == 1) ? "PRMC" : "NPMC";
					c_eff_pt_for->SaveAs(Form("./figs_ctau3D_v2/eff_ptBin_for_%s_pt%.1f-%.1f.pdf", stateStr.Data(), ptBin_for[ipt], ptBin_for[ipt+1]));
					c_eff_pt_for->SaveAs(Form("./figs_ctau3D_v2/eff_ptBin_for_%s_pt%.1f-%.1f.png", stateStr.Data(), ptBin_for[ipt], ptBin_for[ipt+1]));

					delete h_eff_pr_pt_for;
					delete h_eff_np_pt_for;
					delete c_eff_pt_for;
					delete line_pt_for;
					delete leg_pt_for;
				} else {
					cout << Form("pt %.1f-%.1f: Not enough entries (PR=%f, NP=%f)", ptBin_for[ipt], ptBin_for[ipt+1], totalPR_pt, totalNP_pt) << endl;
				}
				delete h_decay_pr_pt_for[ipt];
				delete h_decay_np_pt_for[ipt];
			}

			// Calculate ctau3D cut for each pt bin (mid rapidity)
			cout << "\n--- Mid rapidity (0-1.6) pt bins ---" << endl;
			for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
				// Use Integral() for weighted histograms
				double totalPR_pt = h_decay_pr_pt_mid[ipt]->Integral();
				double totalNP_pt = h_decay_np_pt_mid[ipt]->Integral();

				if (totalPR_pt > 0 && totalNP_pt > 0) {
					for (int bins = 0; bins < h_decay_pr_pt_mid[ipt]->GetNbinsX(); bins++) {
						if (state == 1) {
							if (h_decay_pr_pt_mid[ipt]->Integral(0, bins) < (0.899 * totalPR_pt) || h_decay_pr_pt_mid[ipt]->Integral(0, bins) > (totalPR_pt * 0.9099)) continue;
							ctau3D_cut_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->GetBinCenter(bins);
							pr_eff_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR_pt;
							np_res_ptBin_mid[ipt] = h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP_pt;
						} else {
							if (h_decay_pr_pt_mid[ipt]->Integral(0, bins) < (0.97 * totalPR_pt) || h_decay_pr_pt_mid[ipt]->Integral(0, bins) > (totalPR_pt * 0.98)) continue;
							ctau3D_cut_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->GetBinCenter(bins);
							pr_eff_ptBin_mid[ipt] = 1 - h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR_pt;
							np_res_ptBin_mid[ipt] = 1 - h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP_pt;
						}
					}
					cout << Form("pt %.1f-%.1f: ctau3D cut = %.4f, eff = %.3f, res = %.3f, integral PR=%.1f NP=%.1f", 
							ptBin_mid[ipt], ptBin_mid[ipt+1], ctau3D_cut_ptBin_mid[ipt], pr_eff_ptBin_mid[ipt], np_res_ptBin_mid[ipt], totalPR_pt, totalNP_pt) << endl;

					// Draw efficiency plot for this pt bin (mid)
					TH1D* h_eff_pr_pt_mid = (TH1D*)h_decay_pr_pt_mid[ipt]->Clone(Form("h_eff_pr_pt_mid_%d", ipt));
					TH1D* h_eff_np_pt_mid = (TH1D*)h_decay_np_pt_mid[ipt]->Clone(Form("h_eff_np_pt_mid_%d", ipt));

					for (int ib = 1; ib <= h_eff_pr_pt_mid->GetNbinsX(); ib++) {
						h_eff_pr_pt_mid->SetBinContent(ib, 1.0 - h_decay_pr_pt_mid[ipt]->Integral(0, ib) / totalPR_pt);
						h_eff_np_pt_mid->SetBinContent(ib, h_decay_np_pt_mid[ipt]->Integral(0, ib) / totalNP_pt);
					}

					TCanvas* c_eff_pt_mid = new TCanvas(Form("c_eff_pt_mid_%d", ipt), "", 800, 600);
					c_eff_pt_mid->SetGrid();

					h_eff_np_pt_mid->SetLineColor(kRed);
					h_eff_np_pt_mid->SetLineWidth(2);
					h_eff_np_pt_mid->GetXaxis()->SetTitle("ctau3D (mm)");
					h_eff_np_pt_mid->GetYaxis()->SetTitle("Efficiency");
					h_eff_np_pt_mid->GetYaxis()->SetRangeUser(0, 1.1);
					h_eff_np_pt_mid->GetXaxis()->SetRangeUser(-1, 3);
					h_eff_np_pt_mid->Draw("hist");

					h_eff_pr_pt_mid->SetLineColor(kBlue);
					h_eff_pr_pt_mid->SetLineWidth(2);
					h_eff_pr_pt_mid->Draw("hist same");

					TLine* line_pt_mid = new TLine(ctau3D_cut_ptBin_mid[ipt], 0, ctau3D_cut_ptBin_mid[ipt], 1);
					line_pt_mid->SetLineStyle(2);
					line_pt_mid->SetLineWidth(2);
					line_pt_mid->Draw();

					TLegend* leg_pt_mid = new TLegend(0.15, 0.75, 0.45, 0.88);
					leg_pt_mid->SetBorderSize(0);
					leg_pt_mid->SetFillStyle(0);
					leg_pt_mid->AddEntry(h_eff_np_pt_mid, "Non-Prompt J/#psi (CDF)", "l");
					leg_pt_mid->AddEntry(h_eff_pr_pt_mid, "Prompt J/#psi (1-CDF)", "l");
					leg_pt_mid->Draw();

					drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptBin_mid[ipt], ptBin_mid[ipt+1]), 0.50, 0.85, 1, 18);
					drawText("|y^{#mu#mu}| < 1.6", 0.50, 0.80, 1, 16);
					drawText(Form("Centrality %.0f-%.0f%%", cLow/2, cHigh/2), 0.50, 0.75, 1, 16);
					drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_ptBin_mid[ipt]), 0.50, 0.65, 1, 16);
					if (state == 1) {
						drawText(Form("NP J/#psi eff: %.3f", np_res_ptBin_mid[ipt]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_ptBin_mid[ipt]), 0.50, 0.55, 2, 16);
					} else {
						drawText(Form("NP J/#psi eff: %.3f", 1-np_res_ptBin_mid[ipt]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", pr_eff_ptBin_mid[ipt]), 0.50, 0.55, 2, 16);
					}

					TString stateStr = (state == 1) ? "PRMC" : "NPMC";
					c_eff_pt_mid->SaveAs(Form("./figs_ctau3D_v2/eff_ptBin_mid_%s_pt%.1f-%.1f.pdf", stateStr.Data(), ptBin_mid[ipt], ptBin_mid[ipt+1]));
					c_eff_pt_mid->SaveAs(Form("./figs_ctau3D_v2/eff_ptBin_mid_%s_pt%.1f-%.1f.png", stateStr.Data(), ptBin_mid[ipt], ptBin_mid[ipt+1]));

					delete h_eff_pr_pt_mid;
					delete h_eff_np_pt_mid;
					delete c_eff_pt_mid;
					delete line_pt_mid;
					delete leg_pt_mid;
				} else {
					cout << Form("pt %.1f-%.1f: Not enough entries (PR=%f, NP=%f)", ptBin_mid[ipt], ptBin_mid[ipt+1], totalPR_pt, totalNP_pt) << endl;
				}
				delete h_decay_pr_pt_mid[ipt];
				delete h_decay_np_pt_mid[ipt];
			}
		}

		if(type == 2 || type == 3) {
			//============================================================================
			// Centrality bin-wise ctau3D cut calculation (single pass)
			//============================================================================
			cout << "\n=== Calculating ctau3D cut for each centrality bin ===" << endl;


			// Create histograms for all centrality bins
			TH1D* h_decay_pr_cent_for[nCentBins_for];
			TH1D* h_decay_np_cent_for[nCentBins_for];
			TH1D* h_decay_pr_cent_mid[nCentBins_mid];
			TH1D* h_decay_np_cent_mid[nCentBins_mid];

			// Integrated histograms (forward: pt 3.5-40, cent 0-180; mid: pt 6.5-40, cent 0-180)
			TH1D* h_decay_pr_int_for = new TH1D("h_decay_pr_int_for", "", nbins, xmin, xmax);
			TH1D* h_decay_np_int_for = new TH1D("h_decay_np_int_for", "", nbins, xmin, xmax);
			TH1D* h_decay_pr_int_mid = new TH1D("h_decay_pr_int_mid", "", nbins, xmin, xmax);
			TH1D* h_decay_np_int_mid = new TH1D("h_decay_np_int_mid", "", nbins, xmin, xmax);

			for (int icent = 0; icent < nCentBins_for; icent++) {
				h_decay_pr_cent_for[icent] = new TH1D(Form("h_decay_pr_cent_for_%d", icent), "", nbins, xmin, xmax);
				h_decay_np_cent_for[icent] = new TH1D(Form("h_decay_np_cent_for_%d", icent), "", nbins, xmin, xmax);
			}
			for (int icent = 0; icent < nCentBins_mid; icent++) {
				h_decay_pr_cent_mid[icent] = new TH1D(Form("h_decay_pr_cent_mid_%d", icent), "", nbins, xmin, xmax);
				h_decay_np_cent_mid[icent] = new TH1D(Form("h_decay_np_cent_mid_%d", icent), "", nbins, xmin, xmax);
			}

			// Single pass over PRMC for all centrality bins
			cout << "Filling PRMC histograms for centrality bins..." << endl;
			for (int iev = 0; iev < treePRMC->GetEntries(); ++iev) {
				treePRMC->GetEntry(iev);
				bool HLTPass_pr = ((HLTriggers_pr & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

				for (Int_t irqq = 0; irqq < Reco_QQ_size_pr; irqq++) {
					TLorentzVector *JP_pr = (TLorentzVector*)Reco_QQ_4mom_pr->At(irqq);
					TLorentzVector *mupl_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mupl_idx_pr[irqq]);
					TLorentzVector *mumi_pr = (TLorentzVector*)Reco_mu_4mom_pr->At(Reco_QQ_mumi_idx_pr[irqq]);

					bool HLTFilterPass_pr = ((Reco_QQ_trig_pr[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
					if (!HLTFilterPass_pr) continue;
					if (Reco_mu_whichGen_pr[Reco_QQ_mupl_idx_pr[irqq]] == -1) continue;
					if (Reco_mu_whichGen_pr[Reco_QQ_mumi_idx_pr[irqq]] == -1) continue;
					if (Reco_QQ_whichGen_pr[irqq] == -1) continue;
					if (Reco_QQ_sign_pr[irqq] != 0) continue;
					//if (Reco_QQ_VtxProb_pr[irqq] < 0.005) continue;
					if (JP_pr->M() < massLow || JP_pr->M() > massHigh) continue;
					if (!(IsAcceptanceQQ(mupl_pr->Pt(), mupl_pr->Eta()) && IsAcceptanceQQ(mumi_pr->Pt(), mumi_pr->Eta()))) continue;

					bool passMuonTypePl_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,1))) && 
						(Reco_mu_SelectionType_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((int)pow(2,3)));
					bool passMuonTypeMi_pr = (Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,1))) && 
						(Reco_mu_SelectionType_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((int)pow(2,3)));
					if (!passMuonTypePl_pr || !passMuonTypeMi_pr) continue;

					bool muplSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 5) &&
						(Reco_mu_nPixWMea_pr[Reco_QQ_mupl_idx_pr[irqq]] > 0) &&
						(fabs(Reco_mu_dxy_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 0.3) &&
						(fabs(Reco_mu_dz_pr[Reco_QQ_mupl_idx_pr[irqq]]) < 20.) && passMuonTypePl_pr;
					bool mumiSoft_pr = (Reco_mu_nTrkWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 5) &&
						(Reco_mu_nPixWMea_pr[Reco_QQ_mumi_idx_pr[irqq]] > 0) &&
						(fabs(Reco_mu_dxy_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 0.3) &&
						(fabs(Reco_mu_dz_pr[Reco_QQ_mumi_idx_pr[irqq]]) < 20.) && passMuonTypeMi_pr;
					if (!(muplSoft_pr && mumiSoft_pr)) continue;

					if (!((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
								(Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

					if (isTnP) {
						bool mupl_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3Filter_pr = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_pr = (mupl_L2Filter_pr && !mupl_L3Filter_pr);
						bool mupl_isL3_pr = (mupl_L2Filter_pr && mupl_L3Filter_pr);
						bool mumi_isL2_pr = (mumi_L2Filter_pr && !mumi_L3Filter_pr);
						bool mumi_isL3_pr = (mumi_L2Filter_pr && mumi_L3Filter_pr);
						bool SelDone_pr = ((mupl_isL2_pr && mumi_isL3_pr) || (mupl_isL3_pr && mumi_isL2_pr) || (mupl_isL3_pr && mumi_isL3_pr));
						if (!SelDone_pr) continue;
					}

					if (!HLTPass_pr || !HLTFilterPass_pr) continue;

					double rapidity = fabs(JP_pr->Rapidity());
					double pt_cent = JP_pr->Pt();
					double ctau3D_val = Reco_QQ_ctau3D_pr[irqq];

					// Calculate weight for centrality bins
					double weight_pr_cent = findNcoll(Centrality_pr) * Gen_weight_pr;
					double pt_weight_pr_cent = 1.0;
					if (rapidity < 1.6 && pt_cent >= 6.5)
						pt_weight_pr_cent = fptw1sys->Eval(pt_cent);
					if (rapidity >= 1.6 && rapidity < 2.4)
						pt_weight_pr_cent = fptw2sys->Eval(pt_cent);

					double tnp_weight_pr_cent = 1.0;
					if (isTnP) {
						double pt1_pr = mupl_pr->Pt(); double pt2_pr = mumi_pr->Pt();
						double eta1_pr = mupl_pr->Eta(); double eta2_pr = mumi_pr->Eta();
						tnp_weight_pr_cent = tnp_weight_muid_pbpb(pt1_pr, eta1_pr, 0) * tnp_weight_muid_pbpb(pt2_pr, eta2_pr, 0);
						tnp_weight_pr_cent *= tnp_weight_trk_pbpb(eta1_pr, 0) * tnp_weight_trk_pbpb(eta2_pr, 0);

						bool mupl_L2_c = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3_c = ((Reco_mu_trig_pr[Reco_QQ_mupl_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2_c = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3_c = ((Reco_mu_trig_pr[Reco_QQ_mumi_idx_pr[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_c = (mupl_L2_c && !mupl_L3_c); bool mupl_isL3_c = (mupl_L2_c && mupl_L3_c);
						bool mumi_isL2_c = (mumi_L2_c && !mumi_L3_c); bool mumi_isL3_c = (mumi_L2_c && mumi_L3_c);

						double tnp_trig_mupl_c = 1.0, tnp_trig_mumi_c = 1.0;
						if (mupl_isL2_c && mumi_isL3_c) {
							tnp_trig_mupl_c = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
							tnp_trig_mumi_c = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
						} else if (mupl_isL3_c && mumi_isL2_c) {
							tnp_trig_mupl_c = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
							tnp_trig_mumi_c = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
						} else if (mupl_isL3_c && mumi_isL3_c) {
							int t_c[2] = {-1, 1}; int l_c = rand() % 2;
							if (t_c[l_c] == -1) {
								tnp_trig_mupl_c = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 2, 0);
								tnp_trig_mumi_c = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 3, 0);
							} else {
								tnp_trig_mupl_c = tnp_weight_trg_pbpb(mupl_pr->Pt(), mupl_pr->Eta(), 3, 0);
								tnp_trig_mumi_c = tnp_weight_trg_pbpb(mumi_pr->Pt(), mumi_pr->Eta(), 2, 0);
							}
						}
						tnp_weight_pr_cent *= tnp_trig_mupl_c * tnp_trig_mumi_c;
					}
					double total_weight_pr_cent = weight_pr_cent * pt_weight_pr_cent * tnp_weight_pr_cent;

					// Forward rapidity (1.6-2.4): check pt and fill cent bins + integrated
					if (rapidity >= 1.6 && rapidity < 2.4 && JP_pr->Pt() >= 3.5 && JP_pr->Pt() < 40) {
						// Fill integrated histogram (cent 0-180)
						if (Centrality_pr >= 0 && Centrality_pr < 180) {
							h_decay_pr_int_for->Fill(ctau3D_val, total_weight_pr_cent);
						}
						// Fill cent bins
						for (int icent = 0; icent < nCentBins_for; icent++) {
							if (Centrality_pr >= centBin_for[icent] && Centrality_pr < centBin_for[icent+1]) {
								h_decay_pr_cent_for[icent]->Fill(ctau3D_val, total_weight_pr_cent);
								break;
							}
						}
					}
					// Mid rapidity (0-1.6): check pt and fill cent bins + integrated
					if (rapidity < 1.6 && JP_pr->Pt() >= 6.5 && JP_pr->Pt() < 40) {
						// Fill integrated histogram (cent 0-180)
						if (Centrality_pr >= 0 && Centrality_pr < 180) {
							h_decay_pr_int_mid->Fill(ctau3D_val, total_weight_pr_cent);
						}
						// Fill cent bins
						for (int icent = 0; icent < nCentBins_mid; icent++) {
							if (Centrality_pr >= centBin_mid[icent] && Centrality_pr < centBin_mid[icent+1]) {
								h_decay_pr_cent_mid[icent]->Fill(ctau3D_val, total_weight_pr_cent);
								break;
							}
						}
					}
				}
			}

			// Single pass over NPMC for all centrality bins
			cout << "Filling NPMC histograms for centrality bins..." << endl;
			for (int iev = 0; iev < treeNPMC->GetEntries(); ++iev) {
				treeNPMC->GetEntry(iev);
				bool HLTPass_np = ((HLTriggers_np & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));

				for (Int_t irqq = 0; irqq < Reco_QQ_size_np; irqq++) {
					TLorentzVector *JP_np = (TLorentzVector*)Reco_QQ_4mom_np->At(irqq);
					TLorentzVector *mupl_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mupl_idx_np[irqq]);
					TLorentzVector *mumi_np = (TLorentzVector*)Reco_mu_4mom_np->At(Reco_QQ_mumi_idx_np[irqq]);

					bool HLTFilterPass_np = ((Reco_QQ_trig_np[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)));
					if (!HLTFilterPass_np) continue;
					if (Reco_mu_whichGen_np[Reco_QQ_mupl_idx_np[irqq]] == -1) continue;
					if (Reco_mu_whichGen_np[Reco_QQ_mumi_idx_np[irqq]] == -1) continue;
					if (Reco_QQ_whichGen_np[irqq] == -1) continue;
					if (Reco_QQ_sign_np[irqq] != 0) continue;
					//if (Reco_QQ_VtxProb_np[irqq] < 0.005) continue;
					if (JP_np->M() < massLow || JP_np->M() > massHigh) continue;
					if (!(IsAcceptanceQQ(mupl_np->Pt(), mupl_np->Eta()) && IsAcceptanceQQ(mumi_np->Pt(), mumi_np->Eta()))) continue;

					bool passMuonTypePl_np = (Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,1))) && 
						(Reco_mu_SelectionType_np[Reco_QQ_mupl_idx_np[irqq]] & ((int)pow(2,3)));
					bool passMuonTypeMi_np = (Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,1))) && 
						(Reco_mu_SelectionType_np[Reco_QQ_mumi_idx_np[irqq]] & ((int)pow(2,3)));
					if (!passMuonTypePl_np || !passMuonTypeMi_np) continue;

					bool muplSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mupl_idx_np[irqq]] > 5) &&
						(Reco_mu_nPixWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 0) &&
						(fabs(Reco_mu_dxy_np[Reco_QQ_mupl_idx_np[irqq]]) < 0.3) &&
						(fabs(Reco_mu_dz_np[Reco_QQ_mupl_idx_np[irqq]]) < 20.) && passMuonTypePl_np;
					bool mumiSoft_np = (Reco_mu_nTrkWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 5) &&
						(Reco_mu_nPixWMea_np[Reco_QQ_mumi_idx_np[irqq]] > 0) &&
						(fabs(Reco_mu_dxy_np[Reco_QQ_mumi_idx_np[irqq]]) < 0.3) &&
						(fabs(Reco_mu_dz_np[Reco_QQ_mumi_idx_np[irqq]]) < 20.) && passMuonTypeMi_np;
					if (!(muplSoft_np && mumiSoft_np)) continue;

					if (!((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && 
								(Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)))) continue;

					if (isTnP) {
						bool mupl_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3Filter_np = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_np = (mupl_L2Filter_np && !mupl_L3Filter_np);
						bool mupl_isL3_np = (mupl_L2Filter_np && mupl_L3Filter_np);
						bool mumi_isL2_np = (mumi_L2Filter_np && !mumi_L3Filter_np);
						bool mumi_isL3_np = (mumi_L2Filter_np && mumi_L3Filter_np);
						bool SelDone_np = ((mupl_isL2_np && mumi_isL3_np) || (mupl_isL3_np && mumi_isL2_np) || (mupl_isL3_np && mumi_isL3_np));
						if (!SelDone_np) continue;
					}

					if (!HLTPass_np || !HLTFilterPass_np) continue;

					double rapidity = fabs(JP_np->Rapidity());
					double pt_cent_np = JP_np->Pt();
					double ctau3D_val = Reco_QQ_ctau3D_np[irqq];

					// Calculate weight for centrality bins
					double weight_np_cent = findNcoll(Centrality_np) * Gen_weight_np;
					double pt_weight_np_cent = 1.0;
					if (rapidity < 1.6 && pt_cent_np >= 6.5)
						pt_weight_np_cent = fptw1sys->Eval(pt_cent_np);
					if (rapidity >= 1.6 && rapidity < 2.4)
						pt_weight_np_cent = fptw2sys->Eval(pt_cent_np);

					double tnp_weight_np_cent = 1.0;
					if (isTnP) {
						double pt1_np = mupl_np->Pt(); double pt2_np = mumi_np->Pt();
						double eta1_np = mupl_np->Eta(); double eta2_np = mumi_np->Eta();
						tnp_weight_np_cent = tnp_weight_muid_pbpb(pt1_np, eta1_np, 0) * tnp_weight_muid_pbpb(pt2_np, eta2_np, 0);
						tnp_weight_np_cent *= tnp_weight_trk_pbpb(eta1_np, 0) * tnp_weight_trk_pbpb(eta2_np, 0);

						bool mupl_L2_cn = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mupl_L3_cn = ((Reco_mu_trig_np[Reco_QQ_mupl_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mumi_L2_cn = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)));
						bool mumi_L3_cn = ((Reco_mu_trig_np[Reco_QQ_mumi_idx_np[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)));
						bool mupl_isL2_cn = (mupl_L2_cn && !mupl_L3_cn); bool mupl_isL3_cn = (mupl_L2_cn && mupl_L3_cn);
						bool mumi_isL2_cn = (mumi_L2_cn && !mumi_L3_cn); bool mumi_isL3_cn = (mumi_L2_cn && mumi_L3_cn);

						double tnp_trig_mupl_cn = 1.0, tnp_trig_mumi_cn = 1.0;
						if (mupl_isL2_cn && mumi_isL3_cn) {
							tnp_trig_mupl_cn = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
							tnp_trig_mumi_cn = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
						} else if (mupl_isL3_cn && mumi_isL2_cn) {
							tnp_trig_mupl_cn = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
							tnp_trig_mumi_cn = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
						} else if (mupl_isL3_cn && mumi_isL3_cn) {
							int t_cn[2] = {-1, 1}; int l_cn = rand() % 2;
							if (t_cn[l_cn] == -1) {
								tnp_trig_mupl_cn = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 2, 0);
								tnp_trig_mumi_cn = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 3, 0);
							} else {
								tnp_trig_mupl_cn = tnp_weight_trg_pbpb(mupl_np->Pt(), mupl_np->Eta(), 3, 0);
								tnp_trig_mumi_cn = tnp_weight_trg_pbpb(mumi_np->Pt(), mumi_np->Eta(), 2, 0);
							}
						}
						tnp_weight_np_cent *= tnp_trig_mupl_cn * tnp_trig_mumi_cn;
					}
					double total_weight_np_cent = weight_np_cent * pt_weight_np_cent * tnp_weight_np_cent;

					// Forward rapidity (1.6-2.4): check pt and fill cent bins + integrated
					if (rapidity >= 1.6 && rapidity < 2.4 && JP_np->Pt() >= 3.5 && JP_np->Pt() < 40) {
						// Fill integrated histogram (cent 0-180)
						if (Centrality_np >= 0 && Centrality_np < 180) {
							h_decay_np_int_for->Fill(ctau3D_val, total_weight_np_cent);
						}
						// Fill cent bins
						for (int icent = 0; icent < nCentBins_for; icent++) {
							if (Centrality_np >= centBin_for[icent] && Centrality_np < centBin_for[icent+1]) {
								h_decay_np_cent_for[icent]->Fill(ctau3D_val, total_weight_np_cent);
								break;
							}
						}
					}
					// Mid rapidity (0-1.6): check pt and fill cent bins + integrated
					if (rapidity < 1.6 && JP_np->Pt() >= 6.5 && JP_np->Pt() < 40) {
						// Fill integrated histogram (cent 0-180)
						if (Centrality_np >= 0 && Centrality_np < 180) {
							h_decay_np_int_mid->Fill(ctau3D_val, total_weight_np_cent);
						}
						// Fill cent bins
						for (int icent = 0; icent < nCentBins_mid; icent++) {
							if (Centrality_np >= centBin_mid[icent] && Centrality_np < centBin_mid[icent+1]) {
								h_decay_np_cent_mid[icent]->Fill(ctau3D_val, total_weight_np_cent);
								break;
							}
						}
					}
				}
			}

			//============================================================================
			// Calculate integrated ctau3D cut (forward: pt 3.5-40, cent 0-180; mid: pt 6.5-40, cent 0-180)
			//============================================================================
			cout << "\n=== Integrated ctau3D cut calculation ===" << endl;

			// Forward rapidity integrated
			cout << "\n--- Forward rapidity (1.6-2.4) integrated (pt 3.5-40, cent 0-90%) ---" << endl;
			// Use Integral() for weighted histograms
			double totalPR_int_for = h_decay_pr_int_for->Integral();
			double totalNP_int_for = h_decay_np_int_for->Integral();

			if (totalPR_int_for > 0 && totalNP_int_for > 0) {
				for (int bins = 0; bins < h_decay_pr_int_for->GetNbinsX(); bins++) {
					if (state == 1) {
						if (h_decay_pr_int_for->Integral(0, bins) < (0.899 * totalPR_int_for) || 
								h_decay_pr_int_for->Integral(0, bins) > (totalPR_int_for * 0.9099)) continue;
						ctau3D_cut_int_for = h_decay_pr_int_for->GetBinCenter(bins);
						pr_eff_int_for = h_decay_pr_int_for->Integral(0, bins) / totalPR_int_for;
						np_res_int_for = h_decay_np_int_for->Integral(0, bins) / totalNP_int_for;
					} else {
						if (h_decay_pr_int_for->Integral(0, bins) < (0.97 * totalPR_int_for) || 
								h_decay_pr_int_for->Integral(0, bins) > (totalPR_int_for * 0.98)) continue;
						ctau3D_cut_int_for = h_decay_pr_int_for->GetBinCenter(bins);
						pr_eff_int_for = 1 - h_decay_pr_int_for->Integral(0, bins) / totalPR_int_for;
						np_res_int_for = 1 - h_decay_np_int_for->Integral(0, bins) / totalNP_int_for;
					}
				}
				cout << Form("Integrated: ctau3D cut = %.4f, PR eff = %.3f, NP res = %.3f, integral PR=%.1f NP=%.1f", 
						ctau3D_cut_int_for, pr_eff_int_for, np_res_int_for, totalPR_int_for, totalNP_int_for) << endl;

				// Draw efficiency plot for integrated forward
				TH1D* h_eff_pr_int_for = (TH1D*)h_decay_pr_int_for->Clone("h_eff_pr_int_for_plot");
				TH1D* h_eff_np_int_for = (TH1D*)h_decay_np_int_for->Clone("h_eff_np_int_for_plot");

				for (int ib = 1; ib <= h_eff_pr_int_for->GetNbinsX(); ib++) {
					h_eff_pr_int_for->SetBinContent(ib, 1.0 - h_decay_pr_int_for->Integral(0, ib) / totalPR_int_for);
					h_eff_np_int_for->SetBinContent(ib, h_decay_np_int_for->Integral(0, ib) / totalNP_int_for);
				}

				TCanvas* c_eff_int_for = new TCanvas("c_eff_int_for", "", 800, 600);
				c_eff_int_for->SetGrid();

				h_eff_np_int_for->SetLineColor(kRed);
				h_eff_np_int_for->SetLineWidth(2);
				h_eff_np_int_for->GetXaxis()->SetTitle("ctau3D (mm)");
				h_eff_np_int_for->GetYaxis()->SetTitle("Efficiency");
				h_eff_np_int_for->GetYaxis()->SetRangeUser(0, 1.1);
				h_eff_np_int_for->GetXaxis()->SetRangeUser(-1, 3);
				h_eff_np_int_for->Draw("hist");

				h_eff_pr_int_for->SetLineColor(kBlue);
				h_eff_pr_int_for->SetLineWidth(2);
				h_eff_pr_int_for->Draw("hist same");

				TLine* line_int_for = new TLine(ctau3D_cut_int_for, 0, ctau3D_cut_int_for, 1);
				line_int_for->SetLineStyle(2);
				line_int_for->SetLineWidth(2);
				line_int_for->Draw();

				TLegend* leg_int_for = new TLegend(0.15, 0.75, 0.45, 0.88);
				leg_int_for->SetBorderSize(0);
				leg_int_for->SetFillStyle(0);
				leg_int_for->AddEntry(h_eff_np_int_for, "Non-Prompt J/#psi (CDF)", "l");
				leg_int_for->AddEntry(h_eff_pr_int_for, "Prompt J/#psi (1-CDF)", "l");
				leg_int_for->Draw();

				drawText("3.5 < p_{T}^{#mu#mu} < 40.0 GeV/c", 0.50, 0.85, 1, 18);
				drawText("1.6 < |y^{#mu#mu}| < 2.4", 0.50, 0.80, 1, 16);
				drawText("Centrality 0-90%", 0.50, 0.75, 1, 16);
				drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_int_for), 0.50, 0.65, 1, 16);
				if (state == 1) {
					drawText(Form("NP J/#psi eff: %.3f", np_res_int_for), 0.50, 0.60, 4, 16);
					drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_int_for), 0.50, 0.55, 2, 16);
				} else {
					drawText(Form("NP J/#psi eff: %.3f", 1-np_res_int_for), 0.50, 0.60, 4, 16);
					drawText(Form("PR J/#psi res: %.3f", pr_eff_int_for), 0.50, 0.55, 2, 16);
				}

				TString stateStr = (state == 1) ? "PRMC" : "NPMC";
				c_eff_int_for->SaveAs(Form("./figs_ctau3D_v2/eff_integrated_for_%s_pt3.5-40_cent0-90.pdf", stateStr.Data()));
				c_eff_int_for->SaveAs(Form("./figs_ctau3D_v2/eff_integrated_for_%s_pt3.5-40_cent0-90.png", stateStr.Data()));

				delete h_eff_pr_int_for;
				delete h_eff_np_int_for;
				delete c_eff_int_for;
				delete line_int_for;
				delete leg_int_for;
			}

			// Mid rapidity integrated
			cout << "\n--- Mid rapidity (0-1.6) integrated (pt 6.5-40, cent 0-90%) ---" << endl;
			// Use Integral() for weighted histograms
			double totalPR_int_mid = h_decay_pr_int_mid->Integral();
			double totalNP_int_mid = h_decay_np_int_mid->Integral();

			if (totalPR_int_mid > 0 && totalNP_int_mid > 0) {
				for (int bins = 0; bins < h_decay_pr_int_mid->GetNbinsX(); bins++) {
					if (state == 1) {
						if (h_decay_pr_int_mid->Integral(0, bins) < (0.899 * totalPR_int_mid) || 
								h_decay_pr_int_mid->Integral(0, bins) > (totalPR_int_mid * 0.9099)) continue;
						ctau3D_cut_int_mid = h_decay_pr_int_mid->GetBinCenter(bins);
						pr_eff_int_mid = h_decay_pr_int_mid->Integral(0, bins) / totalPR_int_mid;
						np_res_int_mid = h_decay_np_int_mid->Integral(0, bins) / totalNP_int_mid;
					} else {
						if (h_decay_pr_int_mid->Integral(0, bins) < (0.97 * totalPR_int_mid) || 
								h_decay_pr_int_mid->Integral(0, bins) > (totalPR_int_mid * 0.98)) continue;
						ctau3D_cut_int_mid = h_decay_pr_int_mid->GetBinCenter(bins);
						pr_eff_int_mid = 1 - h_decay_pr_int_mid->Integral(0, bins) / totalPR_int_mid;
						np_res_int_mid = 1 - h_decay_np_int_mid->Integral(0, bins) / totalNP_int_mid;
					}
				}
				cout << Form("Integrated: ctau3D cut = %.4f, PR eff = %.3f, NP res = %.3f, integral PR=%.1f NP=%.1f", 
						ctau3D_cut_int_mid, pr_eff_int_mid, np_res_int_mid, totalPR_int_mid, totalNP_int_mid) << endl;

				// Draw efficiency plot for integrated mid
				TH1D* h_eff_pr_int_mid = (TH1D*)h_decay_pr_int_mid->Clone("h_eff_pr_int_mid_plot");
				TH1D* h_eff_np_int_mid = (TH1D*)h_decay_np_int_mid->Clone("h_eff_np_int_mid_plot");

				for (int ib = 1; ib <= h_eff_pr_int_mid->GetNbinsX(); ib++) {
					h_eff_pr_int_mid->SetBinContent(ib, 1.0 - h_decay_pr_int_mid->Integral(0, ib) / totalPR_int_mid);
					h_eff_np_int_mid->SetBinContent(ib, h_decay_np_int_mid->Integral(0, ib) / totalNP_int_mid);
				}

				TCanvas* c_eff_int_mid = new TCanvas("c_eff_int_mid", "", 800, 600);
				c_eff_int_mid->SetGrid();

				h_eff_np_int_mid->SetLineColor(kRed);
				h_eff_np_int_mid->SetLineWidth(2);
				h_eff_np_int_mid->GetXaxis()->SetTitle("ctau3D (mm)");
				h_eff_np_int_mid->GetYaxis()->SetTitle("Efficiency");
				h_eff_np_int_mid->GetYaxis()->SetRangeUser(0, 1.1);
				h_eff_np_int_mid->GetXaxis()->SetRangeUser(-1, 3);
				h_eff_np_int_mid->Draw("hist");

				h_eff_pr_int_mid->SetLineColor(kBlue);
				h_eff_pr_int_mid->SetLineWidth(2);
				h_eff_pr_int_mid->Draw("hist same");

				TLine* line_int_mid = new TLine(ctau3D_cut_int_mid, 0, ctau3D_cut_int_mid, 1);
				line_int_mid->SetLineStyle(2);
				line_int_mid->SetLineWidth(2);
				line_int_mid->Draw();

				TLegend* leg_int_mid = new TLegend(0.15, 0.75, 0.45, 0.88);
				leg_int_mid->SetBorderSize(0);
				leg_int_mid->SetFillStyle(0);
				leg_int_mid->AddEntry(h_eff_np_int_mid, "Non-Prompt J/#psi (CDF)", "l");
				leg_int_mid->AddEntry(h_eff_pr_int_mid, "Prompt J/#psi (1-CDF)", "l");
				leg_int_mid->Draw();

				drawText("6.5 < p_{T}^{#mu#mu} < 40.0 GeV/c", 0.50, 0.85, 1, 18);
				drawText("|y^{#mu#mu}| < 1.6", 0.50, 0.80, 1, 16);
				drawText("Centrality 0-90%", 0.50, 0.75, 1, 16);
				drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_int_mid), 0.50, 0.65, 1, 16);
				if (state == 1) {
					drawText(Form("NP J/#psi eff: %.3f", np_res_int_mid), 0.50, 0.60, 4, 16);
					drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_int_mid), 0.50, 0.55, 2, 16);
				} else {
					drawText(Form("NP J/#psi eff: %.3f", 1-np_res_int_mid), 0.50, 0.60, 4, 16);
					drawText(Form("PR J/#psi res: %.3f", pr_eff_int_mid), 0.50, 0.55, 2, 16);
				}

				TString stateStr = (state == 1) ? "PRMC" : "NPMC";
				c_eff_int_mid->SaveAs(Form("./figs_ctau3D_v2/eff_integrated_mid_%s_pt6.5-40_cent0-90.pdf", stateStr.Data()));
				c_eff_int_mid->SaveAs(Form("./figs_ctau3D_v2/eff_integrated_mid_%s_pt6.5-40_cent0-90.png", stateStr.Data()));

				delete h_eff_pr_int_mid;
				delete h_eff_np_int_mid;
				delete c_eff_int_mid;
				delete line_int_mid;
				delete leg_int_mid;
			}

			// Save integrated values to ROOT file
			TString rootFileName_int;
			if (state == 1) {
				rootFileName_int = "./roots_ctau3D_v2/ctau3D_cut_integrated_PRMC.root";
			} else {
				rootFileName_int = "./roots_ctau3D_v2/ctau3D_cut_integrated_NPMC.root";
			}
			TFile *wf_int = new TFile(rootFileName_int, "RECREATE");
			wf_int->cd();
			h_decay_pr_int_for->Write();
			h_decay_np_int_for->Write();
			h_decay_pr_int_mid->Write();
			h_decay_np_int_mid->Write();

			// Save ctau3D cut values as TNamed
			TNamed* n_ctau_for = new TNamed("ctau3D_cut_for", Form("%.6f", ctau3D_cut_int_for));
			TNamed* n_ctau_mid = new TNamed("ctau3D_cut_mid", Form("%.6f", ctau3D_cut_int_mid));
			n_ctau_for->Write();
			n_ctau_mid->Write();

			wf_int->Close();
			cout << "Saved integrated ctau3D cuts to: " << rootFileName_int << endl;

			// Cleanup integrated histograms
			delete h_decay_pr_int_for;
			delete h_decay_np_int_for;
			delete h_decay_pr_int_mid;
			delete h_decay_np_int_mid;

			// Calculate ctau3D cut for each centrality bin (forward)
			cout << "\n--- Forward rapidity (1.6-2.4) centrality bins ---" << endl;
			for (int icent = 0; icent < nCentBins_for; icent++) {
				double centLow_bin = centBin_for[icent];
				double centHigh_bin = centBin_for[icent + 1];

				// Use Integral() for weighted histograms
				double totalPR = h_decay_pr_cent_for[icent]->Integral();
				double totalNP = h_decay_np_cent_for[icent]->Integral();

				if (totalPR > 0 && totalNP > 0) {
					for (int bins = 0; bins < h_decay_pr_cent_for[icent]->GetNbinsX(); bins++) {
						if (state == 1) {
							if (h_decay_pr_cent_for[icent]->Integral(0, bins) < (0.899 * totalPR) || 
									h_decay_pr_cent_for[icent]->Integral(0, bins) > (totalPR * 0.9099)) continue;
							ctau3D_cut_centBin_for[icent] = h_decay_pr_cent_for[icent]->GetBinCenter(bins);
							pr_eff_centBin_for[icent] = h_decay_pr_cent_for[icent]->Integral(0, bins) / totalPR;
							np_res_centBin_for[icent] = h_decay_np_cent_for[icent]->Integral(0, bins) / totalNP;
						} else {
							if (h_decay_pr_cent_for[icent]->Integral(0, bins) < (0.97 * totalPR) || 
									h_decay_pr_cent_for[icent]->Integral(0, bins) > (totalPR * 0.98)) continue;
							ctau3D_cut_centBin_for[icent] = h_decay_pr_cent_for[icent]->GetBinCenter(bins);
							pr_eff_centBin_for[icent] = 1 - h_decay_pr_cent_for[icent]->Integral(0, bins) / totalPR;
							np_res_centBin_for[icent] = 1 - h_decay_np_cent_for[icent]->Integral(0, bins) / totalNP;
						}
					}
					cout << Form("cent %.0f-%.0f: ctau3D cut = %.4f, PR eff = %.3f, NP res = %.3f, integral PR=%.1f NP=%.1f", 
							centLow_bin/2, centHigh_bin/2, ctau3D_cut_centBin_for[icent], 
							pr_eff_centBin_for[icent], np_res_centBin_for[icent], totalPR, totalNP) << endl;

					// Draw efficiency plot for this centrality bin (forward)
					TH1D* h_eff_pr_cent_for = (TH1D*)h_decay_pr_cent_for[icent]->Clone(Form("h_eff_pr_cent_for_%d", icent));
					TH1D* h_eff_np_cent_for = (TH1D*)h_decay_np_cent_for[icent]->Clone(Form("h_eff_np_cent_for_%d", icent));

					for (int ib = 1; ib <= h_eff_pr_cent_for->GetNbinsX(); ib++) {
						h_eff_pr_cent_for->SetBinContent(ib, 1.0 - h_decay_pr_cent_for[icent]->Integral(0, ib) / totalPR);
						h_eff_np_cent_for->SetBinContent(ib, h_decay_np_cent_for[icent]->Integral(0, ib) / totalNP);
					}

					TCanvas* c_eff_cent_for = new TCanvas(Form("c_eff_cent_for_%d", icent), "", 800, 600);
					c_eff_cent_for->SetGrid();

					h_eff_np_cent_for->SetLineColor(kRed);
					h_eff_np_cent_for->SetLineWidth(2);
					h_eff_np_cent_for->GetXaxis()->SetTitle("ctau3D (mm)");
					h_eff_np_cent_for->GetYaxis()->SetTitle("Efficiency");
					h_eff_np_cent_for->GetYaxis()->SetRangeUser(0, 1.1);
					h_eff_np_cent_for->GetXaxis()->SetRangeUser(-1, 3);
					h_eff_np_cent_for->Draw("hist");

					h_eff_pr_cent_for->SetLineColor(kBlue);
					h_eff_pr_cent_for->SetLineWidth(2);
					h_eff_pr_cent_for->Draw("hist same");

					TLine* line_cent_for = new TLine(ctau3D_cut_centBin_for[icent], 0, ctau3D_cut_centBin_for[icent], 1);
					line_cent_for->SetLineStyle(2);
					line_cent_for->SetLineWidth(2);
					line_cent_for->Draw();

					TLegend* leg_cent_for = new TLegend(0.15, 0.75, 0.45, 0.88);
					leg_cent_for->SetBorderSize(0);
					leg_cent_for->SetFillStyle(0);
					leg_cent_for->AddEntry(h_eff_np_cent_for, "Non-Prompt J/#psi (CDF)", "l");
					leg_cent_for->AddEntry(h_eff_pr_cent_for, "Prompt J/#psi (1-CDF)", "l");
					leg_cent_for->Draw();

					drawText("3.5 < p_{T}^{#mu#mu} < 40.0 GeV/c", 0.50, 0.85, 1, 18);
					drawText("1.6 < |y^{#mu#mu}| < 2.4", 0.50, 0.80, 1, 16);
					drawText(Form("Centrality %.0f-%.0f%%", centLow_bin/2, centHigh_bin/2), 0.50, 0.75, 1, 16);
					drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_centBin_for[icent]), 0.50, 0.65, 1, 16);
					if (state == 1) {
						drawText(Form("NP J/#psi eff: %.3f", np_res_centBin_for[icent]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_centBin_for[icent]), 0.50, 0.55, 2, 16);
					} else {
						drawText(Form("NP J/#psi eff: %.3f", 1-np_res_centBin_for[icent]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", pr_eff_centBin_for[icent]), 0.50, 0.55, 2, 16);
					}

					TString stateStr = (state == 1) ? "PRMC" : "NPMC";
					c_eff_cent_for->SaveAs(Form("./figs_ctau3D_v2/eff_centBin_for_%s_cent%.0f-%.0f.pdf", stateStr.Data(), centLow_bin/2, centHigh_bin/2));
					c_eff_cent_for->SaveAs(Form("./figs_ctau3D_v2/eff_centBin_for_%s_cent%.0f-%.0f.png", stateStr.Data(), centLow_bin/2, centHigh_bin/2));

					delete h_eff_pr_cent_for;
					delete h_eff_np_cent_for;
					delete c_eff_cent_for;
					delete line_cent_for;
					delete leg_cent_for;
				}
			}

			// Calculate ctau3D cut for each centrality bin (mid)
			cout << "\n--- Mid rapidity (0-1.6) centrality bins ---" << endl;
			for (int icent = 0; icent < nCentBins_mid; icent++) {
				double centLow_bin = centBin_mid[icent];
				double centHigh_bin = centBin_mid[icent + 1];

				// Use Integral() for weighted histograms
				double totalPR = h_decay_pr_cent_mid[icent]->Integral();
				double totalNP = h_decay_np_cent_mid[icent]->Integral();

				if (totalPR > 0 && totalNP > 0) {
					for (int bins = 0; bins < h_decay_pr_cent_mid[icent]->GetNbinsX(); bins++) {
						if (state == 1) {
							if (h_decay_pr_cent_mid[icent]->Integral(0, bins) < (0.899 * totalPR) || 
									h_decay_pr_cent_mid[icent]->Integral(0, bins) > (totalPR * 0.9099)) continue;
							ctau3D_cut_centBin_mid[icent] = h_decay_pr_cent_mid[icent]->GetBinCenter(bins);
							pr_eff_centBin_mid[icent] = h_decay_pr_cent_mid[icent]->Integral(0, bins) / totalPR;
							np_res_centBin_mid[icent] = h_decay_np_cent_mid[icent]->Integral(0, bins) / totalNP;
						} else {
							if (h_decay_pr_cent_mid[icent]->Integral(0, bins) < (0.97 * totalPR) || 
									h_decay_pr_cent_mid[icent]->Integral(0, bins) > (totalPR * 0.98)) continue;
							ctau3D_cut_centBin_mid[icent] = h_decay_pr_cent_mid[icent]->GetBinCenter(bins);
							pr_eff_centBin_mid[icent] = 1 - h_decay_pr_cent_mid[icent]->Integral(0, bins) / totalPR;
							np_res_centBin_mid[icent] = 1 - h_decay_np_cent_mid[icent]->Integral(0, bins) / totalNP;
						}
					}
					cout << Form("cent %.0f-%.0f: ctau3D cut = %.4f, PR eff = %.3f, NP res = %.3f, integral PR=%.1f NP=%.1f", 
							centLow_bin/2, centHigh_bin/2, ctau3D_cut_centBin_mid[icent], 
							pr_eff_centBin_mid[icent], np_res_centBin_mid[icent], totalPR, totalNP) << endl;

					// Draw efficiency plot for this centrality bin (mid)
					TH1D* h_eff_pr_cent_mid = (TH1D*)h_decay_pr_cent_mid[icent]->Clone(Form("h_eff_pr_cent_mid_%d", icent));
					TH1D* h_eff_np_cent_mid = (TH1D*)h_decay_np_cent_mid[icent]->Clone(Form("h_eff_np_cent_mid_%d", icent));

					for (int ib = 1; ib <= h_eff_pr_cent_mid->GetNbinsX(); ib++) {
						h_eff_pr_cent_mid->SetBinContent(ib, 1.0 - h_decay_pr_cent_mid[icent]->Integral(0, ib) / totalPR);
						h_eff_np_cent_mid->SetBinContent(ib, h_decay_np_cent_mid[icent]->Integral(0, ib) / totalNP);
					}

					TCanvas* c_eff_cent_mid = new TCanvas(Form("c_eff_cent_mid_%d", icent), "", 800, 600);
					c_eff_cent_mid->SetGrid();

					h_eff_np_cent_mid->SetLineColor(kRed);
					h_eff_np_cent_mid->SetLineWidth(2);
					h_eff_np_cent_mid->GetXaxis()->SetTitle("ctau3D (mm)");
					h_eff_np_cent_mid->GetYaxis()->SetTitle("Efficiency");
					h_eff_np_cent_mid->GetYaxis()->SetRangeUser(0, 1.1);
					h_eff_np_cent_mid->GetXaxis()->SetRangeUser(-1, 3);
					h_eff_np_cent_mid->Draw("hist");

					h_eff_pr_cent_mid->SetLineColor(kBlue);
					h_eff_pr_cent_mid->SetLineWidth(2);
					h_eff_pr_cent_mid->Draw("hist same");

					TLine* line_cent_mid = new TLine(ctau3D_cut_centBin_mid[icent], 0, ctau3D_cut_centBin_mid[icent], 1);
					line_cent_mid->SetLineStyle(2);
					line_cent_mid->SetLineWidth(2);
					line_cent_mid->Draw();

					TLegend* leg_cent_mid = new TLegend(0.15, 0.75, 0.45, 0.88);
					leg_cent_mid->SetBorderSize(0);
					leg_cent_mid->SetFillStyle(0);
					leg_cent_mid->AddEntry(h_eff_np_cent_mid, "Non-Prompt J/#psi (CDF)", "l");
					leg_cent_mid->AddEntry(h_eff_pr_cent_mid, "Prompt J/#psi (1-CDF)", "l");
					leg_cent_mid->Draw();

					drawText("6.5 < p_{T}^{#mu#mu} < 40.0 GeV/c", 0.50, 0.85, 1, 18);
					drawText("|y^{#mu#mu}| < 1.6", 0.50, 0.80, 1, 16);
					drawText(Form("Centrality %.0f-%.0f%%", centLow_bin/2, centHigh_bin/2), 0.50, 0.75, 1, 16);
					drawText(Form("ctau3D cut: %.4f mm", ctau3D_cut_centBin_mid[icent]), 0.50, 0.65, 1, 16);
					if (state == 1) {
						drawText(Form("NP J/#psi eff: %.3f", np_res_centBin_mid[icent]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", 1-pr_eff_centBin_mid[icent]), 0.50, 0.55, 2, 16);
					} else {
						drawText(Form("NP J/#psi eff: %.3f", 1-np_res_centBin_mid[icent]), 0.50, 0.60, 4, 16);
						drawText(Form("PR J/#psi res: %.3f", pr_eff_centBin_mid[icent]), 0.50, 0.55, 2, 16);
					}

					TString stateStr = (state == 1) ? "PRMC" : "NPMC";
					c_eff_cent_mid->SaveAs(Form("./figs_ctau3D_v2/eff_centBin_mid_%s_cent%.0f-%.0f.pdf", stateStr.Data(), centLow_bin/2, centHigh_bin/2));
					c_eff_cent_mid->SaveAs(Form("./figs_ctau3D_v2/eff_centBin_mid_%s_cent%.0f-%.0f.png", stateStr.Data(), centLow_bin/2, centHigh_bin/2));

					delete h_eff_pr_cent_mid;
					delete h_eff_np_cent_mid;
					delete c_eff_cent_mid;
					delete line_cent_mid;
					delete leg_cent_mid;
				}
			}

			// Cleanup centrality histograms
			for (int icent = 0; icent < nCentBins_for; icent++) {
				delete h_decay_pr_cent_for[icent];
				delete h_decay_np_cent_for[icent];
			}
			for (int icent = 0; icent < nCentBins_mid; icent++) {
				delete h_decay_pr_cent_mid[icent];
				delete h_decay_np_cent_mid[icent];
			}
		}
		//============================================================================
		// Draw centrality vs ctau3D cut plots
		//============================================================================
		cout << "\n--- Drawing centrality vs ctau3D cut plots ---" << endl;

		// Create TGraphErrors for centrality (forward)
		double centCenter_for[nCentBins_for], centErr_for[nCentBins_for], ctauCentErr_for[nCentBins_for];
		for (int i = 0; i < nCentBins_for; i++) {
			centCenter_for[i] = (centBin_for[i] + centBin_for[i+1]) / 4.0;  // /2 for percentage, /2 for center
			centErr_for[i] = (centBin_for[i+1] - centBin_for[i]) / 4.0;
			ctauCentErr_for[i] = 0.001;
		}
		TGraphErrors *gr_cent_ctau_for = new TGraphErrors(nCentBins_for, centCenter_for, ctau3D_cut_centBin_for, centErr_for, ctauCentErr_for);
		gr_cent_ctau_for->SetName("gr_cent_ctau_for");
		gr_cent_ctau_for->SetTitle(";Centrality (%);ctau3D cut (mm)");
		gr_cent_ctau_for->SetMarkerStyle(20);
		gr_cent_ctau_for->SetMarkerColor(kBlue);
		gr_cent_ctau_for->SetLineColor(kBlue);

		// Create TGraphErrors for centrality (mid)
		double centCenter_mid_arr[nCentBins_mid], centErr_mid_arr[nCentBins_mid], ctauCentErr_mid[nCentBins_mid];
		for (int i = 0; i < nCentBins_mid; i++) {
			centCenter_mid_arr[i] = (centBin_mid[i] + centBin_mid[i+1]) / 4.0;
			centErr_mid_arr[i] = (centBin_mid[i+1] - centBin_mid[i]) / 4.0;
			ctauCentErr_mid[i] = 0.001;
		}
		TGraphErrors *gr_cent_ctau_mid = new TGraphErrors(nCentBins_mid, centCenter_mid_arr, ctau3D_cut_centBin_mid, centErr_mid_arr, ctauCentErr_mid);
		gr_cent_ctau_mid->SetName("gr_cent_ctau_mid");
		gr_cent_ctau_mid->SetTitle(";Centrality (%);ctau3D cut (mm)");
		gr_cent_ctau_mid->SetMarkerStyle(21);
		gr_cent_ctau_mid->SetMarkerColor(kRed);
		gr_cent_ctau_mid->SetLineColor(kRed);

		// Draw centrality plot
		TCanvas *c_cent_ctau = new TCanvas("c_cent_ctau", "", 800, 600);
		c_cent_ctau->SetGrid();

		TH2F *hframe_cent = new TH2F("hframe_cent", ";Centrality (%);ctau3D cut (mm)", 100, 0, 100, 100, -0.1, 0.5);
		hframe_cent->Draw();

		gr_cent_ctau_for->Draw("P same");
		gr_cent_ctau_mid->Draw("P same");

		TLegend *leg_cent = new TLegend(0.55, 0.7, 0.88, 0.88);
		leg_cent->SetBorderSize(0);
		leg_cent->SetFillStyle(0);
		leg_cent->AddEntry(gr_cent_ctau_for, "Forward (1.6 < |y| < 2.4)", "lep");
		leg_cent->AddEntry(gr_cent_ctau_mid, "Mid (|y| < 1.6)", "lep");
		leg_cent->Draw();

		if (state == 1) {
			drawText("Prompt J/#psi selection", 0.2, 0.85, 1, 18);
			drawText("(90% PR efficiency)", 0.2, 0.80, 1, 16);
		} else {
			drawText("Non-Prompt J/#psi selection", 0.2, 0.85, 1, 18);
			drawText("(97-98% PR rejection)", 0.2, 0.80, 1, 16);
		}
		drawText(Form("p_{T} %.1f-%.1f GeV/c", ptLow, ptHigh), 0.2, 0.75, 1, 16);

		c_cent_ctau->Update();
		if (state == 1) {
			c_cent_ctau->SaveAs("./figs_ctau3D_v2/cent_vs_ctau3D_cut_PRMC.pdf");
			c_cent_ctau->SaveAs("./figs_ctau3D_v2/cent_vs_ctau3D_cut_PRMC.png");
		} else {
			c_cent_ctau->SaveAs("./figs_ctau3D_v2/cent_vs_ctau3D_cut_NPMC.pdf");
			c_cent_ctau->SaveAs("./figs_ctau3D_v2/cent_vs_ctau3D_cut_NPMC.png");
		}

		// Save centrality bin ctau3D cut values to ROOT file (separate for forward and mid)
		// Forward rapidity uses pt 3.5-40, Mid rapidity uses pt 6.5-40
		TString rootFileName_centBin_for, rootFileName_centBin_mid;
		if (state == 1) {
			rootFileName_centBin_for = "./roots_ctau3D_v2/ctau3D_cut_centBin_for_PRMC_pt3.5-40.0.root";
			rootFileName_centBin_mid = "./roots_ctau3D_v2/ctau3D_cut_centBin_mid_PRMC_pt6.5-40.0.root";
		} else {
			rootFileName_centBin_for = "./roots_ctau3D_v2/ctau3D_cut_centBin_for_NPMC_pt3.5-40.0.root";
			rootFileName_centBin_mid = "./roots_ctau3D_v2/ctau3D_cut_centBin_mid_NPMC_pt6.5-40.0.root";
		}

		// Save forward rapidity file
		TFile *wf_centBin_for = new TFile(rootFileName_centBin_for, "RECREATE");
		wf_centBin_for->cd();
		gr_cent_ctau_for->Write();
		TH1D *h_ctau_centBin_for = new TH1D("h_ctau_centBin_for", "ctau3D cut vs cent (forward, pt 3.5-40)", nCentBins_for, centBin_for);
		for (int i = 0; i < nCentBins_for; i++) h_ctau_centBin_for->SetBinContent(i+1, ctau3D_cut_centBin_for[i]);
		h_ctau_centBin_for->Write();
		wf_centBin_for->Close();
		cout << "Saved forward centrality bin ctau3D cuts to: " << rootFileName_centBin_for << endl;

		// Save mid rapidity file
		TFile *wf_centBin_mid = new TFile(rootFileName_centBin_mid, "RECREATE");
		wf_centBin_mid->cd();
		gr_cent_ctau_mid->Write();
		TH1D *h_ctau_centBin_mid = new TH1D("h_ctau_centBin_mid", "ctau3D cut vs cent (mid, pt 6.5-40)", nCentBins_mid, centBin_mid);
		for (int i = 0; i < nCentBins_mid; i++) h_ctau_centBin_mid->SetBinContent(i+1, ctau3D_cut_centBin_mid[i]);
		h_ctau_centBin_mid->Write();
		wf_centBin_mid->Close();
		cout << "Saved mid centrality bin ctau3D cuts to: " << rootFileName_centBin_mid << endl;

		//============================================================================
		// Draw pt vs ctau3D cut plots
		//============================================================================
		cout << "\n--- Drawing pt vs ctau3D cut plots ---" << endl;

		// Create TGraphErrors for forward rapidity
		double ptCenter_for[5], ptErr_for[5], ctauErr_for[5];
		for (int i = 0; i < nPtBins_for; i++) {
			ptCenter_for[i] = (ptBin_for[i] + ptBin_for[i+1]) / 2.0;
			ptErr_for[i] = (ptBin_for[i+1] - ptBin_for[i]) / 2.0;
			ctauErr_for[i] = 0.001;  // small error for visualization
		}
		TGraphErrors *gr_pt_ctau_for = new TGraphErrors(nPtBins_for, ptCenter_for, ctau3D_cut_ptBin_for, ptErr_for, ctauErr_for);
		gr_pt_ctau_for->SetName("gr_pt_ctau_for");
		gr_pt_ctau_for->SetTitle(";p_{T} (GeV/c);ctau3D cut (mm)");
		gr_pt_ctau_for->SetMarkerStyle(20);
		gr_pt_ctau_for->SetMarkerColor(kBlue);
		gr_pt_ctau_for->SetLineColor(kBlue);

		// Create TGraphErrors for mid rapidity
		double ptCenter_mid[7], ptErr_mid[7], ctauErr_mid[7];
		for (int i = 0; i < nPtBins_mid; i++) {
			ptCenter_mid[i] = (ptBin_mid[i] + ptBin_mid[i+1]) / 2.0;
			ptErr_mid[i] = (ptBin_mid[i+1] - ptBin_mid[i]) / 2.0;
			ctauErr_mid[i] = 0.001;
		}
		TGraphErrors *gr_pt_ctau_mid = new TGraphErrors(nPtBins_mid, ptCenter_mid, ctau3D_cut_ptBin_mid, ptErr_mid, ctauErr_mid);
		gr_pt_ctau_mid->SetName("gr_pt_ctau_mid");
		gr_pt_ctau_mid->SetTitle(";p_{T} (GeV/c);ctau3D cut (mm)");
		gr_pt_ctau_mid->SetMarkerStyle(21);
		gr_pt_ctau_mid->SetMarkerColor(kRed);
		gr_pt_ctau_mid->SetLineColor(kRed);

		// Draw combined plot
		TCanvas *c_pt_ctau = new TCanvas("c_pt_ctau", "", 800, 600);
		c_pt_ctau->SetGrid();

		TH2F *hframe = new TH2F("hframe", ";p_{T} (GeV/c);ctau3D cut (mm)", 100, 0, 45, 100, -0.1, 0.5);
		hframe->Draw();

		gr_pt_ctau_for->Draw("P same");
		gr_pt_ctau_mid->Draw("P same");

		TLegend *leg_pt = new TLegend(0.55, 0.7, 0.88, 0.88);
		leg_pt->SetBorderSize(0);
		leg_pt->SetFillStyle(0);
		leg_pt->AddEntry(gr_pt_ctau_for, "Forward (1.6 < |y| < 2.4)", "lep");
		leg_pt->AddEntry(gr_pt_ctau_mid, "Mid (|y| < 1.6)", "lep");
		leg_pt->Draw();

		if (state == 1) {
			drawText("Prompt J/#psi selection", 0.2, 0.85, 1, 18);
			drawText("(90% PR efficiency)", 0.2, 0.80, 1, 16);
		} else {
			drawText("Non-Prompt J/#psi selection", 0.2, 0.85, 1, 18);
			drawText("(97-98% PR rejection)", 0.2, 0.80, 1, 16);
		}
		drawText(Form("Centrality %.0f-%.0f%s", cLow/2, cHigh/2, perc.Data()), 0.2, 0.75, 1, 16);

		c_pt_ctau->Update();
		if (state == 1) {
			c_pt_ctau->SaveAs("./figs_ctau3D_v2/pt_vs_ctau3D_cut_PRMC.pdf");
			c_pt_ctau->SaveAs("./figs_ctau3D_v2/pt_vs_ctau3D_cut_PRMC.png");
		} else {
			c_pt_ctau->SaveAs("./figs_ctau3D_v2/pt_vs_ctau3D_cut_NPMC.pdf");
			c_pt_ctau->SaveAs("./figs_ctau3D_v2/pt_vs_ctau3D_cut_NPMC.png");
		}

		// Save pt bin ctau3D cut values to ROOT file
		TString rootFileName_ptBin;
		if (state == 1) {
			rootFileName_ptBin = Form("./roots_ctau3D_v2/ctau3D_cut_ptBin_PRMC_cent%.0f-%.0f.root", cLow, cHigh);
		} else {
			rootFileName_ptBin = Form("./roots_ctau3D_v2/ctau3D_cut_ptBin_NPMC_cent%.0f-%.0f.root", cLow, cHigh);
		}
		TFile *wf_ptBin = new TFile(rootFileName_ptBin, "RECREATE");
		wf_ptBin->cd();
		gr_pt_ctau_for->Write();
		gr_pt_ctau_mid->Write();

		// Save as histograms too
		TH1D *h_ctau_ptBin_for = new TH1D("h_ctau_ptBin_for", "ctau3D cut vs pt (forward)", nPtBins_for, ptBin_for);
		TH1D *h_ctau_ptBin_mid = new TH1D("h_ctau_ptBin_mid", "ctau3D cut vs pt (mid)", nPtBins_mid, ptBin_mid);
		for (int i = 0; i < nPtBins_for; i++) h_ctau_ptBin_for->SetBinContent(i+1, ctau3D_cut_ptBin_for[i]);
		for (int i = 0; i < nPtBins_mid; i++) h_ctau_ptBin_mid->SetBinContent(i+1, ctau3D_cut_ptBin_mid[i]);
		h_ctau_ptBin_for->Write();
		h_ctau_ptBin_mid->Write();

		wf_ptBin->Close();
		cout << "Saved pt bin ctau3D cuts to: " << rootFileName_ptBin << endl;

		//============================================================================
		// End of ctau3D cut calculation
		//============================================================================


		//cout << "count " << count << endl;
		//cout << "counttnp " << counttnp << endl;

		// Divide
		// draw same
		// gROOT->Macro("~/rootlogon.C");
	}
}
