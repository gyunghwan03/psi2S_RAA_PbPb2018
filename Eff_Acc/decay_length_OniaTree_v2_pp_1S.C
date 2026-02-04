#include <iostream>
#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tnp_weight_pp.h"
#include <TAttMarker.h>

using namespace std;

void decay_length_OniaTree_v2_pp_1S(
		int state = 1,
		int type = 1, // 1 : pt dependent ctau cut, 2 : cent dependent ctau cut, 3 : whole pt and cent dependent ctau cut
		bool isTnP = true, int isPtWeight = 0, // 0 : nominal , 1 : up, -1 : down
		bool isMC = true,
		float ptLow = 3.0, float ptHigh = 40.0,
		float yLow = 0.0, float yHigh = 2.4,
		int PRMC = 0 // 0: Prompt, 1: Non-Prompt
		)
{

	gStyle->SetOptStat(0);
	int kTrigSel = 3; // pp Jpsi=3

	float muPtCut = 0;

	TString ptSys;
	if (isPtWeight == 0)
		ptSys = "nomi";
	else if (isPtWeight == 1)
		ptSys = "up";
	else if (isPtWeight == -1)
		ptSys = "down";

	// jpsi mass range
	float massLow = 2.9;
	float massHigh = 3.3;

	double min = 0;
	double max = ptHigh;

	TString kineLabel = Form("state%d_pt%.1f-%.1f_y%.1f-%.1f_m%.1f-%.1f", 
			state, ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

	float pos_x_mass = 0.55;
	float pos_y = 0.65;
	float pos_y_diff = 0.071;
	int text_color = 1;
	float text_size = 16;
	TString perc = "%";

	// Create directories
	gSystem->mkdir("./roots_1S_pp/");
	gSystem->mkdir("./roots_1S_pp/decayL");
	gSystem->mkdir("./roots_1S_pp/decayL/PRMC");
	gSystem->mkdir("./roots_1S_pp/decayL/NPMC");

	gSystem->mkdir("./figs_1S_pp/");
	gSystem->mkdir("./figs_1S_pp/decayL");
	gSystem->mkdir("./figs_1S_pp/decayL/PRMC");
	gSystem->mkdir("./figs_1S_pp/decayL/NPMC");

	// Input files - pp collision
	TString fPRMC_path = "/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root"; // pp prompt
	TString fNPMC_path = "/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root"; // pp non-prompt

	// Separate TChains for PRMC and NPMC (for ctau3D cut calculation)
	TChain *treePRMC = new TChain("hionia/myTree");
	TChain *treeNPMC = new TChain("hionia/myTree");
	treePRMC->Add(fPRMC_path.Data());
	treeNPMC->Add(fNPMC_path.Data());

	// Combined tree for main analysis
	TChain *mytree = new TChain("hionia/myTree");
	mytree->Add(fPRMC_path.Data());
	if(state==2) mytree->Add(fNPMC_path.Data());

	cout << "Total entries : " << mytree->GetEntries() << endl;
	cout << "PRMC entries : " << treePRMC->GetEntries() << endl;
	cout << "NPMC entries : " << treeNPMC->GetEntries() << endl;

	// pT reweighting function
	TFile *fPtW1;
	TFile *fPtW2;
	if(state==1){
		fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
		fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
	}
	else if(state==2){
		fPtW1 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251103.root", "read");
		fPtW2 = new TFile("../compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251103.root", "read");
	}

	TF1 *fptw1 = (TF1 *)fPtW1->Get("dataMC_Ratio1");
	TF1 *fptw2 = (TF1 *)fPtW2->Get("dataMC_Ratio1");

	double ptBin_all[18] = {0, 3., 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 22.5, 25, 27.5, 30, 50};
	double ptBin_for[6] = {0, 3.5, 6.5, 9, 12, 40};
	double ptBin_mid[8] = {0, 6.5, 9, 12, 15, 20, 25, 40};

	float nbins = 4000;
	float xmin = -1; 
	float xmax = 3;

	TH1D *hpt_ctau_1 = new TH1D("hpt_ctau_1", "hpt_ctau_1", 5, ptBin_for);
	TH1D *hpt_ctau_2 = new TH1D("hpt_ctau_2", "hpt_ctau_2", 7, ptBin_mid);

	TH1D *hInt_ctau_1 = new TH1D("hInt_ctau_1", "", 1, 0, 50);
	TH1D *hInt_ctau_2 = new TH1D("hInt_ctau_2", "", 1, 0, 50);

	TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
	TH1D* h_massCut = new TH1D("h_massCut",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
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
	////////////////////////////// for systematics  ////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	// parameters for mid-rapidity
	double p10 = 0.0, p11 = 0.0, p12 = 0.0, p13 = 0.0;
	// parameters for forward
	double p20 = 0.0, p21 = 0.0, p22 = 0.0, p23 = 0.0;
	// errors for mid-rapidity
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

	TF1 *fptw1sys = new TF1("fptw1sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 3.0, 50.0);
	TF1 *fptw2sys = new TF1("fptw2sys", "( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 3.0, 50.0);

	cout << "p10 : " << p10 << ", e10 : " << e10 << ", p20 : " << p20 << ", e20 : " << e20 << endl;

	// sys for nominal
	if (isPtWeight == 0)
	{
		p10 = p10; p11 = p11; p12 = p12; p13 = p13;
		p20 = p20; p21 = p21; p22 = p22; p23 = p23;
	}

	// sys for up
	if (isPtWeight == 1)
	{
		p10 = p10 + e10; p11 = p11 + e11; p12 = p12 + e12; p13 = p13 + e13;
		p20 = p20 + e20; p21 = p21 + e21; p22 = p22 + e22; p23 = p23 + e23;
	}

	// sys for down
	if (isPtWeight == -1)
	{
		p10 = p10 - e10; p11 = p11 - e11; p12 = p12 - e12; p13 = p13 - e13;
		p20 = p20 - e20; p21 = p21 - e21; p22 = p22 - e22; p23 = p23 - e23;
	}

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

	//============================================================================
	// Calculate ctau3D cut (lcutv) from PRMC and NPMC
	//============================================================================
	cout << "============================================" << endl;
	cout << "Calculating ctau3D cut from PRMC and NPMC..." << endl;
	cout << "============================================" << endl;

	const int maxBranchSize_ctau = 1000;

	// Set up branches for PRMC tree
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

	// Set up branches for NPMC tree
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

	// Fill decay length histograms for PRMC
	TLorentzVector *JP_Reco_pr = new TLorentzVector;
	TLorentzVector *mupl_Reco_pr = new TLorentzVector;
	TLorentzVector *mumi_Reco_pr = new TLorentzVector;

	cout << "Filling PRMC histogram..." << endl;
	Long64_t nentries_pr = treePRMC->GetEntries();
	for (Long64_t i = 0; i < nentries_pr; i++)
	{
		if (i % 10000000 == 0) cout << "PRMC: " << i << "/" << nentries_pr << endl;
		treePRMC->GetEntry(i);

		if (!(HLTriggers_pr & ((ULong64_t)pow(2, kTrigSel)))) continue;

		for (int iQQ = 0; iQQ < Reco_QQ_size_pr; iQQ++)
		{
			JP_Reco_pr = (TLorentzVector *)Reco_QQ_4mom_pr->At(iQQ);
			
			double mass_pr = JP_Reco_pr->M();
			double pt_pr = JP_Reco_pr->Pt();
			double y_pr = JP_Reco_pr->Rapidity();

			if (!(mass_pr >= massLow && mass_pr <= massHigh)) continue;
			if (!(pt_pr >= ptLow && pt_pr <= ptHigh)) continue;
			if (!(TMath::Abs(y_pr) >= yLow && TMath::Abs(y_pr) <= yHigh)) continue;
			if (Reco_QQ_sign_pr[iQQ] != 0) continue;
			if (Reco_QQ_VtxProb_pr[iQQ] < 0.01) continue;

			int iMupl_pr = Reco_QQ_mupl_idx_pr[iQQ];
			int iMumi_pr = Reco_QQ_mumi_idx_pr[iQQ];
			mupl_Reco_pr = (TLorentzVector *)Reco_mu_4mom_pr->At(iMupl_pr);
			mumi_Reco_pr = (TLorentzVector *)Reco_mu_4mom_pr->At(iMumi_pr);

			double eta1_pr = mupl_Reco_pr->Eta();
			double eta2_pr = mumi_Reco_pr->Eta();

			// Muon acceptance cuts (pp)
			if (!IsAcceptanceQQ(mupl_Reco_pr->Pt(), TMath::Abs(eta1_pr))) continue;
			if (!IsAcceptanceQQ(mumi_Reco_pr->Pt(), TMath::Abs(eta2_pr))) continue;
			if (!(TMath::Abs(eta1_pr) < 2.4 && TMath::Abs(eta2_pr) < 2.4)) continue;

			// Muon ID cuts
			bool passMuonTypePl_pr = true;
			passMuonTypePl_pr = passMuonTypePl_pr && (Reco_mu_SelectionType_pr[iMupl_pr] & ((int)pow(2, 1)));
			passMuonTypePl_pr = passMuonTypePl_pr && (Reco_mu_SelectionType_pr[iMupl_pr] & ((int)pow(2, 3)));

			bool passMuonTypeMi_pr = true;
			passMuonTypeMi_pr = passMuonTypeMi_pr && (Reco_mu_SelectionType_pr[iMumi_pr] & ((int)pow(2, 1)));
			passMuonTypeMi_pr = passMuonTypeMi_pr && (Reco_mu_SelectionType_pr[iMumi_pr] & ((int)pow(2, 3)));

			bool muplSoft_pr = (
				(Reco_mu_nTrkWMea_pr[iMupl_pr] > 5) &&
				(Reco_mu_nPixWMea_pr[iMupl_pr] > 0) &&
				(TMath::Abs(Reco_mu_dxy_pr[iMupl_pr]) < 0.3) &&
				(TMath::Abs(Reco_mu_dz_pr[iMupl_pr]) < 20.) &&
				passMuonTypePl_pr
			);

			bool mumiSoft_pr = (
				(Reco_mu_nTrkWMea_pr[iMumi_pr] > 5) &&
				(Reco_mu_nPixWMea_pr[iMumi_pr] > 0) &&
				(TMath::Abs(Reco_mu_dxy_pr[iMumi_pr]) < 0.3) &&
				(TMath::Abs(Reco_mu_dz_pr[iMumi_pr]) < 20.) &&
				passMuonTypeMi_pr
			);

			if (!(muplSoft_pr && mumiSoft_pr)) continue;

			h_decayPRMC->Fill(Reco_QQ_ctau3D_pr[iQQ]);
		}
	}

	// Fill decay length histograms for NPMC
	TLorentzVector *JP_Reco_np = new TLorentzVector;
	TLorentzVector *mupl_Reco_np = new TLorentzVector;
	TLorentzVector *mumi_Reco_np = new TLorentzVector;

	cout << "Filling NPMC histogram..." << endl;
	Long64_t nentries_np = treeNPMC->GetEntries();
	for (Long64_t i = 0; i < nentries_np; i++)
	{
		if (i % 10000000 == 0) cout << "NPMC: " << i << "/" << nentries_np << endl;
		treeNPMC->GetEntry(i);

		if (!(HLTriggers_np & ((ULong64_t)pow(2, kTrigSel)))) continue;

		for (int iQQ = 0; iQQ < Reco_QQ_size_np; iQQ++)
		{
			JP_Reco_np = (TLorentzVector *)Reco_QQ_4mom_np->At(iQQ);
			
			double mass_np = JP_Reco_np->M();
			double pt_np = JP_Reco_np->Pt();
			double y_np = JP_Reco_np->Rapidity();

			if (!(mass_np >= massLow && mass_np <= massHigh)) continue;
			if (!(pt_np >= ptLow && pt_np <= ptHigh)) continue;
			if (!(TMath::Abs(y_np) >= yLow && TMath::Abs(y_np) <= yHigh)) continue;
			if (Reco_QQ_sign_np[iQQ] != 0) continue;
			if (Reco_QQ_VtxProb_np[iQQ] < 0.01) continue;

			int iMupl_np = Reco_QQ_mupl_idx_np[iQQ];
			int iMumi_np = Reco_QQ_mumi_idx_np[iQQ];
			mupl_Reco_np = (TLorentzVector *)Reco_mu_4mom_np->At(iMupl_np);
			mumi_Reco_np = (TLorentzVector *)Reco_mu_4mom_np->At(iMumi_np);

			double eta1_np = mupl_Reco_np->Eta();
			double eta2_np = mumi_Reco_np->Eta();

			// Muon acceptance cuts (pp)
			if (!IsAcceptanceQQ(mupl_Reco_np->Pt(), TMath::Abs(eta1_np))) continue;
			if (!IsAcceptanceQQ(mumi_Reco_np->Pt(), TMath::Abs(eta2_np))) continue;
			if (!(TMath::Abs(eta1_np) < 2.4 && TMath::Abs(eta2_np) < 2.4)) continue;

			// Muon ID cuts
			bool passMuonTypePl_np = true;
			passMuonTypePl_np = passMuonTypePl_np && (Reco_mu_SelectionType_np[iMupl_np] & ((int)pow(2, 1)));
			passMuonTypePl_np = passMuonTypePl_np && (Reco_mu_SelectionType_np[iMupl_np] & ((int)pow(2, 3)));

			bool passMuonTypeMi_np = true;
			passMuonTypeMi_np = passMuonTypeMi_np && (Reco_mu_SelectionType_np[iMumi_np] & ((int)pow(2, 1)));
			passMuonTypeMi_np = passMuonTypeMi_np && (Reco_mu_SelectionType_np[iMumi_np] & ((int)pow(2, 3)));

			bool muplSoft_np = (
				(Reco_mu_nTrkWMea_np[iMupl_np] > 5) &&
				(Reco_mu_nPixWMea_np[iMupl_np] > 0) &&
				(TMath::Abs(Reco_mu_dxy_np[iMupl_np]) < 0.3) &&
				(TMath::Abs(Reco_mu_dz_np[iMupl_np]) < 20.) &&
				passMuonTypePl_np
			);

			bool mumiSoft_np = (
				(Reco_mu_nTrkWMea_np[iMumi_np] > 5) &&
				(Reco_mu_nPixWMea_np[iMumi_np] > 0) &&
				(TMath::Abs(Reco_mu_dxy_np[iMumi_np]) < 0.3) &&
				(TMath::Abs(Reco_mu_dz_np[iMumi_np]) < 20.) &&
				passMuonTypeMi_np
			);

			if (!(muplSoft_np && mumiSoft_np)) continue;

			h_decayNPMC->Fill(Reco_QQ_ctau3D_np[iQQ]);
		}
	}

	cout << "PRMC Tree: " << treePRMC->GetEntries() << ", by total cut: " << h_decayPRMC->GetEntries() << ", " << h_decayPRMC->GetEntries() / treePRMC->GetEntries() * 100 << "%" << endl;
	cout << "NPMC Tree: " << treeNPMC->GetEntries() << ", by total cut: " << h_decayNPMC->GetEntries() << ", " << h_decayNPMC->GetEntries() / treeNPMC->GetEntries() * 100 << "%" << endl;

	// Calculate efficiency curves and find cut value
	int totalPRMC = h_decayPRMC->GetEntries();
	int totalNPMC = h_decayNPMC->GetEntries();

	for(int bins=0; bins<h_decayPRMC->GetNbinsX(); bins++){
		if(PRMC==0){
			h_deffPRMC->SetBinContent(bins,h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
			h_deffNPMC->SetBinContent(bins,1-(h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries()));

			if(h_decayPRMC->Integral(0,bins)<(0.899*totalPRMC)||h_decayPRMC->Integral(0,bins)>(totalPRMC*0.9099)) continue;
			cout<<bins<<"th bin"<<", ljpsi: "<<h_decayPRMC->GetBinCenter(bins)<<", PR eff: "<<h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries()*100<<"%"<<endl;
			lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
			lcuth = new TLine(xmin, h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries(), h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
			lresi = new TLine(xmin, h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries(), h_decayNPMC->GetBinCenter(bins), h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
		}
		else if(PRMC==1){
			h_deffNPMC->SetBinContent(bins,h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
			h_deffPRMC->SetBinContent(bins,1-(h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries()));

			if(h_decayPRMC->Integral(0,bins)<(0.97*totalPRMC)||h_decayPRMC->Integral(0,bins)>(totalPRMC*0.98)) continue;
			cout<<bins<<"th bin"<<", ljpsi: "<<h_decayPRMC->GetBinCenter(bins)<<", NP eff: "<<(1-h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries())*100<<"%"<<endl;
			lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
			lcuth = new TLine(xmin, h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries(), h_decayNPMC->GetBinCenter(bins), h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
			lresi = new TLine(xmin, h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries(), h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
		}
	}

	ctau3D_cut = lcutv->GetX1();
	if(PRMC==0){
		pr_eff = lcuth->GetY1();
		np_res = 1 - lresi->GetY1();
	}
	else if(PRMC==1){
		pr_eff = 1-lresi->GetY1();
		np_res = 1-lcuth->GetY1();
	}

	cout << "============================================" << endl;
	cout << "ctau3D_cut: " << ctau3D_cut << endl;
	cout << "PR eff: " << pr_eff << endl;
	cout << "NP res: " << np_res << endl;
	cout << "============================================" << endl;

	// Draw and save plots
	TCanvas* c_decayL = new TCanvas("c_decayL","",600,600);
	h_deffPRMC->Draw("l");
	h_deffNPMC->Draw("same");
	h_deffNPMC->SetLineColor(kRed+2);
	jumSun(xmin,1,xmax,1,1,1);
	lcutv->Draw("same");
	if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
	if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff*0.5,text_color,text_size);
	else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
	drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*1.5,text_color,text_size);
	if(PRMC==0){
		drawText(Form("L_{J/psi} cut: %.4f", lcutv->GetX1()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
		drawText(Form("PR J/psi eff: %.3f", lcuth->GetY1()),pos_x_mass,pos_y-pos_y_diff*3.5,text_color,text_size);
		drawText(Form("NP J/psi res: %.3f", 1 - lresi->GetY1()),pos_x_mass,pos_y-pos_y_diff*4, kRed+2,text_size);
		c_decayL->Update();
		c_decayL->SaveAs(Form("figs_1S_pp/decayL/PRMC/decay_%s.pdf",kineLabel.Data()));
	}
	else if(PRMC==1){
		drawText(Form("L_{J/psi} cut: %.4f", lcutv->GetX1()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
		drawText(Form("NP J/psi eff: %.3f", 1-lcuth->GetY1()),pos_x_mass,pos_y-pos_y_diff*3.5,kRed+2,text_size);
		drawText(Form("PR J/psi res: %.3f", 1-lresi->GetY1()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
		c_decayL->Update();
		c_decayL->SaveAs(Form("figs_1S_pp/decayL/NPMC/decay_%s.pdf",kineLabel.Data()));
	}

	// Save histograms
	TFile *wf;
	if(PRMC==0){
		wf = new TFile(Form("roots_1S_pp/decayL/PRMC/decay_hist_%s.root", kineLabel.Data()),"recreate");
	}
	else if(PRMC==1){
		wf = new TFile(Form("roots_1S_pp/decayL/NPMC/decay_hist_%s.root", kineLabel.Data()),"recreate");
	}
	wf->cd();

	h_deffPRMC->Write();
	h_deffNPMC->Write();
	h_decayPRMC->Write();
	h_decayNPMC->Write();

	// save the l_jpsi_cut
	auto h_Lcut = new TH1D("h_Lcut", "l_jpsi_cut", 1, 0, 1);
	h_Lcut->SetBinContent(1, lcutv->GetX1());
	h_Lcut->Write();

	auto h_eff = new TH1D("h_eff", "eff", 1, 0, 1);
	if(PRMC==0) h_eff->SetBinContent(1, lcuth->GetY1());
	else if(PRMC==1) h_eff->SetBinContent(1, 1-lcuth->GetY1());
	h_eff->Write();

	auto h_res = new TH1D("h_res", "res", 1, 0, 1);
	if(PRMC==0) h_res->SetBinContent(1, lresi->GetY1());
	else if(PRMC==1) h_res->SetBinContent(1, 1-lresi->GetY1());
	h_res->Write();

	wf->Close();

	cout << "============================================" << endl;
	cout << "Decay length analysis completed!" << endl;
	cout << "Output saved to: " << wf->GetName() << endl;
	cout << "============================================" << endl;

	//============================================================================
	// Calculate ctau3D cut for each pt bin
	//============================================================================
	if (type == 1 || type == 3) {
		cout << endl;
		cout << "============================================" << endl;
		cout << "Calculating ctau3D cut for each pt bin..." << endl;
		cout << "============================================" << endl;

		// Define pt bins based on rapidity
		int nPtBins_for = 5;  // for forward (1.6-2.4)
		int nPtBins_mid = 7;  // for mid (0-1.6)

		// Arrays to store ctau3D cut values for each pt bin
		double ctau3D_cut_ptBin_for[5] = {0};
		double ctau3D_cut_ptBin_mid[7] = {0};
		double pr_eff_ptBin_for[5] = {0};
		double pr_eff_ptBin_mid[7] = {0};
		double np_res_ptBin_for[5] = {0};
		double np_res_ptBin_mid[7] = {0};

		// Integrated-bin results (forward/mid rapidity)
		double ctau3D_cut_int_for = 0;
		double ctau3D_cut_int_mid = 0;
		double pr_eff_int_for = 0;
		double pr_eff_int_mid = 0;
		double np_res_int_for = 0;
		double np_res_int_mid = 0;

		// Create histograms for each pt bin
		TH1D* h_decay_pr_pt_for[5];
		TH1D* h_decay_np_pt_for[5];
		TH1D* h_decay_pr_pt_mid[7];
		TH1D* h_decay_np_pt_mid[7];

		// Integrated histograms (fixed pT ranges)
		TH1D* h_decay_pr_int_for = new TH1D("h_decay_pr_int_for", "", nbins, xmin, xmax);
		TH1D* h_decay_np_int_for = new TH1D("h_decay_np_int_for", "", nbins, xmin, xmax);
		TH1D* h_decay_pr_int_mid = new TH1D("h_decay_pr_int_mid", "", nbins, xmin, xmax);
		TH1D* h_decay_np_int_mid = new TH1D("h_decay_np_int_mid", "", nbins, xmin, xmax);

		for (int ipt = 0; ipt < nPtBins_for; ipt++) {
			h_decay_pr_pt_for[ipt] = new TH1D(Form("h_decay_pr_pt_for_%d", ipt), "", nbins, xmin, xmax);
			h_decay_np_pt_for[ipt] = new TH1D(Form("h_decay_np_pt_for_%d", ipt), "", nbins, xmin, xmax);
		}
		for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
			h_decay_pr_pt_mid[ipt] = new TH1D(Form("h_decay_pr_pt_mid_%d", ipt), "", nbins, xmin, xmax);
			h_decay_np_pt_mid[ipt] = new TH1D(Form("h_decay_np_pt_mid_%d", ipt), "", nbins, xmin, xmax);
		}

		// Fill PRMC histograms for ALL pt bins
		cout << "Filling PRMC histograms for all pt bins..." << endl;
		for (Long64_t iev = 0; iev < treePRMC->GetEntries(); ++iev) {
			if (iev % 100000 == 0) cout << "PRMC pt-binned: " << iev << " / " << treePRMC->GetEntries() << endl;
			treePRMC->GetEntry(iev);

			if (!(HLTriggers_pr & ((ULong64_t)pow(2, kTrigSel)))) continue;

			for (int iQQ = 0; iQQ < Reco_QQ_size_pr; iQQ++) {
				TLorentzVector *JP_pr = (TLorentzVector *)Reco_QQ_4mom_pr->At(iQQ);
				TLorentzVector *mupl_pr = (TLorentzVector *)Reco_mu_4mom_pr->At(Reco_QQ_mupl_idx_pr[iQQ]);
				TLorentzVector *mumi_pr = (TLorentzVector *)Reco_mu_4mom_pr->At(Reco_QQ_mumi_idx_pr[iQQ]);

				double mass_pr = JP_pr->M();
				double pt_pr = JP_pr->Pt();
				double y_pr = JP_pr->Rapidity();

				if (!(mass_pr >= massLow && mass_pr <= massHigh)) continue;
				if (!(TMath::Abs(y_pr) >= yLow && TMath::Abs(y_pr) <= yHigh)) continue;
				if (Reco_QQ_sign_pr[iQQ] != 0) continue;
				if (Reco_QQ_VtxProb_pr[iQQ] < 0.01) continue;

				int iMupl_pr = Reco_QQ_mupl_idx_pr[iQQ];
				int iMumi_pr = Reco_QQ_mumi_idx_pr[iQQ];

				double eta1_pr = mupl_pr->Eta();
				double eta2_pr = mumi_pr->Eta();

				// Muon acceptance cuts (pp)
				if (!IsAcceptanceQQ(mupl_pr->Pt(), TMath::Abs(eta1_pr))) continue;
				if (!IsAcceptanceQQ(mumi_pr->Pt(), TMath::Abs(eta2_pr))) continue;
				if (!(TMath::Abs(eta1_pr) < 2.4 && TMath::Abs(eta2_pr) < 2.4)) continue;

				// Muon ID cuts
				bool passMuonTypePl_pr = true;
				passMuonTypePl_pr = passMuonTypePl_pr && (Reco_mu_SelectionType_pr[iMupl_pr] & ((int)pow(2, 1)));
				passMuonTypePl_pr = passMuonTypePl_pr && (Reco_mu_SelectionType_pr[iMupl_pr] & ((int)pow(2, 3)));

				bool passMuonTypeMi_pr = true;
				passMuonTypeMi_pr = passMuonTypeMi_pr && (Reco_mu_SelectionType_pr[iMumi_pr] & ((int)pow(2, 1)));
				passMuonTypeMi_pr = passMuonTypeMi_pr && (Reco_mu_SelectionType_pr[iMumi_pr] & ((int)pow(2, 3)));

				bool muplSoft_pr = (
					(Reco_mu_nTrkWMea_pr[iMupl_pr] > 5) &&
					(Reco_mu_nPixWMea_pr[iMupl_pr] > 0) &&
					(TMath::Abs(Reco_mu_dxy_pr[iMupl_pr]) < 0.3) &&
					(TMath::Abs(Reco_mu_dz_pr[iMupl_pr]) < 20.) &&
					passMuonTypePl_pr
				);

				bool mumiSoft_pr = (
					(Reco_mu_nTrkWMea_pr[iMumi_pr] > 5) &&
					(Reco_mu_nPixWMea_pr[iMumi_pr] > 0) &&
					(TMath::Abs(Reco_mu_dxy_pr[iMumi_pr]) < 0.3) &&
					(TMath::Abs(Reco_mu_dz_pr[iMumi_pr]) < 20.) &&
					passMuonTypeMi_pr
				);

				if (!(muplSoft_pr && mumiSoft_pr)) continue;

				double rap = TMath::Abs(y_pr);
				double ctau = Reco_QQ_ctau3D_pr[iQQ];

				// Fill forward rapidity histograms (1.6-2.4)
				if (rap >= 1.6 && rap < 2.4) {
					for (int ipt = 0; ipt < nPtBins_for; ipt++) {
						if (pt_pr >= ptBin_for[ipt] && pt_pr < ptBin_for[ipt+1]) {
							h_decay_pr_pt_for[ipt]->Fill(ctau);
						}
					}
					if (pt_pr > 3.5 && pt_pr < 40.0) {
						h_decay_pr_int_for->Fill(ctau);
					}
				}
				// Fill mid rapidity histograms (0-1.6)
				if (rap < 1.6) {
					for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
						if (pt_pr >= ptBin_mid[ipt] && pt_pr < ptBin_mid[ipt+1]) {
							h_decay_pr_pt_mid[ipt]->Fill(ctau);
						}
					}
					if (pt_pr > 6.5 && pt_pr < 40.0) {
						h_decay_pr_int_mid->Fill(ctau);
					}
				}
			}
		}

		// Fill NPMC histograms for ALL pt bins
		cout << "Filling NPMC histograms for all pt bins..." << endl;
		for (Long64_t iev = 0; iev < treeNPMC->GetEntries(); ++iev) {
			if (iev % 100000 == 0) cout << "NPMC pt-binned: " << iev << " / " << treeNPMC->GetEntries() << endl;
			treeNPMC->GetEntry(iev);

			if (!(HLTriggers_np & ((ULong64_t)pow(2, kTrigSel)))) continue;

			for (int iQQ = 0; iQQ < Reco_QQ_size_np; iQQ++) {
				TLorentzVector *JP_np = (TLorentzVector *)Reco_QQ_4mom_np->At(iQQ);
				TLorentzVector *mupl_np = (TLorentzVector *)Reco_mu_4mom_np->At(Reco_QQ_mupl_idx_np[iQQ]);
				TLorentzVector *mumi_np = (TLorentzVector *)Reco_mu_4mom_np->At(Reco_QQ_mumi_idx_np[iQQ]);

				double mass_np = JP_np->M();
				double pt_np = JP_np->Pt();
				double y_np = JP_np->Rapidity();

				if (!(mass_np >= massLow && mass_np <= massHigh)) continue;
				if (!(TMath::Abs(y_np) >= yLow && TMath::Abs(y_np) <= yHigh)) continue;
				if (Reco_QQ_sign_np[iQQ] != 0) continue;
				if (Reco_QQ_VtxProb_np[iQQ] < 0.01) continue;

				int iMupl_np = Reco_QQ_mupl_idx_np[iQQ];
				int iMumi_np = Reco_QQ_mumi_idx_np[iQQ];

				double eta1_np = mupl_np->Eta();
				double eta2_np = mumi_np->Eta();

				// Muon acceptance cuts (pp)
				if (!IsAcceptanceQQ(mupl_np->Pt(), TMath::Abs(eta1_np))) continue;
				if (!IsAcceptanceQQ(mumi_np->Pt(), TMath::Abs(eta2_np))) continue;
				if (!(TMath::Abs(eta1_np) < 2.4 && TMath::Abs(eta2_np) < 2.4)) continue;

				// Muon ID cuts
				bool passMuonTypePl_np = true;
				passMuonTypePl_np = passMuonTypePl_np && (Reco_mu_SelectionType_np[iMupl_np] & ((int)pow(2, 1)));
				passMuonTypePl_np = passMuonTypePl_np && (Reco_mu_SelectionType_np[iMupl_np] & ((int)pow(2, 3)));

				bool passMuonTypeMi_np = true;
				passMuonTypeMi_np = passMuonTypeMi_np && (Reco_mu_SelectionType_np[iMumi_np] & ((int)pow(2, 1)));
				passMuonTypeMi_np = passMuonTypeMi_np && (Reco_mu_SelectionType_np[iMumi_np] & ((int)pow(2, 3)));

				bool muplSoft_np = (
					(Reco_mu_nTrkWMea_np[iMupl_np] > 5) &&
					(Reco_mu_nPixWMea_np[iMupl_np] > 0) &&
					(TMath::Abs(Reco_mu_dxy_np[iMupl_np]) < 0.3) &&
					(TMath::Abs(Reco_mu_dz_np[iMupl_np]) < 20.) &&
					passMuonTypePl_np
				);

				bool mumiSoft_np = (
					(Reco_mu_nTrkWMea_np[iMumi_np] > 5) &&
					(Reco_mu_nPixWMea_np[iMumi_np] > 0) &&
					(TMath::Abs(Reco_mu_dxy_np[iMumi_np]) < 0.3) &&
					(TMath::Abs(Reco_mu_dz_np[iMumi_np]) < 20.) &&
					passMuonTypeMi_np
				);

				if (!(muplSoft_np && mumiSoft_np)) continue;

				double rap = TMath::Abs(y_np);
				double ctau = Reco_QQ_ctau3D_np[iQQ];

				// Fill forward rapidity histograms (1.6-2.4)
				if (rap >= 1.6 && rap < 2.4) {
					for (int ipt = 0; ipt < nPtBins_for; ipt++) {
						if (pt_np >= ptBin_for[ipt] && pt_np < ptBin_for[ipt+1]) {
							h_decay_np_pt_for[ipt]->Fill(ctau);
						}
					}
					if (pt_np > 3.5 && pt_np < 40.0) {
						h_decay_np_int_for->Fill(ctau);
					}
				}
				// Fill mid rapidity histograms (0-1.6)
				if (rap < 1.6) {
					for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
						if (pt_np >= ptBin_mid[ipt] && pt_np < ptBin_mid[ipt+1]) {
							h_decay_np_pt_mid[ipt]->Fill(ctau);
						}
					}
					if (pt_np > 6.5 && pt_np < 40.0) {
						h_decay_np_int_mid->Fill(ctau);
					}
				}
			}
		}

		// Calculate ctau3D cut for integrated bins
		cout << endl << "Integrated bins:" << endl;
		{
			// Forward rapidity integrated: 1.6<|y|<2.4, 3.5<pT<40
			double totalPR = h_decay_pr_int_for->GetEntries();
			double totalNP = h_decay_np_int_for->GetEntries();
			if (totalPR < 10) {
				cout << "  Forward integrated: Not enough PRMC entries (" << (int)totalPR << ")" << endl;
			} else {
				TH1D* h_deff_pr_int_for = new TH1D("h_deff_pr_int_for", ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
				TH1D* h_deff_np_int_for = new TH1D("h_deff_np_int_for", ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
				for (int bins = 0; bins < h_decay_pr_int_for->GetNbinsX(); bins++) {
					if (state == 1) {
						h_deff_pr_int_for->SetBinContent(bins, h_decay_pr_int_for->Integral(0, bins) / totalPR);
						h_deff_np_int_for->SetBinContent(bins, 1 - (h_decay_np_int_for->Integral(0, bins) / totalNP));

						if (h_decay_pr_int_for->Integral(0, bins) < (0.899 * totalPR) ||
							h_decay_pr_int_for->Integral(0, bins) > (totalPR * 0.9099)) continue;

						ctau3D_cut_int_for = h_decay_pr_int_for->GetBinCenter(bins);
						pr_eff_int_for = h_decay_pr_int_for->Integral(0, bins) / totalPR;
						np_res_int_for = (h_decay_np_int_for->Integral(0, bins) / totalNP);
					} else if (state == 2) {
						h_deff_np_int_for->SetBinContent(bins, h_decay_np_int_for->Integral(0, bins) / totalNP);
						h_deff_pr_int_for->SetBinContent(bins, 1 - (h_decay_pr_int_for->Integral(0, bins) / totalPR));

						if (h_decay_pr_int_for->Integral(0, bins) < (0.97 * totalPR) ||
							h_decay_pr_int_for->Integral(0, bins) > (totalPR * 0.98)) continue;
						ctau3D_cut_int_for = h_decay_pr_int_for->GetBinCenter(bins);
						np_res_int_for = 1 - (h_decay_np_int_for->Integral(0, bins) / totalNP);
						pr_eff_int_for = 1 - (h_decay_pr_int_for->Integral(0, bins) / totalPR);
					}
				}
				cout << Form("  Forward integrated (3.5-40): ctau cut = %.4f, PR eff = %.3f, NP res = %.3f",
					ctau3D_cut_int_for, pr_eff_int_for, np_res_int_for) << endl;
			}
		}

		{
			// Mid rapidity integrated: |y|<1.6, 6.5<pT<40
			double totalPR = h_decay_pr_int_mid->GetEntries();
			double totalNP = h_decay_np_int_mid->GetEntries();
			if (totalPR < 10) {
				cout << "  Mid integrated: Not enough PRMC entries (" << (int)totalPR << ")" << endl;
			} else {
				TH1D* h_deff_pr_int_mid = new TH1D("h_deff_pr_int_mid", ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
				TH1D* h_deff_np_int_mid = new TH1D("h_deff_np_int_mid", ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
				for (int bins = 0; bins < h_decay_pr_int_mid->GetNbinsX(); bins++) {
					if (state == 1) {
						h_deff_pr_int_mid->SetBinContent(bins, h_decay_pr_int_mid->Integral(0, bins) / totalPR);
						h_deff_np_int_mid->SetBinContent(bins, 1 - (h_decay_np_int_mid->Integral(0, bins) / totalNP));

						if (h_decay_pr_int_mid->Integral(0, bins) < (0.899 * totalPR) ||
							h_decay_pr_int_mid->Integral(0, bins) > (totalPR * 0.9099)) continue;

						ctau3D_cut_int_mid = h_decay_pr_int_mid->GetBinCenter(bins);
						pr_eff_int_mid = h_decay_pr_int_mid->Integral(0, bins) / totalPR;
						np_res_int_mid = (h_decay_np_int_mid->Integral(0, bins) / totalNP);
					} else if (state == 2) {
						h_deff_np_int_mid->SetBinContent(bins, h_decay_np_int_mid->Integral(0, bins) / totalNP);
						h_deff_pr_int_mid->SetBinContent(bins, 1 - (h_decay_pr_int_mid->Integral(0, bins) / totalPR));

						if (h_decay_pr_int_mid->Integral(0, bins) < (0.97 * totalPR) ||
							h_decay_pr_int_mid->Integral(0, bins) > (totalPR * 0.98)) continue;
						ctau3D_cut_int_mid = h_decay_pr_int_mid->GetBinCenter(bins);
						np_res_int_mid = 1 - (h_decay_np_int_mid->Integral(0, bins) / totalNP);
						pr_eff_int_mid = 1 - (h_decay_pr_int_mid->Integral(0, bins) / totalPR);
					}
				}
				cout << Form("  Mid integrated (6.5-40): ctau cut = %.4f, PR eff = %.3f, NP res = %.3f",
					ctau3D_cut_int_mid, pr_eff_int_mid, np_res_int_mid) << endl;
			}
		}

		// Calculate ctau3D cut for each forward pt bin
		cout << endl << "Forward rapidity (1.6-2.4):" << endl;
		for (int ipt = 0; ipt < nPtBins_for; ipt++) {
			double totalPR = h_decay_pr_pt_for[ipt]->GetEntries();
			double totalNP = h_decay_np_pt_for[ipt]->GetEntries();
			
			if (totalPR < 10) {
				cout << Form("pt bin %.1f-%.1f: Not enough PRMC entries (%d)", ptBin_for[ipt], ptBin_for[ipt+1], (int)totalPR) << endl;
				continue;
			}

			// Create efficiency histograms
			TH1D* h_deff_pr_for = new TH1D(Form("h_deff_pr_for_%d", ipt), ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
			TH1D* h_deff_np_for = new TH1D(Form("h_deff_np_for_%d", ipt), ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
			TLine *lcutv_for = nullptr;
			TLine *lcuth_for = nullptr;
			TLine *lresi_for = nullptr;

			for (int bins = 0; bins < h_decay_pr_pt_for[ipt]->GetNbinsX(); bins++) {
				if (state == 1) {
					h_deff_pr_for->SetBinContent(bins, h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR);
					h_deff_np_for->SetBinContent(bins, 1 - (h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP));

					if (h_decay_pr_pt_for[ipt]->Integral(0, bins) < (0.899 * totalPR) || 
						h_decay_pr_pt_for[ipt]->Integral(0, bins) > (totalPR * 0.9099)) continue;
					
					ctau3D_cut_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->GetBinCenter(bins);
					pr_eff_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR;
					np_res_ptBin_for[ipt] = (h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP);
					
					lcutv_for = new TLine(h_decay_pr_pt_for[ipt]->GetBinCenter(bins), 0, h_decay_pr_pt_for[ipt]->GetBinCenter(bins), 1);
					lcuth_for = new TLine(xmin, pr_eff_ptBin_for[ipt], h_decay_pr_pt_for[ipt]->GetBinCenter(bins), pr_eff_ptBin_for[ipt]);
					lresi_for = new TLine(xmin, np_res_ptBin_for[ipt], h_decay_pr_pt_for[ipt]->GetBinCenter(bins), np_res_ptBin_for[ipt]);
				} else if (state == 2) {
					h_deff_np_for->SetBinContent(bins, h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP);
					h_deff_pr_for->SetBinContent(bins, 1 - (h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR));

					if (h_decay_pr_pt_for[ipt]->Integral(0, bins) < (0.97 * totalPR) || 
						h_decay_pr_pt_for[ipt]->Integral(0, bins) > (totalPR * 0.98)) continue;
					
					ctau3D_cut_ptBin_for[ipt] = h_decay_pr_pt_for[ipt]->GetBinCenter(bins);
					np_res_ptBin_for[ipt] = 1 - (h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP);
					pr_eff_ptBin_for[ipt] = 1 - (h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR);
					
					lcutv_for = new TLine(h_decay_pr_pt_for[ipt]->GetBinCenter(bins), 0, h_decay_pr_pt_for[ipt]->GetBinCenter(bins), 1);
					lcuth_for = new TLine(xmin, h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP, h_decay_pr_pt_for[ipt]->GetBinCenter(bins), h_decay_np_pt_for[ipt]->Integral(0, bins) / totalNP);
					lresi_for = new TLine(xmin, h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR, h_decay_pr_pt_for[ipt]->GetBinCenter(bins), h_decay_pr_pt_for[ipt]->Integral(0, bins) / totalPR);
				}
			}
			
			cout << Form("  pt %.1f-%.1f: ctau cut = %.4f, PR eff = %.3f, NP res = %.3f", 
				ptBin_for[ipt], ptBin_for[ipt+1], ctau3D_cut_ptBin_for[ipt], 
				pr_eff_ptBin_for[ipt], np_res_ptBin_for[ipt]) << endl;
			
			hpt_ctau_1->SetBinContent(ipt+1, ctau3D_cut_ptBin_for[ipt]);

			// Draw and save plots for this pt bin
			TCanvas* c_decayL_for = new TCanvas(Form("c_decayL_for_%d", ipt), "", 600, 600);
			h_deff_pr_for->Draw("l");
			h_deff_np_for->Draw("same");
			h_deff_np_for->SetLineColor(kRed+2);
			jumSun(xmin, 1, xmax, 1, 1, 1);
			if (lcutv_for) lcutv_for->Draw("same");
			
			drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptBin_for[ipt], ptBin_for[ipt+1]), pos_x_mass, pos_y, text_color, text_size);
			drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", 1.6, 2.4), pos_x_mass, pos_y - pos_y_diff, text_color, text_size);
			drawText("|#eta^{#mu}| < 2.4", pos_x_mass, pos_y - pos_y_diff * 1.5, text_color, text_size);
			
			if (state == 1) {
				drawText(Form("L_{J/psi} cut: %.4f", ctau3D_cut_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
				drawText(Form("PR J/psi eff: %.3f", pr_eff_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 3.5, text_color, text_size);
				drawText(Form("NP J/psi res: %.3f", np_res_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 4, kRed + 2, text_size);
				c_decayL_for->Update();
				c_decayL_for->SaveAs(Form("figs_1S_pp/decayL/PRMC/decay_state%d_pt%.1f-%.1f_y1.6-2.4_m%.1f-%.1f.pdf", state, ptBin_for[ipt], ptBin_for[ipt+1], massLow, massHigh));
			} else if (state == 2) {
				drawText(Form("L_{J/psi} cut: %.4f", ctau3D_cut_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
				drawText(Form("NP J/psi eff: %.3f", np_res_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 3.5, kRed + 2, text_size);
				drawText(Form("PR J/psi res: %.3f", pr_eff_ptBin_for[ipt]), pos_x_mass, pos_y - pos_y_diff * 4, text_color, text_size);
				c_decayL_for->Update();
				c_decayL_for->SaveAs(Form("figs_1S_pp/decayL/NPMC/decay_state%d_pt%.1f-%.1f_y1.6-2.4_m%.1f-%.1f.pdf", state, ptBin_for[ipt], ptBin_for[ipt+1], massLow, massHigh));
			}
			delete c_decayL_for;
		}

		// Calculate ctau3D cut for each mid pt bin
		cout << endl << "Mid rapidity (0-1.6):" << endl;
		for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
			double totalPR = h_decay_pr_pt_mid[ipt]->GetEntries();
			double totalNP = h_decay_np_pt_mid[ipt]->GetEntries();
			
			if (totalPR < 10) {
				cout << Form("pt bin %.1f-%.1f: Not enough PRMC entries (%d)", ptBin_mid[ipt], ptBin_mid[ipt+1], (int)totalPR) << endl;
				continue;
			}

			// Create efficiency histograms
			TH1D* h_deff_pr_mid = new TH1D(Form("h_deff_pr_mid_%d", ipt), ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
			TH1D* h_deff_np_mid = new TH1D(Form("h_deff_np_mid_%d", ipt), ";l_{J/#psi};Efficiency", nbins, xmin, xmax);
			TLine *lcutv_mid = nullptr;
			TLine *lcuth_mid = nullptr;
			TLine *lresi_mid = nullptr;

			for (int bins = 0; bins < h_decay_pr_pt_mid[ipt]->GetNbinsX(); bins++) {
				if (state == 1) {
					h_deff_pr_mid->SetBinContent(bins, h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR);
					h_deff_np_mid->SetBinContent(bins, 1 - (h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP));

					if (h_decay_pr_pt_mid[ipt]->Integral(0, bins) < (0.899 * totalPR) || 
						h_decay_pr_pt_mid[ipt]->Integral(0, bins) > (totalPR * 0.9099)) continue;
					
					ctau3D_cut_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->GetBinCenter(bins);
					pr_eff_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR;
					np_res_ptBin_mid[ipt] = (h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP);
					
					lcutv_mid = new TLine(h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), 0, h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), 1);
					lcuth_mid = new TLine(xmin, pr_eff_ptBin_mid[ipt], h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), pr_eff_ptBin_mid[ipt]);
					lresi_mid = new TLine(xmin, np_res_ptBin_mid[ipt], h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), np_res_ptBin_mid[ipt]);
				} else if (state == 2) {
					h_deff_np_mid->SetBinContent(bins, h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP);
					h_deff_pr_mid->SetBinContent(bins, 1 - (h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR));

					if (h_decay_pr_pt_mid[ipt]->Integral(0, bins) < (0.97 * totalPR) || 
						h_decay_pr_pt_mid[ipt]->Integral(0, bins) > (totalPR * 0.98)) continue;
					
					ctau3D_cut_ptBin_mid[ipt] = h_decay_pr_pt_mid[ipt]->GetBinCenter(bins);
					np_res_ptBin_mid[ipt] = 1 - (h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP);
					pr_eff_ptBin_mid[ipt] = 1 - (h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR);
					
					lcutv_mid = new TLine(h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), 0, h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), 1);
					lcuth_mid = new TLine(xmin, h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP, h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), h_decay_np_pt_mid[ipt]->Integral(0, bins) / totalNP);
					lresi_mid = new TLine(xmin, h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR, h_decay_pr_pt_mid[ipt]->GetBinCenter(bins), h_decay_pr_pt_mid[ipt]->Integral(0, bins) / totalPR);
				}
			}
			
			cout << Form("  pt %.1f-%.1f: ctau cut = %.4f, PR eff = %.3f, NP res = %.3f", 
				ptBin_mid[ipt], ptBin_mid[ipt+1], ctau3D_cut_ptBin_mid[ipt], 
				pr_eff_ptBin_mid[ipt], np_res_ptBin_mid[ipt]) << endl;
			
			hpt_ctau_2->SetBinContent(ipt+1, ctau3D_cut_ptBin_mid[ipt]);

			// Draw and save plots for this pt bin
			TCanvas* c_decayL_mid = new TCanvas(Form("c_decayL_mid_%d", ipt), "", 600, 600);
			h_deff_pr_mid->Draw("l");
			h_deff_np_mid->Draw("same");
			h_deff_np_mid->SetLineColor(kRed+2);
			jumSun(xmin, 1, xmax, 1, 1, 1);
			if (lcutv_mid) lcutv_mid->Draw("same");
			
			drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptBin_mid[ipt], ptBin_mid[ipt+1]), pos_x_mass, pos_y, text_color, text_size);
			drawText(Form("|y^{#mu#mu}| < %.1f", 1.6), pos_x_mass, pos_y - pos_y_diff, text_color, text_size);
			drawText("|#eta^{#mu}| < 2.4", pos_x_mass, pos_y - pos_y_diff * 1.5, text_color, text_size);
			
			if (state == 1) {
				drawText(Form("L_{J/psi} cut: %.4f", ctau3D_cut_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
				drawText(Form("PR J/psi eff: %.3f", pr_eff_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 3.5, text_color, text_size);
				drawText(Form("NP J/psi res: %.3f", np_res_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 4, kRed + 2, text_size);
				c_decayL_mid->Update();
				c_decayL_mid->SaveAs(Form("figs_1S_pp/decayL/PRMC/decay_state%d_pt%.1f-%.1f_y0.0-1.6_m%.1f-%.1f.pdf", state, ptBin_mid[ipt], ptBin_mid[ipt+1], massLow, massHigh));
			} else if (state == 2) {
				drawText(Form("L_{J/psi} cut: %.4f", ctau3D_cut_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 3, text_color, text_size);
				drawText(Form("NP J/psi eff: %.3f", np_res_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 3.5, kRed + 2, text_size);
				drawText(Form("PR J/psi res: %.3f", pr_eff_ptBin_mid[ipt]), pos_x_mass, pos_y - pos_y_diff * 4, text_color, text_size);
				c_decayL_mid->Update();
				c_decayL_mid->SaveAs(Form("figs_1S_pp/decayL/NPMC/decay_state%d_pt%.1f-%.1f_y0.0-1.6_m%.1f-%.1f.pdf", state, ptBin_mid[ipt], ptBin_mid[ipt+1], massLow, massHigh));
			}
			delete c_decayL_mid;
		}

		// Save pt-binned results
		TString rootFileNamePt;
		if (state == 1) {
			rootFileNamePt = Form("./roots_1S_pp/ctau3D_cut_ptBin_PRMC_y%.1f-%.1f.root", yLow, yHigh);
		} else {
			rootFileNamePt = Form("./roots_1S_pp/ctau3D_cut_ptBin_NPMC_y%.1f-%.1f.root", yLow, yHigh);
		}
		TFile *wfPt = new TFile(rootFileNamePt, "RECREATE");
		wfPt->cd();
		
		hpt_ctau_1->Write();
		hpt_ctau_2->Write();

		// Save integrated-bin results
		auto h_ctau_int_for = new TH1D("h_ctau_int_for", "ctau cut (forward integrated)", 1, 0, 1);
		auto h_ctau_int_mid = new TH1D("h_ctau_int_mid", "ctau cut (mid integrated)", 1, 0, 1);
		auto h_pr_eff_int_for = new TH1D("h_pr_eff_int_for", "PR eff (forward integrated)", 1, 0, 1);
		auto h_pr_eff_int_mid = new TH1D("h_pr_eff_int_mid", "PR eff (mid integrated)", 1, 0, 1);
		auto h_np_res_int_for = new TH1D("h_np_res_int_for", "NP res (forward integrated)", 1, 0, 1);
		auto h_np_res_int_mid = new TH1D("h_np_res_int_mid", "NP res (mid integrated)", 1, 0, 1);

		h_ctau_int_for->SetBinContent(1, ctau3D_cut_int_for);
		h_ctau_int_mid->SetBinContent(1, ctau3D_cut_int_mid);
		h_pr_eff_int_for->SetBinContent(1, pr_eff_int_for);
		h_pr_eff_int_mid->SetBinContent(1, pr_eff_int_mid);
		h_np_res_int_for->SetBinContent(1, np_res_int_for);
		h_np_res_int_mid->SetBinContent(1, np_res_int_mid);

		h_ctau_int_for->Write();
		h_ctau_int_mid->Write();
		h_pr_eff_int_for->Write();
		h_pr_eff_int_mid->Write();
		h_np_res_int_for->Write();
		h_np_res_int_mid->Write();

		h_decay_pr_int_for->Write();
		h_decay_np_int_for->Write();
		h_decay_pr_int_mid->Write();
		h_decay_np_int_mid->Write();
		
		// Save histograms for each pt bin
		for (int ipt = 0; ipt < nPtBins_for; ipt++) {
			h_decay_pr_pt_for[ipt]->Write();
			h_decay_np_pt_for[ipt]->Write();
		}
		for (int ipt = 0; ipt < nPtBins_mid; ipt++) {
			h_decay_pr_pt_mid[ipt]->Write();
			h_decay_np_pt_mid[ipt]->Write();
		}
		
		wfPt->Close();
		cout << endl << "Saved pt-binned ROOT file: " << rootFileNamePt << endl;
	}

	cout << endl << "============================================" << endl;
	cout << "All decay length analysis completed!" << endl;
	cout << "============================================" << endl;
}
