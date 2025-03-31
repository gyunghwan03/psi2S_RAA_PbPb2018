#include "../../Eff_Acc/TreeSetting.h"
#include "../../cutsAndBin.h"

std::pair<string, string> get_histnames(bool get_pbpb, bool isPR, bool isFwd, float ptLow, float ptHigh, int cbinLow, int cbinHigh ){
	string name_ref ="";
	if( isPR && get_pbpb ) name_ref = "../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtW1_tnp1_20250202.root";
	if(!isPR && get_pbpb) name_ref = "../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtW1_tnp1_20250202.root";
	string name_cmp="";
	double yl, yh;
	if( isFwd){ yl=1.6; yh=2.4;}
	if(!isFwd){ yl=0.0; yh=1.6;}
	if(get_pbpb) name_cmp = Form("../../Eff_Acc/roots/TnP_Jpsi/mc_eff_PbPb_%s_pt%.1f_%.1f_y%.1f_%.1f_accYes_mass2.6_3.5_cent%d_%d_isTnP1_isPtWeight1_SYSTNP.root", (isPR)? "PR" : "NP", ptLow, ptHigh, yl, yh, cbinLow, cbinHigh);
	return std::make_pair(name_ref, name_cmp);
};

void TNPUnc(bool incl_pp, bool incl_pbpb, bool dopt){
	string tag_kine = (dopt) ? "pt" : "cent"; 
	string tag_coll = (incl_pp) ? (incl_pbpb) ? "" : "_ppOnly" : "_PbPbOnly";
	TFile* output = new TFile(Form("./syst_roots/syst_%s_TNP%s_test.root", tag_kine.c_str(), tag_coll.c_str()),"recreate");

	const size_t numbin_for =  (dopt) ? 5 : 4;
	const size_t numbin_mid =  (dopt) ? 7 : 6;
	double ptBin_for[6] = {0,3.5,6.5,9,12,40};
	double ptBin_mid[8] = {0,6.5,9,12,15,20,25,40};
	double centBin_for[5] = {0,10,30,50,90};
	double centBin_mid[7] = {0,10,20,30,40,50,90};
	float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4}; // Not Used *Yet*

	string name_hist_for_PR = "fwd_PR";
	string name_hist_mid_PR = "mid_PR";
	string name_hist_for_NP = "fwd_NP";
	string name_hist_mid_NP = "mid_NP";
	TH1D* hist_for_PR = new TH1D(name_hist_for_PR.c_str(), "", numbin_for, (dopt) ? ptBin_for : centBin_for);
	TH1D* hist_mid_PR = new TH1D(name_hist_mid_PR.c_str(), "", numbin_mid, (dopt) ? ptBin_mid : centBin_mid);
	TH1D* hist_for_NP = new TH1D(name_hist_for_NP.c_str(), "", numbin_for, (dopt) ? ptBin_for : centBin_for);
	TH1D* hist_mid_NP = new TH1D(name_hist_mid_NP.c_str(), "", numbin_mid, (dopt) ? ptBin_mid : centBin_mid);

	auto get_input_config = [&] (int cvr, int mod)
	{
		TString _input_config, _idstr, _trkstr, _trgstr;
		if( cvr == 1 ){ _trkstr = "Nom"; _trgstr = "Nom"; _idstr = Form("%s",mode_str[mod-1].c_str()); }
		if( cvr == 2 ){ _idstr = "Nom"; _trgstr = "Nom"; _trkstr = Form("%s",mode_str[mod-1].c_str()); }
		if( cvr == 3 ){ _idstr = "Nom"; _trkstr = "Nom"; _trgstr = Form("%s",mode_str[mod-1].c_str()); }
		_input_config = Form("id%s_trk%s_trg%s",_idstr.Data(), _trkstr.Data(), _trgstr.Data()); 
		return _input_config;
	};

	auto getUncBin = [&](bool isPbPb, bool isPR, bool isFwd, bool isPT, float ptLow, float ptHigh, int cbinLow, int cbinHigh){
		auto name_pair = get_histnames( isPbPb, isPR, isFwd, ptLow, ptHigh, cbinLow, cbinHigh);
		TFile *f_hreco_ref, *f_hreco;
		string name_reco_ref = "";
		if( isFwd && isPT) name_reco_ref = "hpt_reco_1";
		if(!isFwd && isPT) name_reco_ref = "hpt_reco_2";
		if( isFwd &&!isPT) name_reco_ref = "hcent_reco_1";
		if(!isFwd &&!isPT) name_reco_ref = "hcent_reco_2";
		std::cout << name_pair.first.c_str() << std::endl;
		std::cout << name_pair.second.c_str() << std::endl;
		f_hreco_ref = TFile::Open(Form("%s", name_pair.first.c_str() ) );
		f_hreco = TFile::Open(Form("%s", name_pair.second.c_str() ) );

		TH1D* hreco_ref;
		std::map<TString, TH1D*> map_hreco_tnp;
		double bincenter = (isPT) ? (double) (ptLow + ptHigh) /2 : (double) (cbinLow + cbinHigh) /2;
		hreco_ref = (TH1D*) f_hreco_ref->Get(name_reco_ref.c_str());
		double COUNT_REF = hreco_ref->GetBinContent(hreco_ref->FindBin(bincenter));
		double ERR_CVAR[3] = {0,0,0};
		for( int cv : {1,2,3}){
			double COUNT_TAG, COUNT_STAT, COUNT_SYS;
			for( int md : {1,2,3,4,5}){
				TString the_string = get_input_config(cv, md);
				TH1D* clone_hreco = (TH1D*) ((TH1D*) f_hreco->Get(Form("hreco_%s",the_string.Data())))->Clone();	
				clone_hreco->SetDirectory(0);

				map_hreco_tnp[the_string] = clone_hreco;

				if(md == 1) COUNT_TAG = map_hreco_tnp[the_string]->GetSum();
			}
			std::vector<std::pair<int,int> > udp = {{2,3}, {4,5}};
			for( auto mdp : udp){
				double COUNT_UP, COUNT_DOWN;
				TString cfg_up, cfg_down;
				cfg_up = get_input_config(cv, mdp.first);
				cfg_down = get_input_config(cv, mdp.second);
				COUNT_UP = map_hreco_tnp[cfg_up]->GetSum();
				COUNT_DOWN = map_hreco_tnp[cfg_down]->GetSum();
				double _this = ( (fabs(COUNT_UP - COUNT_REF) - fabs(COUNT_DOWN - COUNT_REF)) >0) ? COUNT_UP : COUNT_DOWN;
				if( mdp.first == 2) COUNT_STAT = _this;
				if( mdp.first == 4) COUNT_SYS = _this;
			}
			std::cout <<"COUNT_REF " << COUNT_REF <<  ", COUNT_TAG " << COUNT_TAG - COUNT_REF << ", COUNT_STAT " << COUNT_STAT - COUNT_REF << ", COUNT_SYS " << COUNT_SYS - COUNT_REF << std::endl;
			ERR_CVAR[cv-1] = TMath::Sqrt( TMath::Power(COUNT_TAG -COUNT_REF, 2) + TMath::Power(COUNT_STAT- COUNT_REF, 2) + TMath::Power(COUNT_SYS - COUNT_REF, 2) );
			std::cout << "ER_CVAR[cv-1] " << ERR_CVAR[cv-1] << std::endl;
		}
		f_hreco->Close();
		f_hreco_ref->Close();
		double COUNT_TNP = TMath::Sqrt( TMath::Power(ERR_CVAR[0], 2) + TMath::Power(ERR_CVAR[1], 2) + TMath::Power(ERR_CVAR[2], 2) ); 

		std::cout << Form("Efficiency for ref, TNP : %.5f, %.5f", COUNT_REF, COUNT_TNP ) << std::endl;

		//double tnpUnc = (COUNT_TNP )/ COUNT_REF;
		double tnpUnc = COUNT_REF/ (COUNT_TNP);
		std::cout << "\033[1;31m" << "tnpUnc : " << tnpUnc << "\033[0m" << std::endl;
		return tnpUnc;
	};

	//(bool isPbPb, bool isPR, bool isFwd, bool isPT, float ptLow, float ptHigh, int cbinLow, int cbinHigh)
	if( dopt ){
		for( auto idx : ROOT::TSeqI(numbin_for)){
			hist_for_PR->SetBinContent(idx+1, getUncBin(incl_pbpb, true , true, dopt, (float) ptBin_for[idx], (float) ptBin_for[idx+1] , 0, 180)) ;
			hist_for_NP->SetBinContent(idx+1, getUncBin(incl_pbpb, false, true, dopt, (float) ptBin_for[idx], (float) ptBin_for[idx+1] , 0, 180)) ;
		}
		for( auto idx : ROOT::TSeqI(numbin_mid)){
			hist_mid_PR->SetBinContent(idx+1, getUncBin(incl_pbpb, true , false, dopt, (float) ptBin_mid[idx], (float) ptBin_mid[idx+1] , 0, 180)) ;
			hist_mid_NP->SetBinContent(idx+1, getUncBin(incl_pbpb, false, false, dopt, (float) ptBin_mid[idx], (float) ptBin_mid[idx+1] , 0, 180)) ;
		}
	}
	if( !dopt ){
		for( auto idx : ROOT::TSeqI(numbin_for)){
			hist_for_PR->SetBinContent(idx+1, getUncBin(incl_pbpb, true , true, dopt, 3.5, 40 , (int) centBin_for[idx]*2, (int) centBin_for[idx+1]*2 )) ;
			hist_for_NP->SetBinContent(idx+1, getUncBin(incl_pbpb, false, true, dopt, 3.5, 40 , (int) centBin_for[idx]*2, (int) centBin_for[idx+1]*2 )) ;
		}
		for( auto idx : ROOT::TSeqI(numbin_mid)){
			hist_mid_PR->SetBinContent(idx+1, getUncBin(incl_pbpb, true , false, dopt, 6.5, 40 , (int) centBin_mid[idx]*2, (int) centBin_mid[idx+1]*2 )) ;
			hist_mid_NP->SetBinContent(idx+1, getUncBin(incl_pbpb, false, false, dopt, 6.5, 40 , (int) centBin_mid[idx]*2, (int) centBin_mid[idx+1]*2 )) ;
		}
	}

	output->cd();
	hist_for_PR->Write();
	hist_mid_PR->Write();
	hist_for_NP->Write();
	hist_mid_NP->Write();
	output->Close();
}

