#include "../../Eff_Acc/TreeSetting.h"
#include "../../cutsAndBin.h"

std::pair<string, string> get_histnames(bool get_pp, bool isPR, bool isFwd, float ptLow, float ptHigh){
	string name_ref ="";
	string name_cmp="";
	double yl, yh;
	if( isFwd){ yl=1.6; yh=2.4;}
	if(!isFwd){ yl=0.0; yh=1.6;}
	if(get_pp) name_ref = Form("../../Eff_Acc/roots/TnP/mc_eff_pp_%s_pt%.1f_%.1f_y%.1f_%.1f_accYes_mass3.4_4.1_isTnP1_isPtWeight1_NOM.root", (isPR)? "PR" : "NP", ptLow, ptHigh, yl, yh);
	if(get_pp) name_cmp = Form("../../Eff_Acc/roots/TnP/mc_eff_pp_%s_pt%.1f_%.1f_y%.1f_%.1f_accYes_mass3.4_4.1_isTnP1_isPtWeight1_SYSTNP.root", (isPR)? "PR" : "NP", ptLow, ptHigh, yl, yh);
	return std::make_pair(name_ref, name_cmp);
};

void TNPUnc_pp(bool incl_pp, bool incl_pbpb, bool dopt){
	string tag_kine = (dopt) ? "pt" : "cent"; 
	string tag_coll = (incl_pp) ? (incl_pbpb) ? "" : "_ppOnly" : "_PbPbOnly";
	TFile* output = new TFile(Form("./syst_roots/syst_%s_TNP%s_pp.root", tag_kine.c_str(), tag_coll.c_str()),"recreate");

	const size_t numbin_for =  (dopt) ? 5 : 1;
	const size_t numbin_mid =  (dopt) ? 7 : 1;
	double ptBin_for[6] = {0,3.5,6.5,9,12,50};
	double ptBin_mid[8] = {0,6.5,9,12,15,20,25,50};
	double centBin[2] = {0,180};
	float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4}; // Not Used *Yet*

	string name_hist_for_PR = "fwd_PR";
	string name_hist_mid_PR = "mid_PR";
	string name_hist_for_NP = "fwd_NP";
	string name_hist_mid_NP = "mid_NP";
	TH1D* hist_for_PR = new TH1D(name_hist_for_PR.c_str(), "", numbin_for, (dopt) ? ptBin_for : centBin);
	TH1D* hist_mid_PR = new TH1D(name_hist_mid_PR.c_str(), "", numbin_mid, (dopt) ? ptBin_mid : centBin);
	TH1D* hist_for_NP = new TH1D(name_hist_for_NP.c_str(), "", numbin_for, (dopt) ? ptBin_for : centBin);
	TH1D* hist_mid_NP = new TH1D(name_hist_mid_NP.c_str(), "", numbin_mid, (dopt) ? ptBin_mid : centBin);


	auto getUncBin = [&](bool isPP, bool isPR, bool isFwd, bool isPT, float ptLow, float ptHigh){
		auto name_pair = get_histnames( isPP, isPR, isFwd, ptLow, ptHigh);
		TFile *f_hreco_ref, *f_hreco;
		string name_reco_ref = "hreco_tnpNom";
		std::cout << name_pair.first.c_str() << std::endl;
		std::cout << name_pair.second.c_str() << std::endl;
		f_hreco_ref = TFile::Open(Form("%s", name_pair.first.c_str() ) );
		f_hreco = TFile::Open(Form("%s", name_pair.second.c_str() ) );

		TH1D* hreco_ref;
		std::map<TString, TH1D*> map_hreco_tnp;
		double bincenter = (isPT) ? (double) (ptLow + ptHigh) /2 : (double) (0 + 180) /2 ;
		hreco_ref = (TH1D*) f_hreco_ref->Get(name_reco_ref.c_str());
		double COUNT_REF;
		//COUNT_REF = hreco_ref->GetBinContent(hreco_ref->FindBin(bincenter));
		COUNT_REF = hreco_ref->GetBinContent(1);
		double COUNT_SYS;
		double ERR_CVAR = 0;
		TString the_string = "tnp_sys";
		TH1D* clone_hreco = (TH1D*) ((TH1D*) f_hreco->Get(Form("hreco_%s",the_string.Data())))->Clone();	
		clone_hreco->SetDirectory(0);

		map_hreco_tnp[the_string] = clone_hreco;

		COUNT_SYS = map_hreco_tnp[the_string]->GetSum();
		std::cout <<"NAME : " << the_string.Data() << std::endl;

		std::cout <<"COUNT_REF " << COUNT_REF << ", COUNT_SYS " << COUNT_SYS - COUNT_REF << std::endl;
		ERR_CVAR = TMath::Sqrt( TMath::Power(COUNT_SYS - COUNT_REF, 2) );
		std::cout << "ER_CVAR " << ERR_CVAR << std::endl;

		f_hreco->Close();
		f_hreco_ref->Close();
		double COUNT_TNP = TMath::Sqrt( TMath::Power(ERR_CVAR, 2) ); 

		std::cout << Form("Efficiency for ref, TNP : %.5f, %.5f", COUNT_REF, COUNT_TNP ) << std::endl;

		double tnpUnc = (COUNT_TNP )/ COUNT_REF;
		//double tnpUnc = COUNT_REF / (COUNT_TNP);
		std::cout << "\033[1;31m" << "tnpUnc : " << tnpUnc << "\033[0m" << std::endl;
		return tnpUnc;
	};

	//(bool isPP, bool isPR, bool isFwd, bool isPT, float ptLow, float ptHigh, int cbinLow, int cbinHigh)
	if( dopt ){
		for( auto idx : ROOT::TSeqI(numbin_for)){
			hist_for_PR->SetBinContent(idx+1, getUncBin(incl_pp, true , true, dopt, (float) ptBin_for[idx], (float) ptBin_for[idx+1])) ;
			hist_for_NP->SetBinContent(idx+1, getUncBin(incl_pp, false, true, dopt, (float) ptBin_for[idx], (float) ptBin_for[idx+1])) ;
		}
		for( auto idx : ROOT::TSeqI(numbin_mid)){
			hist_mid_PR->SetBinContent(idx+1, getUncBin(incl_pp, true , false, dopt, (float) ptBin_mid[idx], (float) ptBin_mid[idx+1])) ;
			hist_mid_NP->SetBinContent(idx+1, getUncBin(incl_pp, false, false, dopt, (float) ptBin_mid[idx], (float) ptBin_mid[idx+1])) ;
		}
	}
	if( !dopt ){
		for( auto idx : ROOT::TSeqI(numbin_for)){
			hist_for_PR->SetBinContent(idx+1, getUncBin(incl_pp, true , true, dopt, 3.5, 50 )) ;
			hist_for_NP->SetBinContent(idx+1, getUncBin(incl_pp, false, true, dopt, 3.5, 50 )) ;
		}
		for( auto idx : ROOT::TSeqI(numbin_mid)){
			hist_mid_PR->SetBinContent(idx+1, getUncBin(incl_pp, true , false, dopt, 6.5, 50)) ;
			hist_mid_NP->SetBinContent(idx+1, getUncBin(incl_pp, false, false, dopt, 6.5, 50)) ;
		}
	}

	output->cd();
	hist_for_PR->Write();
	hist_mid_PR->Write();
	hist_for_NP->Write();
	hist_mid_NP->Write();
	output->Close();
}

