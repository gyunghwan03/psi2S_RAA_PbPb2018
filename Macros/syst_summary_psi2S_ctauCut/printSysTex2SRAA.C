#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TMath.h>
#include <TChain.h>
#include <fstream>
#include <iostream>

void printSysTex2SRAA()
{
	using namespace std;

	/*
		 TString Acc = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Acceptance/Acceptance_sys_v3.root";
		 TString bkgPdf = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/BackgroundPDFVariation/AllParmFreeFit/BackgroundPDFVar_sys_v3.root";
		 TString bkgv2 = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Background_v2Variation/SigAllFreeFitFix/V2BkgFuncVar_sys_v3.root";
		 TString EvntSel = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/EventSelection/EventSelection_sys_v3.root";
		 TString SigPdf = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalPDFVariation/AllParFreeFit/SignalPDFVar_sys_v3.root";
		 TString SigPar = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalParVariation/SignalParVar_sys_v3.root";
		 TString TnP = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/TnP/TnP_sys_v3.root";
		 */

	const int nSyst = 8;
	const int nCent = 2;
	const char* RapName[2] ={ "fwd","mid" };
	const char* RapTxt[2] = {"1.6 < |y| <2.4","|y| < 1.6"};

	const char* fnamePt[nSyst]; // nCent==0 : pT & nCent==1 : Cent
	const char* fnameCent[nSyst]; // nCent==0 : pT & nCent==1 : Cent


	fnamePt[0] = "syst_roots/syst_pt_sigPDF.root";
	fnamePt[1] = "syst_roots/syst_pt_sigPAR.root";
	fnamePt[2] = "syst_roots/syst_pt_bkgPDF.root";
	fnamePt[3] = "syst_roots/syst_pt_HF.root";
//	fnamePt[4] = "syst_roots/syst_pt_bFrac.root";
	fnamePt[4] = "syst_roots/syst_pt_TNP.root";
	fnamePt[5] = "syst_roots/syst_pt_eff.root";
	fnamePt[6] = "syst_roots/syst_pt_acc.root";
	fnamePt[7] = "syst_roots/total_syst.root";

	fnameCent[0] = "syst_roots/syst_cent_sigPDF.root";
	fnameCent[1] = "syst_roots/syst_cent_sigPAR.root";
	fnameCent[2] = "syst_roots/syst_cent_bkgPDF.root";
	fnameCent[3] = "syst_roots/syst_cent_HF.root";
	//fnameCent[4] = "syst_roots/syst_cent_bFrac.root";
	fnameCent[4] = "syst_roots/syst_cent_TNP.root";
	fnameCent[5] = "syst_roots/syst_cent_eff.root";	
	fnameCent[6] = "syst_roots/syst_cent_acc.root";
	fnameCent[7] = "syst_roots/total_syst.root"; 


	const int nPtmid = 6;
	const int nPtfwd = 4;
	const int nCentBinmid = 6;
	const int nCentBinfwd = 4;
	const char* ptBinmidName[nPtmid] = { "6.5 -- 9", "9 -- 12", "12 -- 15", "15 -- 20", "20 -- 25",  "25 -- 40" };
	const char* ptBinfwdName[nPtfwd] = { "3.5 -- 6.5", "6.5 -- 9", "9 -- 12", "12 -- 40" };
	const char* CentBinmidName[nCentBinmid] = { "0 -- 10", "10 -- 20", "20 -- 30", "30 -- 40", "40 -- 50", "50 -- 90"};
	const char* CentBinfwdName[nCentBinfwd] = { "0 -- 10", "10 -- 30", "30 -- 50", "50 -- 90" };

	TString SystName[nSyst+1] = { "SigPDF", "sigPar", "bkgPDF", "EventSel", "TNP", "Eff", "Acc", "Total"};

	TFile *filePt[nSyst];
	TFile *fileCent[nSyst];
	TH1D *hpt_mid_PR[nSyst];
	TH1D *hpt_mid_NP[nSyst];
	TH1D *hpt_fwd_PR[nSyst];
	TH1D *hpt_fwd_NP[nSyst];
	TH1D *hCent_mid_PR[nSyst];
	TH1D *hCent_mid_NP[nSyst];
	TH1D *hCent_fwd_PR[nSyst];
	TH1D *hCent_fwd_NP[nSyst];
	TH1D *hpt_mid_PR_merged;
	TH1D *hpt_mid_NP_merged;
	TH1D *hpt_fwd_PR_merged;
	TH1D *hpt_fwd_NP_merged;
	TH1D *hCent_mid_PR_merged;
	TH1D *hCent_mid_NP_merged;
	TH1D *hCent_fwd_PR_merged;
	TH1D *hCent_fwd_NP_merged;

		for(int j=0; j<nSyst; j++){
		filePt[j] = new TFile(fnamePt[j]);
		fileCent[j] = new TFile(fnameCent[j]);
			//if(j==4) { //TnP
			//	continue; // skip TnP for now
			//}
			if(j<nSyst-2){
				hpt_mid_PR[j] = (TH1D*) filePt[j]->Get("mid_PR"); 
				hpt_mid_NP[j] = (TH1D*) filePt[j]->Get("mid_NP"); 
				hpt_fwd_PR[j] = (TH1D*) filePt[j]->Get("fwd_PR"); 
				hpt_fwd_NP[j] = (TH1D*) filePt[j]->Get("fwd_NP"); 
				hCent_mid_PR[j] = (TH1D*) fileCent[j]->Get("mid_PR"); 
				hCent_mid_NP[j] = (TH1D*) fileCent[j]->Get("mid_NP"); 
				hCent_fwd_PR[j] = (TH1D*) fileCent[j]->Get("fwd_PR"); 
				hCent_fwd_NP[j] = (TH1D*) fileCent[j]->Get("fwd_NP"); 
			}
			else if(j==nSyst-2){
				hpt_mid_PR[j] = (TH1D*) filePt[j]->Get("mid"); 
				hpt_mid_NP[j] = (TH1D*) filePt[j]->Get("mid"); 
				hpt_fwd_PR[j] = (TH1D*) filePt[j]->Get("fwd"); 
				hpt_fwd_NP[j] = (TH1D*) filePt[j]->Get("fwd"); 
				hCent_mid_PR[j] = (TH1D*) fileCent[j]->Get("mid"); 
				hCent_mid_NP[j] = (TH1D*) fileCent[j]->Get("mid"); 
				hCent_fwd_PR[j] = (TH1D*) fileCent[j]->Get("fwd"); 
				hCent_fwd_NP[j] = (TH1D*) fileCent[j]->Get("fwd"); 
			}
			else if(j==nSyst-1){ //total
				hpt_mid_PR[j] = (TH1D*) filePt[j]->Get("mid_pt_PR"); 
				hpt_mid_NP[j] = (TH1D*) filePt[j]->Get("mid_pt_NP"); 
				hpt_fwd_PR[j] = (TH1D*) filePt[j]->Get("fwd_pt_PR"); 
				hpt_fwd_NP[j] = (TH1D*) filePt[j]->Get("fwd_pt_NP"); 
				hCent_mid_PR[j] = (TH1D*) fileCent[j]->Get("mid_cent_PR"); 
				hCent_mid_NP[j] = (TH1D*) fileCent[j]->Get("mid_cent_NP"); 
				hCent_fwd_PR[j] = (TH1D*) fileCent[j]->Get("fwd_cent_PR"); 
				hCent_fwd_NP[j] = (TH1D*) fileCent[j]->Get("fwd_cent_NP"); 
			}
			else if (j==4) continue; // skip TnP for now

		}
		
			hpt_mid_PR[4] = (TH1D*) filePt[nSyst-1]->Get("TNP_mid_pt_PR");
			hpt_mid_NP[4] = (TH1D*) filePt[nSyst-1]->Get("TNP_mid_pt_NP");
			hpt_fwd_PR[4] = (TH1D*) filePt[nSyst-1]->Get("TNP_fwd_pt_PR");
			hpt_fwd_NP[4] = (TH1D*) filePt[nSyst-1]->Get("TNP_fwd_pt_NP");
		
			hCent_mid_PR[4] = (TH1D*) fileCent[nSyst-1]->Get("TNP_mid_cent_PR");
			hCent_mid_NP[4] = (TH1D*) fileCent[nSyst-1]->Get("TNP_mid_cent_NP");
			hCent_fwd_PR[4] = (TH1D*) fileCent[nSyst-1]->Get("TNP_fwd_cent_PR");
			hCent_fwd_NP[4] = (TH1D*) fileCent[nSyst-1]->Get("TNP_fwd_cent_NP");

	

 
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/////////////////// Fwd  Rapidity ////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
	cout << "\\multicolumn{2}{c}{$1.6 < |y| < 2.4$} \\\\ \\hline" << endl;
	cout << "\\multicolumn{2}{c}{Prompt ~\\PGyP{2S}} \\\\ \\hline \\hline" << endl;
	for(int ipt = 0; ipt < nPtfwd; ipt++){
		for(int i=0; i<nSyst; i++){
			if(ipt==0 && i==0) cout << Form("\\multirow{%d}{*}{%s}	& %s	& ",nPtfwd, "0--90\\%",ptBinfwdName[ipt]);
			else if(ipt!=0 &&  i==0) cout << Form("				& %s	& ",ptBinfwdName[ipt]);
			//else if(ipt!=0 && ipt<nPtfwd && i==0) cout << Form("				& %s	& ",CentBinName[ipt]);
			if(i<nSyst-1) {
				cout << Form("%.7f",hpt_fwd_PR[i]->GetBinContent(ipt+1)) << " & ";
			}
			else if(i==nSyst-1) {
				if(ipt<nPtfwd-1) cout << Form("%.7f",hpt_fwd_PR[i]->GetBinContent(ipt+1)) << " \\\\ "<< endl;
				if(ipt==nPtfwd-1) cout << Form("%.7f",hpt_fwd_PR[i]->GetBinContent(ipt+1)) << "\\\\ \\hline " << endl;
			}
		}
	}
	for (int iCent = 0; iCent < nCentBinfwd; iCent++){
		for(int i=0; i<nSyst; i++){
			if(iCent==0 && i==0) cout << "0--10\\% 	&  \\multirow{6}{*}{3.5--40} &";
			else if(iCent!=0 && i==0) cout << Form("  %s \\%s  & & ", CentBinfwdName[iCent], "%");
			if(i<nSyst-1) cout << Form("%.7f",hCent_fwd_PR[i]->GetBinContent(iCent+1)) << " & ";
			else if(i==nSyst-1) {
				if(iCent<nCentBinfwd-1) cout << Form("%.7f",hCent_fwd_PR[i]->GetBinContent(iCent+1)) << " \\\\ " << endl;
				if(iCent==nCentBinfwd-1) cout << Form("%.7f",hCent_fwd_PR[i]->GetBinContent(iCent+1)) << "\\\\ \\hline" << endl;
			}
		}
	}
	cout << "\\multicolumn{2}{c}{NonPrompt ~\\PGyP{2S}} \\\\ \\hline \\hline" << endl;
	for(int ipt = 0; ipt < nPtfwd; ipt++){
		for(int i=0; i<nSyst; i++){
			if(ipt==0 && i==0) cout << Form("\\multirow{%d}{*}{%s}	& %s	& ",nPtfwd, "0--90\\%",ptBinfwdName[ipt]);
			else if(ipt!=0 &&  i==0) cout << Form("				& %s	& ",ptBinfwdName[ipt]);
			//else if(ipt!=0 && ipt<nPtfwd && i==0) cout << Form("				& %s	& ",CentBinName[ipt]);
			if(i<nSyst-1) {
				cout << Form("%.7f",hpt_fwd_NP[i]->GetBinContent(ipt+1)) << " & ";
			}
			else if(i==nSyst-1) {
				if(ipt<nPtfwd-1) cout << Form("%.7f",hpt_fwd_NP[i]->GetBinContent(ipt+1)) << " \\\\ "<< endl;
				if(ipt==nPtfwd-1) cout << Form("%.7f",hpt_fwd_NP[i]->GetBinContent(ipt+1)) << "\\\\ \\hline " << endl;
			}
		}
	}
	for (int iCent = 0; iCent < nCentBinfwd; iCent++){
		for(int i=0; i<nSyst; i++){
			if(iCent==0 && i==0) cout << "0--10\\% 	&  \\multirow{6}{*}{3.5--40} &";
			else if(iCent!=0 && i==0) cout << Form("  %s \\%s  & & ", CentBinfwdName[iCent], "%");
			if(i<nSyst-1) cout << Form("%.7f",hCent_fwd_NP[i]->GetBinContent(iCent+1)) << " & ";
			else if(i==nSyst-1) {
				if(iCent<nCentBinfwd-1) cout << Form("%.7f",hCent_fwd_NP[i]->GetBinContent(iCent+1)) << " \\\\ " << endl;
				if(iCent==nCentBinfwd-1) cout << Form("%.7f",hCent_fwd_NP[i]->GetBinContent(iCent+1)) << "\\\\ \\hline" << endl;
			}
		}
	}
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////// Mid Rapidity ////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
	cout << "\\multicolumn{2}{c}{$|y| < 1.6$} \\\\ \\hline" << endl;
	cout << "\\multicolumn{2}{c}{Prompt ~\\PGyP{2S}} \\\\ \\hline \\hline" << endl;
	for(int ipt = 0; ipt < nPtmid; ipt++){
		for(int i=0; i<nSyst; i++){
			if(ipt==0 && i==0) cout << Form("\\multirow{%d}{*}{%s}	& %s	& ",nPtmid, "0--90\\%",ptBinmidName[ipt]);
			else if(ipt!=0 &&  i==0) cout << Form("				& %s	& ",ptBinmidName[ipt]);
			//else if(ipt!=0 && ipt<nPtmid && i==0) cout << Form("				& %s	& ",CentBinName[ipt]);
			if(i<nSyst-1) {
				cout << Form("%.7f",hpt_mid_PR[i]->GetBinContent(ipt+1)) << " & ";
			}
			else if(i==nSyst-1) {
				if(ipt<nPtmid-1) cout << Form("%.7f",hpt_mid_PR[i]->GetBinContent(ipt+1)) << " \\\\ "<< endl;
				if(ipt==nPtmid-1) cout << Form("%.7f",hpt_mid_PR[i]->GetBinContent(ipt+1)) << "\\\\ \\hline " << endl;
			}
		}
	}
	for (int iCent = 0; iCent < nCentBinmid; iCent++){
		for(int i=0; i<nSyst; i++){
			if(iCent==0 && i==0) cout << "0--10\\% 	&  \\multirow{6}{*}{6.5--40} &";
			else if(iCent!=0 && i==0) cout << Form("  %s \\%s  & & ", CentBinmidName[iCent], "%");
			if(i<nSyst-1) cout << Form("%.7f",hCent_mid_PR[i]->GetBinContent(iCent+1)) << " & ";
			else if(i==nSyst-1) {
				if(iCent<nCentBinmid-1) cout << Form("%.7f",hCent_mid_PR[i]->GetBinContent(iCent+1)) << " \\\\ " << endl;
				if(iCent==nCentBinmid-1) cout << Form("%.7f",hCent_mid_PR[i]->GetBinContent(iCent+1)) << "\\\\ \\hline" << endl;
			}
		}
	}
	cout << "\\multicolumn{2}{c}{NonPrompt ~\\PGyP{2S}} \\\\ \\hline \\hline" << endl;
	for(int ipt = 0; ipt < nPtmid; ipt++){
		for(int i=0; i<nSyst; i++){
			if(ipt==0 && i==0) cout << Form("\\multirow{%d}{*}{%s}	& %s	& ",nPtmid, "0--90\\%",ptBinmidName[ipt]);
			else if(ipt!=0 &&  i==0) cout << Form("				& %s	& ",ptBinmidName[ipt]);
			//else if(ipt!=0 && ipt<nPtmid && i==0) cout << Form("				& %s	& ",CentBinName[ipt]);
			if(i<nSyst-1) {
				cout << Form("%.7f",hpt_mid_NP[i]->GetBinContent(ipt+1)) << " & ";
			}
			else if(i==nSyst-1) {
				if(ipt<nPtmid-1) cout << Form("%.7f",hpt_mid_NP[i]->GetBinContent(ipt+1)) << " \\\\ "<< endl;
				if(ipt==nPtmid-1) cout << Form("%.7f",hpt_mid_NP[i]->GetBinContent(ipt+1)) << "\\\\ \\hline " << endl;
			}
		}
	}
	for (int iCent = 0; iCent < nCentBinmid; iCent++){
		for(int i=0; i<nSyst; i++){
			if(iCent==0 && i==0) cout << "0--10\\% 	&  \\multirow{6}{*}{6.5--40} &";
			else if(iCent!=0 && i==0) cout << Form("  %s \\%s  & & ", CentBinmidName[iCent], "%");
			if(i<nSyst-1) cout << Form("%.7f",hCent_mid_NP[i]->GetBinContent(iCent+1)) << " & ";
			else if(i==nSyst-1) {
				if(iCent<nCentBinmid-1) cout << Form("%.7f",hCent_mid_NP[i]->GetBinContent(iCent+1)) << " \\\\ " << endl;
				if(iCent==nCentBinmid-1) cout << Form("%.7f",hCent_mid_NP[i]->GetBinContent(iCent+1)) << "\\\\ \\hline" << endl;
			}
		}
	}
/*
	cout << "\\multicolumn{2}{c}{Non Prompt ~\\JPsi} \\\\ \\hline \\hline" << endl;
	for(int icent = 0; icent<nCent; icent++){
		for(int ipt = 0; ipt < nPt; ipt++){
			for(int i=0; i<nSyst; i++){
				if(ipt==0 && icent<nCent-1 && i==0) cout << Form("\\multirow{%d}{*}{%s\\%}	& %s	& ",nPt, RapTxt[icent],ptBinName[ipt]);
				else if(ipt==0 && icent==nCent-1 && i==0) cout << Form("\\multirow{}{*}{6.5--50} 	& %s	& ",CentBinName[ipt]);
				else if(ipt!=0 && icent<nCent-1 && i==0) cout << Form("				& %s	& ",ptBinName[ipt]);
				else if(ipt!=0 && ipt<6 && icent==nCent-1 && i==0) cout << Form("				& %s	& ",CentBinName[ipt]);
				if(i<nSyst-1 && icent<nCent-1){ cout << Form("%.7f",hpt_NP[icent][i]->GetBinContent(ipt+1)) << " & ";	}
				else if(i==nSyst-1 && icent<nCent-1){
					if(ipt<nPt-1) cout << Form("%.7f",hpt_NP[icent][i]->GetBinContent(ipt+1)) << " \\\\ " << endl;
					if(ipt==nPt-1) cout << Form("%.7f",hpt_NP[icent][i]->GetBinContent(ipt+1)) << " \\\\ \\hline " << endl;
				}
				else if(icent==1 && ipt<6){
					if(i<nSyst-1 && ipt<5) cout << Form("%.7f",hCent_NP[i]->GetBinContent(ipt+1)) << " & ";
					if(i<nSyst-1 && ipt==5) cout << Form("%.7f",hCent_NP[i]->GetBinContent(ipt+1)) << " & ";
					if(i==nSyst-1 && ipt<6)  cout << Form("%.7f",hCent_NP[i]->GetBinContent(ipt+1)) << "	\\\\ " << endl;
					if(i==nSyst-1 && ipt==6) cout << Form("%.7f",hCent_NP[i]->GetBinContent(ipt+1)) << "	\\\\ \\hline \\hline " << endl;
				}
			}
		}
	}


	const int nPtv3=5;

	cout << " +++++++++++++++++++++ v3 +++++++++++++++++++++++++ " << endl;
	cout << "\\multicolumn{2}{c}{Prompt ~\\JPsi} \\\\ \\hline \\hline" << endl;
	for(int icent = 0; icent<nCent; icent++){
		for(int ipt = 0; ipt < nPtv3; ipt++){
			for(int i=0; i<nSyst; i++){
				if(ipt==0 && icent==0 && i==0) cout << Form("\\multirow{%d}{*}{%s\\%}	& %s	& ",nPtv3, RapTxt[icent],ptBinv3Name[ipt]);
				else if(ipt==0 && icent==nCent-1 && i==0) cout << Form("\\multirow{6}{*}{6.5--50}	& %s	& ",CentBinv3Name[ipt]);
				else if(ipt!=0 && icent<nCent-1 && i==0) cout << Form("				& %s	& ",ptBinv3Name[ipt]);
				else if(ipt!=0 && ipt<3 && icent==nCent-1 && i==0) cout << Form("				& %s	& ",CentBinv3Name[ipt]);
				if(i<nSyst-1 && icent<nCent-1){
					cout << Form("%.7f",hpt_PRv3[icent][i]->GetBinContent(ipt+1)) << " & ";
				}
				else if(i==nSyst-1 && icent<nCent-1){
					if(ipt<nPtv3-1) cout << Form("%.7f",hpt_PRv3[icent][i]->GetBinContent(ipt+1)) << " \\\\ " << endl;
					if(ipt==nPtv3-1) cout << Form("%.7f",hpt_PRv3[icent][i]->GetBinContent(ipt+1)) << " \\\\ \\hline " << endl;
				}
				if(icent==1 && ipt<3){
					if(i<nSyst-1 && ipt<2) cout << Form("%.7f",hCent_PRv3[i]->GetBinContent(ipt+1)) << " & ";
					if(i<nSyst-1 && ipt==2) cout << Form("%.7f",hCent_PRv3[i]->GetBinContent(ipt+1)) << " & ";
					if(i==nSyst-1 && ipt<3)  cout << Form("%.7f",hCent_PRv3[i]->GetBinContent(ipt+1)) << "	\\\\ " << endl;
					if(i==nSyst-1 && ipt==3) cout << Form("%.7f",hCent_PRv3[i]->GetBinContent(ipt+1)) << "	\\\\ \\hline \\hline" << endl;
				}
			}
		}
	}

	cout << "\\multicolumn{2}{c}{Non Prompt ~\\JPsi} \\\\ \\hline \\hline" << endl;
	for(int icent = 0; icent<nCent; icent++){
		for(int ipt = 0; ipt < nPtv3; ipt++){
			for(int i=0; i<nSyst; i++){
				if(ipt==0 && icent<nCent-1 && i==0) cout << Form("\\multirow{%d}{*}{%s\\%}	& %s	& ",nPtv3, RapTxt[icent],ptBinv3Name[ipt]);
				else if(ipt==0 && icent==nCent-1 && i==0) cout << Form("\\multirow{}{*}{6.5--50} 	& %s	& ",CentBinv3Name[ipt]);
				else if(ipt!=0 && icent<nCent-1 && i==0) cout << Form("				& %s	& ",ptBinv3Name[ipt]);
				else if(ipt!=0 && ipt<6 && icent==nCent-1 && i==0) cout << Form("				& %s	& ",CentBinv3Name[ipt]);
				if(i<nSyst-1 && icent<nCent-1){ cout << Form("%.7f",hpt_NPv3[icent][i]->GetBinContent(ipt+1)) << " & ";	}
				else if(i==nSyst-1 && icent<nCent-1){
					if(ipt<nPtv3-1) cout << Form("%.7f",hpt_NPv3[icent][i]->GetBinContent(ipt+1)) << " \\\\ " << endl;
					if(ipt==nPtv3-1) cout << Form("%.7f",hpt_NPv3[icent][i]->GetBinContent(ipt+1)) << " \\\\ \\hline " << endl;
				}
				else if(icent==1 && ipt<3){
					if(i<nSyst-1 && ipt<3) cout << Form("%.7f",hCent_NPv3[i]->GetBinContent(ipt+1)) << " & ";
					if(i<nSyst-1 && ipt==3) cout << Form("%.7f",hCent_NPv3[i]->GetBinContent(ipt+1)) << " & ";
					if(i==nSyst-1 && ipt<3)  cout << Form("%.7f",hCent_NPv3[i]->GetBinContent(ipt+1)) << "	\\\\ " << endl;
					if(i==nSyst-1 && ipt==3) cout << Form("%.7f",hCent_NPv3[i]->GetBinContent(ipt+1)) << "	\\\\ \\hline \\hline " << endl;
				}
			}
		}
	}
*/

	/*
		 for(int icent = 0; icent<nCent+1;icent++){V
//for(int ipt = 0; ipt < nPt-1; ipt++){
for(int ipt = 0; ipt < nPt; ipt++){
for(int i=0; i<nSyst; i++){
if(ipt==0 && icent<nCent && i==0) cout << Form("\\multirow{3}{*}{%.f-%.f\\%%}     &  %s  & ",CentName[icent],CentName[icent+1],ptBinName[ipt]);
else if(ipt!=0 && i==0) cout << Form("                            &  %s  & ",ptBinName[ipt]);
if(i<nSyst-1 && icent<nCent){
if(ipt<nPt-1) cout << Form("%.8f",hpt[icent][i]->GetBinContent(ipt+1)) << " & ";
if(ipt==nPt-1) cout << Form("%.8f",hCent[nCent-1][i]->GetBinContent(icent+1)) << " & ";
}
else if(i==nSyst-1 && icent<nCent){
if(ipt<nPt-1) cout << Form("%.8f",hpt[icent][i]->GetBinContent(ipt+1)) << "  \\\\ " << endl;
if(ipt==nPt-1) cout << Form("%.8f",hCent[nCent-1][i]->GetBinContent(1+icent)) << "  \\\\ \\hline " << endl;
}
if(icent==nCent){
if(i<nSyst-1 && ipt<nPt-1) cout << Form("%.8f",hpt[nCent][i]->GetBinContent(ipt+1)) << " & ";
if(i<nSyst-1 && ipt==nPt-1) cout << Form("%.8f",hInt[i]->GetBinContent(1)) << " & ";
if(i==nSyst-1 && ipt<nPt-1) cout << Form("%.8f",hpt[nCent][i]->GetBinContent(ipt+1)) << "  \\\\ " << endl;
if(i==nSyst-1 && ipt==nPt-1) cout << Form("%.8f",hInt[i]->GetBinContent(1)) << "  \\\\ \\hline " << endl;
}
} 
}
}
*/








}
