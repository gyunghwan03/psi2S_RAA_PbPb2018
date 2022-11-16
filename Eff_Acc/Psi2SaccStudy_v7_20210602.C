#include <iostream>
#include "Style.h"
#include "tdrstyle.C"
#include "rootFitHeaders.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "TLegend.h"
#include "TString.h"

using namespace std;
using namespace RooFit;

bool IsAcceptable(double pt, double eta);

void PrintAcc(double bins[],TH1F* h); 

void Psi2SaccStudy_v7_20210602(int wtopt=0, TString rmk="20220125")
{
    //Basic Setting
    gStyle->SetOptStat(0);
    setTDRStyle();

    //READ Input Skimmed File
    TFile *rf;
    rf = new TFile("../OniaTree_Psi2SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
    TTree *tree = (TTree*) rf -> Get("hionia/myTree");

    TFile* wgtf1 = new TFile("ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","READ");
    TF1* tf1_wgt1 = (TF1*)wgtf1->Get("dataMC_Ratio1");

    //([0]+[1]*x+[2]*x*x+[4]*x*x*x)/((x-[3])*(x-[3])*(x-[3]))
    TFile* wgtf2 = new TFile("ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","READ");
    TF1* tf1_wgt2 = (TF1*)wgtf2->Get("dataMC_Ratio1");


    cout<<" Entry : "<<tree->GetEntries()<<endl;;

    Int_t           Gen_QQ_size;
    TClonesArray    *Gen_QQ_4mom;
    TClonesArray    *Gen_QQ_mupl_4mom;
    TClonesArray    *Gen_QQ_mumi_4mom;
    Float_t         Gen_QQ_ctau3D[1000];   //[Gen_QQ_size]
    Float_t         Gen_QQ_ctau[1000];   //[Gen_QQ_size]

    TBranch        *b_Gen_QQ_size;
    TBranch        *b_Gen_QQ_4mom;   //!
    TBranch        *b_Gen_QQ_mupl_4mom;   //!
    TBranch        *b_Gen_QQ_mumi_4mom;   //!
    TBranch        *b_Gen_QQ_ctau3D;   //!
    TBranch        *b_Gen_QQ_ctau;   //!



    Gen_QQ_4mom = 0;
    Gen_QQ_mupl_4mom = 0;
    Gen_QQ_mumi_4mom = 0;

    tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    tree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
    tree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
    tree->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
    tree->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);

    TLorentzVector* JP= new TLorentzVector;
    TLorentzVector* Mu1= new TLorentzVector;
    TLorentzVector* Mu2= new TLorentzVector;

    //double aDenPt_2021_allybin[] = {6.5,7.5,9.0,10.0,12.0,15.0,50.0};
    double aDenPt_2021_allybin[] = {6.5,7.5,8.5,10.0,12.0,14.0,16.0,18.0,20.0,25.0,50.0};


    TH1F *hDenPt_2021_ally = new TH1F("hDenPt_2021_ally",";p_{T} (GeV/c};",10,aDenPt_2021_allybin);
    TH1F *hNumPt_2021_ally = new TH1F("hNumPt_2021_ally",";p_{T} (GeV/c};",10,aDenPt_2021_allybin);
    TH1F *hAccPt_2021_ally = new TH1F("hAccPt_2021_ally",";p_{T} (GeV/c};",10,aDenPt_2021_allybin);

    //double aDenPt_2021_midybin[] = {6.5,7.5,9.0,10.0,12.0,15.0,50.0};
    double aDenPt_2021_midybin[] = {6.5,7.5,8.5,10.0,12.0,14.0,16.0,18.0,20.0,25.0,50.0};


    TH1F *hDenPt_2021_midy = new TH1F("hDenPt_2021_midy",";p_{T} (GeV/c};",10,aDenPt_2021_midybin);
    TH1F *hNumPt_2021_midy = new TH1F("hNumPt_2021_midy",";p_{T} (GeV/c};",10,aDenPt_2021_midybin);
    TH1F *hAccPt_2021_midy = new TH1F("hAccPt_2021_midy",";p_{T} (GeV/c};",10,aDenPt_2021_midybin);


    //  double aDenPt_2021_Forybin[] = {3.0,4.5,6.5,9.0,12.0,50.0};
    double aDenPt_2021_Forybin[] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,25.0,50.0};


    TH1F *hDenPt_2021_Fory = new TH1F("hDenPt_2021_Fory",";p_{T} (GeV/c};",14,aDenPt_2021_Forybin);
    TH1F *hNumPt_2021_Fory = new TH1F("hNumPt_2021_Fory",";p_{T} (GeV/c};",14,aDenPt_2021_Forybin);
    TH1F *hAccPt_2021_Fory = new TH1F("hAccPt_2021_Fory",";p_{T} (GeV/c};",14,aDenPt_2021_Forybin);


    hDenPt_2021_ally->Sumw2();
    hNumPt_2021_ally->Sumw2();
    hAccPt_2021_ally->Sumw2();

    hDenPt_2021_midy->Sumw2();
    hNumPt_2021_midy->Sumw2();
    hAccPt_2021_midy->Sumw2();

    hDenPt_2021_Fory->Sumw2();
    hNumPt_2021_Fory->Sumw2();
    hAccPt_2021_Fory->Sumw2();

    ///////////////////////////////////////
    hDenPt_2021_ally->SetMarkerStyle(24);
    hNumPt_2021_ally->SetMarkerStyle(24);
    hAccPt_2021_ally->SetMarkerStyle(24);

    hDenPt_2021_midy->SetMarkerStyle(26);
    hNumPt_2021_midy->SetMarkerStyle(26);
    hAccPt_2021_midy->SetMarkerStyle(26);

    hDenPt_2021_Fory->SetMarkerStyle(25);
    hNumPt_2021_Fory->SetMarkerStyle(25);
    hAccPt_2021_Fory->SetMarkerStyle(25);

    TH1F *hNocutY = new TH1F("hNocutY",";y;",200,-10,10);
    TH1F *hNocutPt = new TH1F("hNocutPt",";p_{T} (GeV/c);",200,0.0,100.0);
    TH2F *hNocut_2DPtY = new TH2F("hNocut_2DPtY",";y;p_{T} (GeV/c)",200,-10,10,200,0.0,100.0);

    TH1F *hDenY_Nocut = new TH1F("hDenY_Nocut",";y;",100,-2.5,2.5);
    TH1F *hNumY_Nocut = new TH1F("hNumY_Nocut",";y;",100,-2.5,2.5);
    TH1F *hAccY_Nocut = new TH1F("hAccY_Nocut",";y;",100,-2.5,2.5);

    TH1F *hDenPt_Nocut = new TH1F("hDenPt_Nocut",";p_{T} (GeV/c);",200,0.0,100.0);
    TH1F *hNumPt_Nocut = new TH1F("hNumPt_Nocut",";p_{T} (GeV/c);",200,0.0,100.0);
    TH1F *hAccPt_Nocut = new TH1F("hAccPt_Nocut",";p_{T} (GeV/c);",200,0.0,100.0);

    TH2F *hDen2D_Nocut = new TH2F("hDen2D_Nocut",";y;p_{T} (GeV/c)",100,-2.5,2.5,200,0.0,100.0);
    TH2F *hNum2D_Nocut = new TH2F("hNum2D_Nocut",";y;p_{T} (GeV/c)",100,-2.5,2.5,200,0.0,100.0);
    TH2F *hAcc2D_Nocut = new TH2F("hAcc2D_Nocut",";y;p_{T} (GeV/c)",100,-2.5,2.5,200,0.0,100.0);


    hDenY_Nocut->Sumw2();
    hNumY_Nocut->Sumw2();
    hAccY_Nocut->Sumw2();

    hDenPt_Nocut->Sumw2();
    hNumPt_Nocut->Sumw2();
    hAccPt_Nocut->Sumw2();

    hDen2D_Nocut->Sumw2();
    hNum2D_Nocut->Sumw2();
    hAcc2D_Nocut->Sumw2();

    Int_t nEvt = tree->GetEntries();
    cout << "nEvt : " << nEvt << endl;

    double pt = 0.0, y = 0.0, mu1_pt = 0.0, mu1_eta = 0.0;
    double mu2_pt = 0.0, mu2_eta = 0.0;
    double mass = 0.0, muPtCut = 0.0;

    //nEvt = 1000000;
    //Begin Loop
    for(int i=0; i<nEvt; i++){

        tree->GetEntry(i);
        if (i%100000==0) cout << ">>>>> EVENT " << i << " / " << tree->GetEntries() <<  endl;
        //cout<<"# of Gen QQ : "<<Gen_QQ_size<<endl;
        //if(Gen_QQ_size > 1) continue;
        for(int j = 0; j < Gen_QQ_size; j++){
            JP = (TLorentzVector*) Gen_QQ_4mom->At(j);
            Mu1 = (TLorentzVector*) Gen_QQ_mupl_4mom->At(j);
            Mu2 = (TLorentzVector*) Gen_QQ_mumi_4mom->At(j);
            //cout<<"dimuon pt : "<<JP->Pt()<<", y : "<<JP->Rapidity()<<endl;
            pt = JP->Pt();
            mass = JP->M();
            y = JP->Rapidity();
            mu1_pt = Mu1->Pt();
            mu1_eta = Mu1->Eta();
            mu2_pt = Mu2->Pt();
            mu2_eta = Mu2->Eta();
            hNocutY->Fill(y);
            hNocutPt->Fill(pt);
            hNocut_2DPtY->Fill(y,pt);

            if (fabs(y)>2.4) continue;          

            hDenY_Nocut->Fill(y);
            hDenPt_Nocut->Fill(pt);
            hDen2D_Nocut->Fill(y,pt);


            double wt1 = 1.0;
            double wt2 = 1.0;
            if (wtopt==1) {
                wt1 = tf1_wgt1->Eval(pt);
                wt2 = tf1_wgt2->Eval(pt);
            }


            hDenPt_2021_ally->Fill(pt,wt1); 
            if (1.6<=fabs(y) && fabs(y)<2.4) hDenPt_2021_Fory->Fill(pt,wt2);
            else if (fabs(y)<1.6) hDenPt_2021_midy->Fill(pt,wt1);
            bool mu1pass = IsAcceptable(mu1_pt,mu1_eta);
            bool mu2pass = IsAcceptable(mu2_pt,mu2_eta);

            if (mu1pass!=true || mu2pass!=true) continue;
            hNumPt_2021_ally->Fill(pt,wt1); 
            if (1.6<=fabs(y) && fabs(y)<2.4) hNumPt_2021_Fory->Fill(pt,wt2);
            else if (fabs(y)<1.6) hNumPt_2021_midy->Fill(pt,wt1);

            hNumY_Nocut->Fill(y);
            hNumPt_Nocut->Fill(pt);
            hNum2D_Nocut->Fill(y,pt);

        }

    }

    hAccY_Nocut->Divide(hNumY_Nocut,hDenY_Nocut,1,1,"B");
    hAccPt_Nocut->Divide(hNumPt_Nocut,hDenPt_Nocut,1,1,"B");
    hAcc2D_Nocut->Divide(hNum2D_Nocut,hDen2D_Nocut,1,1,"B");

    hAccPt_2021_ally->Divide(hNumPt_2021_ally,hDenPt_2021_ally,1,1,"B");
    hAccPt_2021_midy->Divide(hNumPt_2021_midy,hDenPt_2021_midy,1,1,"B");
    hAccPt_2021_Fory->Divide(hNumPt_2021_Fory,hDenPt_2021_Fory,1,1,"B");

    TCanvas* cpt = new TCanvas("cpt","",1200,400);
    cpt->Divide(3,1);
    cpt->cd(1);
    hDenPt_2021_ally->Draw("e");

    cpt->cd(2);
    hNumPt_2021_ally->Draw("e");

    cpt->cd(3);
    hAccPt_2021_ally->GetXaxis()->SetRangeUser(6.5,50.0);
    hAccPt_2021_ally->GetYaxis()->SetRangeUser(0.0,1.3);
    hAccPt_2021_ally->Draw("e");

    TLegend* leg_pt = new TLegend(0.65,0.65,0.90,0.90);
    leg_pt->SetBorderSize(0);
    leg_pt->SetFillColor(0);
    leg_pt->AddEntry(hAccPt_2021_ally,"|y| < 2.4","lp");
    leg_pt->Draw();

    cpt->SaveAs(Form("AccPt_AllY_wgt%d_%s.png",wtopt,rmk.Data()));

    TCanvas* cpt3 = new TCanvas("cpt3","",1200,400);
    cpt3->Divide(3,1);
    cpt3->cd(1);
    hDenPt_2021_midy->Draw("e");

    cpt3->cd(2);
    hNumPt_2021_midy->Draw("e");

    cpt3->cd(3);
    hAccPt_2021_midy->GetXaxis()->SetRangeUser(6.5,50.0);
    hAccPt_2021_midy->GetYaxis()->SetRangeUser(0.0,1.3);
    hAccPt_2021_midy->Draw("e");

    TLegend* leg_pt3 = new TLegend(0.65,0.65,0.90,0.90);
    leg_pt3->SetBorderSize(0);
    leg_pt3->SetFillColor(0);
    leg_pt3->AddEntry(hAccPt_2021_midy,"|y| < 1.6","lp");
    leg_pt3->Draw();

    cpt->SaveAs(Form("AccPt_MidY_wgt%d_%s.png",wtopt,rmk.Data()));


    TCanvas* cpt2 = new TCanvas("cpt2","",1200,400);
    cpt2->Divide(3,1);
    cpt2->cd(1);
    hDenPt_2021_Fory->Draw("e");

    cpt2->cd(2);
    hNumPt_2021_Fory->Draw("e");

    cpt2->cd(3);
    hAccPt_2021_Fory->GetXaxis()->SetRangeUser(3.0,50.0);
    hAccPt_2021_Fory->GetYaxis()->SetRangeUser(0.0,1.3);
    hAccPt_2021_Fory->Draw("e");

    TLegend* leg_pt2 = new TLegend(0.65,0.65,0.90,0.90);
    leg_pt2->SetBorderSize(0);
    leg_pt2->SetFillColor(0);
    leg_pt2->AddEntry(hAccPt_2021_Fory,"1.6 < |y| < 2.4","lp");
    leg_pt2->Draw();

    cpt2->SaveAs(Form("AccPt_ForY_wgt%d_%s.png",wtopt,rmk.Data()));


    TCanvas* cnocut = new TCanvas("cnocut","",1200,400);
    cnocut->Divide(3,1);
    cnocut->cd(1);
    hNocut_2DPtY->Draw("colz");
    cnocut->cd(2);
    gPad->SetLogy(1);

    hNocutPt->Draw("e");
    cnocut->cd(3);
    gPad->SetLogy(0);

    hNocutY->Draw("e");
    cnocut->SaveAs(Form("NoCutPtY_NoCut_wgt%d_%s.png",wtopt,rmk.Data()));

    TCanvas* cabsa = new TCanvas("cabsa","",1200,400);
    cabsa->Divide(3,1);
    cabsa->cd(1);
    hAcc2D_Nocut->Draw("colz");
    cabsa->cd(2);
    hAccPt_Nocut->GetYaxis()->SetRangeUser(0.0,1.0);
    hAccPt_Nocut->Draw("e");
    cabsa->cd(3);
    hAccY_Nocut->GetYaxis()->SetRangeUser(0.0,1.0);
    hAccY_Nocut->Draw("e");
    cabsa->SaveAs(Form("AccPtY_NoCut_wgt%d_%s.png",wtopt,rmk.Data()));


    std::cout << "--- p_{T} binning, |y| < 2.4" << std::endl; 
    PrintAcc(aDenPt_2021_allybin,hAccPt_2021_ally);
    std::cout << std::endl;
    std::cout << "--- p_{T} binning, |y| < 1.6" << std::endl; 
    PrintAcc(aDenPt_2021_midybin,hAccPt_2021_midy);
    std::cout << std::endl;
    std::cout << "--- p_{T} binning, 1.6 < |y| < 2.4" << std::endl; 
    PrintAcc(aDenPt_2021_Forybin,hAccPt_2021_Fory);
    std::cout << std::endl;


    TFile *wf = new TFile(Form("acceptance_Prompt_Psi2s_GenOnly_wgt%d_%s.root",wtopt,rmk.Data()),"RECREATE");
    wf->cd();

    hAccPt_2021_ally->Write();
    hAccPt_2021_midy->Write();
    hAccPt_2021_Fory->Write();


    wf->Write();
}

bool IsAcceptable(double pt, double eta){
    if ((fabs(eta) < 1.2 && pt>3.5) || ((fabs(eta) >=1.2 && fabs(eta) <2.1) && pt>(5.47-1.89*fabs(eta))) || ((fabs(eta) >=2.1 && fabs(eta) <2.4) && pt>1.5)) return true;


    else return false;
}

void PrintAcc(double bins[],TH1F* h) { 
    for (int i=0;i<h->GetNbinsX();i++) {
        std::cout << bins[i] << " - " << bins[i+1] << " : " << h->GetBinContent(i+1) << " +- " << h->GetBinError(i+1) << std::endl; 
    }
}

