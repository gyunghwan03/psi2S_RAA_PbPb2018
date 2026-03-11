{

    //pT reweighting function
    TFile *fPtW1 = new TFile("ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","READ");
    TFile *fPtW2 = new TFile("ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
    TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
    TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");

    double corerr = 1.0; //0.85;

    // for systematics
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

    cout<<"forward, p20: "<<p20<<", p21 : "<<p21<<", p22 : "<<p22<<", p23 : "<<p23<<endl; 

    TF1 *fptw1sys = new TF1("fptw1sys","( [0]*x + [1]*x*x +[3]*x*x*x ) / (  (x-[2])*(x-[2])*(x-[2])  )",6.5, 50.0);
    TF1 *fptw2sys = new TF1("fptw2sys","( [0]*x + [1]*x*x +[3]*x*x*x ) / (  (x-[2])*(x-[2])*(x-[2])  )",3.0, 50.0);
    fptw2sys->SetParameter(0, p20);
    fptw2sys->SetParameter(1, p21);
    fptw2sys->SetParameter(2, p22);
    fptw2sys->SetParameter(3, p23);

    fptw2sys->Draw();
    fptw2->Draw("same");


    double np10 = 0.0, np11 = 0.0, np12 = 0.0, np13 = 0.0;
    double np20 = 0.0, np21 = 0.0, np22 = 0.0, np23 = 0.0;


    cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

    TF1 *sysup;
    TF1 *sysdo;
    TF1 *sysno;

    TF1 *sysup1;
    TF1 *sysdo1;
    TF1 *sysno1;

    int isPtWgtUp = 0;
    for(int i = 0; i < 3; i++){
        isPtWgtUp = i;

        // sys for nominal
        if(isPtWgtUp == 0) {
            cout<<"isPt : "<<isPtWgtUp<<endl;
            np10 = p10;
            np11 = p11;
            np12 = p12;
            np13 = p13;

            np20 = p20;
            np21 = p21;
            np22 = p22;
            np23 = p23;
            cout<<"i : "<<isPtWgtUp<<", np10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;
        }

        // sys for up
        if(isPtWgtUp == 1) {
            cout<<"isPt : "<<isPtWgtUp<<endl;
            np10 = p10 + e10;
            np11 = p11 + e11;
            np12 = p12 + e12;
            np13 = p13 + e13;

            np20 = p20 + e20;
            np21 = p21 + e21;
            np22 = p22 + e22;
            np23 = p23 + e23;
            cout<<"i : "<<isPtWgtUp<<", np10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;
        }

        // sys for down
        if(isPtWgtUp == 2) {
            cout<<"isPt : "<<isPtWgtUp<<endl;
            np10 = p10 - e10*corerr;
            np11 = p11 - e11*corerr;
            np12 = p12 - e12*corerr;
            np13 = p13 - e13*corerr;

            np20 = p20 - e20*corerr;
            np21 = p21 - e21*corerr;
            np22 = p22 - e22*corerr;
            np23 = p23 - e23*corerr;
            cout<<"i : "<<isPtWgtUp<<", np10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;
        }

        cout<<"final : "<<isPtWgtUp<<", p10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;

        fptw1sys->SetParameter(0, np10);
        fptw1sys->SetParameter(1, np11);
        fptw1sys->SetParameter(2, np12);
        fptw1sys->SetParameter(3, np13);

        fptw2sys->SetParameter(0, np20);
        fptw2sys->SetParameter(1, np21);
        fptw2sys->SetParameter(2, np22);
        fptw2sys->SetParameter(3, np23);

        np10 = p10;
        np11 = p11;
        np12 = p12;
        np13 = p13;

        np20 = p20;
        np21 = p21;
        np22 = p22;
        np23 = p23;


        cout<<"final 2 : "<<isPtWgtUp<<", p10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;

        if(i == 0) sysno = (TF1*)fptw1sys->Clone();
        if(i == 1) sysup = (TF1*)fptw1sys->Clone();
        if(i == 2) sysdo = (TF1*)fptw1sys->Clone();

        if(i == 0) sysno1 = (TF1*)fptw2sys->Clone();
        if(i == 1) sysup1 = (TF1*)fptw2sys->Clone();
        if(i == 2) sysdo1 = (TF1*)fptw2sys->Clone();


    }

    cout<<""<<endl;
    cout<<"final 3 nominal 1 : "<<sysno->GetParameter(0)<<endl;
    cout<<"final 3 up 1 : "<<sysup->GetParameter(0)<<endl;
    cout<<"final 3 down 1 : "<<sysdo->GetParameter(0)<<endl;

    cout<<""<<endl;

    cout<<"final 3 nominal 2 : "<<sysno1->GetParameter(0)<<endl;
    cout<<"final 3 up 2 : "<<sysup1->GetParameter(0)<<endl;
    cout<<"final 3 down 2 : "<<sysdo1->GetParameter(0)<<endl;


    sysno->SetLineColor(kBlue);
    sysdo->SetLineColor(kRed);
    sysdo->SetLineStyle(2);
    sysup->SetLineStyle(3);
    sysup->SetLineColor(kViolet);

    sysno1->SetLineColor(kBlue);
    sysdo1->SetLineColor(kRed);
    sysup1->SetLineColor(kViolet);
    sysdo1->SetLineStyle(2);
    sysup1->SetLineStyle(3);

    TLegend *l = new TLegend(0.5,0.5,0.82,0.63);
    l->SetFillColor(0);
    l->SetBorderSize(0);
    l->SetMargin(0.2);
    l->SetTextSize(0.029);
    l->SetTextFont(42);
    l->AddEntry(sysno,"Nominal","PL");
    l->AddEntry(sysup,"1 Sigma up","PL");
    l->AddEntry(sysdo,"1 Sigma down","PL");

    TCanvas *c1 = new TCanvas("c1","",1000,500);
    c1->Divide(2,1);

    c1->cd(1);

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();
    lt1->SetTextSize(0.035);
    sysup->GetYaxis()->SetRangeUser(-0.2,10.0);
    sysup->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    sysup->GetXaxis()->CenterTitle();
    sysup->Draw();
    //sysno->Draw();
    sysno->Draw("same");
    sysdo->Draw("same");
    lt1->DrawLatex(0.3,0.7,"|y| < 2.4 && p_{T} > 6.5 GeV/c");
    l->Draw("same");

    c1->cd(2);
    sysup1->GetYaxis()->SetRangeUser(-0.2,2.0);
    sysup1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    sysup1->GetXaxis()->CenterTitle();
    sysup1->Draw();
    //sysno->Draw();
    sysno1->Draw("same");
    sysdo1->Draw("same");
    //fptw2->Draw("same");
    lt1->DrawLatex(0.3,0.7,"1.6 < |y| < 2.4 && p_{T} > 3.0 GeV/c");
    l->Draw("same");


    TLine *l1 = new TLine(3.0, 0.0, 50.0, 0.0);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->Draw("same");


    c1->SaveAs("weighting_factor_sys.png");
    cout<<"pt : 3.07271, "<<"weight : "<<fptw2sys->Eval(3.07271)<<endl;

}
