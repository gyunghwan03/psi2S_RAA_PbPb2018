{

  gStyle->SetOptStat(0);
    //pT reweighting function
    TFile *fPtW1 = new TFile("ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","READ");
    TFile *fPtW2 = new TFile("ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
    TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
    TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");

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

    double np10 = 0.0, np11 = 0.0, np12 = 0.0, np13 = 0.0;
    double np20 = 0.0, np21 = 0.0, np22 = 0.0, np23 = 0.0;


    cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;
    cout<<"p11 : "<<p11<<", e11 : "<<e11<<", p21 : "<<p21<<", e21 : "<<e21<<endl;
    cout<<"p12 : "<<p12<<", e12 : "<<e12<<", p22 : "<<p22<<", e22 : "<<e22<<endl;
    cout<<"p13 : "<<p13<<", e13 : "<<e13<<", p23 : "<<p23<<", e23 : "<<e23<<endl;

    // sys for down
    np10 = p10 - e10;
    np11 = p11 - e11;
    np12 = p12 - e12;
    np13 = p13 - e13;

    np20 = p20 - e20;
    np21 = p21 - e21;
    np22 = p22 - e22;
    np23 = p23 - e23;

    cout<<"np10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;
    cout<<"np11 : "<<np11<<", e11 : "<<e11<<", np21 : "<<np21<<", e21 : "<<e21<<endl;
    cout<<"np12 : "<<np12<<", e12 : "<<e12<<", np22 : "<<np22<<", e22 : "<<e22<<endl;
    cout<<"np13 : "<<np13<<", e13 : "<<e13<<", np23 : "<<np23<<", e23 : "<<e23<<endl;



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


    cout<<"np10 : "<<np10<<", e10 : "<<e10<<", np20 : "<<np20<<", e20 : "<<e20<<endl;
    cout<<"np11 : "<<np11<<", e11 : "<<e11<<", np21 : "<<np21<<", e21 : "<<e21<<endl;
    cout<<"np12 : "<<np12<<", e12 : "<<e12<<", np22 : "<<np22<<", e22 : "<<e22<<endl;
    cout<<"np13 : "<<np13<<", e13 : "<<e13<<", np23 : "<<np23<<", e23 : "<<e23<<endl;


    //cout<<"pt : 3.07271, "<<"weight : "<<fptw1sys->Eval(3.07271)<<endl;
    cout<<"pt : 3.07271, "<<"weight : "<<fptw2sys->Eval(3.07271)<<endl;
    TLegend *l = new TLegend(0.5,0.5,0.82,0.63);
    l->SetFillColor(0);
    l->SetBorderSize(0);
    l->SetMargin(0.2);
    l->SetTextSize(0.029);
    l->SetTextFont(42);
    l->AddEntry(fptw1sys,"Down all","PL");
    l->AddEntry(fptw2sys,"Down forward","PL");

    TH1F *h1 = new TH1F("h1",";p_{T} (GeV/c);",10,3.0,6.0);
    h1->GetYaxis()->SetRangeUser(-0.2,0.5);
    h1->GetXaxis()->SetRangeUser(3,6);

    h1->Draw();
    fptw1sys->SetLineColor(kBlue);
    fptw2sys->SetLineColor(kRed);
    //fptw2sys->GetYaxis()->SetRangeUser(-0.2,0.5);
    //fptw2sys->GetXaxis()->SetRangeUser(3.0,6.0);
    fptw2sys->Draw("same");
    //fptw1sys->Draw("same");
    //l->Draw("same");

    TLine *l1 = new TLine(3.0, 0.0, 6.0, 0.0);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    //l1->Draw("same");
}
