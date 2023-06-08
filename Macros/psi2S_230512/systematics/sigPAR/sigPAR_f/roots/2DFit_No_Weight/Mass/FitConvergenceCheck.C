void FitConvergenceCheck(){
    TString folderPath = "."; // 가져올 폴더 경로로 변경해주세요
    TSystemDirectory dir(".", folderPath);
    TList* fileList = dir.GetListOfFiles();
        TSystemFile* file;
        TIter next(fileList);
            while ((file = (TSystemFile*)next())) {
            TString fileName = file->GetName();
            if(fileName == "." || fileName==".." || fileName=="FitConvergenceCheck.C") continue;
            TFile *f = TFile::Open(fileName.Data());
            auto fitresult = (RooFitResult*)f->Get("fitresult_pdfMASS_Tot_dsAB");
            // fitresult->Print();
            auto data = (RooDataSet*)f->Get("datasetMass");
            auto a = data->get();
            cout << fileName << endl;
            cout << "MINIMIZE ERROR FLAG : " << fitresult->statusCodeHistory(0) << endl;
            cout << "HESSE ERROR FLAG : " << fitresult->statusCodeHistory(1) << endl;
            // bool check = false;
            for(auto i : *a){
                RooRealVar *h;
                h = dynamic_cast<RooRealVar*>(a->find(i->GetName()));
                TString results = Form("Variable [%s] // Final Value : [%f] // Error : [%f] // Lower Limit : [%f] // Upper Limit : [%f]",h->GetName(),h->getVal(),h->getError(),h->getRange().first,h->getRange().second);
                cout << results << endl;
                // cout << "Variable Name : " << h->GetName() << " Final Value : " << h->getVal() << " Upper Limit : " << h->getRange().second << " Lower Limit : "  <<  h->getRange().first << endl;
                // if(h->getVal() == h->getRange().first|| h->getVal() == h->getRange().second) cout << "Var "<<i->GetName() << " reaches the limit" << endl;
                // h->delete();
                // check = true;
            }
            cout << "\n"<< endl;
            // if(check) cout << fileName << endl;
            // cout << fitConvCheck(fitresult,data) << endl;
            // if(fitConvCheck(fitresult,data)) break;
            }
            
}
