static char *__tty_color[] = {
    "\033[0;40;30m",         /* 0   black on black */
    "\033[0;40;31m",         /* 1   red */
    "\033[0;40;32m",         /* 2   green */
    "\033[0;40;33m",         /* 3   brown */
    "\033[0;40;34m",         /* 4   blue */
    "\033[0;40;35m",         /* 5   magenta */
    "\033[0;40;36m",         /* 6   cyan */
    "\033[0;40;37m",         /* 7   light gray */
    "\033[1;40;30m",         /* 0   gray */
    "\033[1;40;31m",         /* 1   brightred */
    "\033[1;40;32m",         /* 2   brightgreen */
    "\033[1;40;33m",         /* 3   yellow */
    "\033[1;40;34m",         /* 4   brightblue */
    "\033[1;40;35m",         /* 5   brighmagenta */
    "\033[1;40;36m",         /* 6   brightcyan */
    "\033[1;40;37m",         /* 7   white */
};
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
            auto fitresult = (RooFitResult*)f->Get("fitresult_pdfCTAUMASS_Tot_dsToFit");
            // fitresult->Print();
            // auto data = (RooDataSet*)f->Get("datasetMass");
            // auto a = data->get();
			if( fitresult->statusCodeHistory(0)!=0 || fitresult->statusCodeHistory(1)!=0) {
                cout <<  __tty_color[1] << "" << endl;
                cout << fileName << endl;
                cout << "MINIMIZE ERROR FLAG : " << fitresult->statusCodeHistory(0) << endl;
                cout << "HESSE ERROR FLAG : " << fitresult->statusCodeHistory(1) << endl;
            }
            else {
                cout <<  __tty_color[7] << "" << endl;
                cout << fileName << endl;
                cout << "MINIMIZE ERROR FLAG : " << fitresult->statusCodeHistory(0) << endl;
                cout << "HESSE ERROR FLAG : " << fitresult->statusCodeHistory(1) << endl;
            }
            // bool check = false;
            // for(auto i : *a){
            //     RooRealVar *h;
            //     h = dynamic_cast<RooRealVar*>(a->find(i->GetName()));
            //     TString results = Form("Variable [%s] // Final Value : [%f] // Lower Limit : [%f] // Upper Limit : [%f]",h->GetName(),h->getVal(),h->getRange().first,h->getRange().second);
            //     cout << results << endl;
                // cout << "Variable Name : " << h->GetName() << " Final Value : " << h->getVal() << " Upper Limit : " << h->getRange().second << " Lower Limit : "  <<  h->getRange().first << endl;
                // if(h->getVal() == h->getRange().first|| h->getVal() == h->getRange().second) cout << "Var "<<i->GetName() << " reaches the limit" << endl;
                // h->delete();
                // check = true;:w

            // }
            cout << "\n"<< endl;
            // if(check) cout << fileName << endl;
            // cout << fitConvCheck(fitresult,data) << endl;
            // if(fitConvCheck(fitresult,data)) break;
            }
            
}
