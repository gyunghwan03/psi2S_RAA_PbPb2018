void check_convergence()
{
	TString folderPath=".";
	TSystemDirectory dir(".", folderPath);
	TList* fileList = dir.GetListOfFiles();
	TSystemFile* file;
	TIter next(fileList);
	while ((file = (TSystemFile*)next())) {
		TString fileName = file->GetName();
		if(fileName == "." || fileName==".." || fileName=="FitConvergenceCheck.C" || fileName==".check_convergence.C.swp" || fileName=="check_convergence.C") continue;
		TFile *f = TFile::Open(fileName.Data());
		auto fit_result = (RooFitResult*)f->Get("fitresult_pdfMASS_Tot_dsAB");
		// fitresult->Print();
		auto data = (RooDataSet*)f->Get("datasetMass");

		int hesse_code = fit_result->status();
		double edm = fit_result->edm();
		double mll = fit_result->minNll();
		//RooRealVar* par_fitresult = (RooRealVar*)fit_result->floatParsFinal().find("N_Jpsi");
		RooRealVar* fit_para = (RooRealVar*)fit_result->floatParsFinal().at(0);
		double val_ = fit_para->getVal();
		double err_ = fit_para->getError();
		double min_ = fit_para->getMin();
		double max_ = fit_para->getMax();

		//cout << "\n###### Fit Convergence Check ######\n";
		//cout << "Caveat: mean was fixed so it's not a matter" << endl;
		//cout << "Hesse: " << hesse_code << "\t edm: " << edm << "\t mll: " << mll << endl;
		int cnt_ = 0;
		for (int idx = 0; idx < fit_result->floatParsFinal().getSize(); idx++) {
			RooRealVar *fit_para = (RooRealVar *)fit_result->floatParsFinal().at(idx);
			double val_ = fit_para->getVal();
			double err_ = fit_para->getError();
			double min_ = fit_para->getMin();
			double max_ = fit_para->getMax();
			if (hesse_code!=0||(val_-err_<min_)||(val_+err_>max_)){
				cout << " " << endl;
				cout << "\033[31m" << endl;
				cout << "File Name : " << fileName << endl;
				cout << "Hesse : " << hesse_code << endl;
                cout << "[Stuck] " << fit_para->GetName() << " Final Value : " << val_ << " (" << min_ << " ~ " << max_ << ")" << endl << endl;
				cout << "\033[0m" << endl;
                cnt_++;
			}
			//if ((val_ - err_ > min_) && (val_ + err_ < max_)) {
			//	// No work is intended
			//}
			//else {
			//	cout << "File Name : " << fileName << endl;
			//	cout << "\033[31m" << "[Stuck] " << fit_para->GetName() << " \033[0m // Final Value : " << val_ << " (" << min_ << " ~ " << max_ << ")" << endl << endl;
			//	cnt_++;
			//}
		}
		if (cnt_ == 0) cout << "[Fit Converged]" << endl;
	}
}
