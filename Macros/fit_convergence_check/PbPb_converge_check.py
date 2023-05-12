import os
import ROOT
from ROOT import TFile, RooFit

#dir_path = '/Users/pjgwak/work/psi2S_Raa_Run2018_at2023/Macros/psi2S/'# need slash at last
dir_path = '../psi2S_230512/'# need slash at last

def make_pt_file_list(rapi_low, rapi_high, cent_low, cent_high, bin_tags, file_list):
    # forward pT dependence
    for pt_low, pt_high in zip(bin_tags, bin_tags[1:]):
        mc_mass = mass = f'mc_MassFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        mass = f'Mass_FixedFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_err = f'CtauErrResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_res= f'CtauResResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_bkg= f'CtauBkgResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_true= f'CtauTrueResult_Inclusive_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}.root'
        final= f'2DFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        file_list.append(mc_mass)
        file_list.append(mass)
        file_list.append(ctau_err)
        file_list.append(ctau_res)
        file_list.append(ctau_bkg)
        file_list.append(ctau_true)
        file_list.append(final)

def make_cent_file_list(rapi_low, rapi_high, pt_low, pt_high, bin_tags, file_list):
    # forward pT dependence
    for cent_low, cent_high in zip(bin_tags, bin_tags[1:]):
        mc_mass = f'mc_Mass_FitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        mass = f'Mass_FixedFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_err = f'CtauErrResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_res= f'CtauResResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_bkg= f'CtauBkgResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_true= f'CtauTrueResult_Inclusive_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}.root'
        final= f'2DFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_centrality{cent_low}-{cent_high}_PRw_Effw0_Accw0_PtW0_TnP0.root'
        file_list.append(mc_mass)
        file_list.append(mass)
        file_list.append(ctau_err)
        file_list.append(ctau_res)
        file_list.append(ctau_bkg)
        file_list.append(ctau_true)
        file_list.append(final)

def check_fit_status(file_name):
    # check fitting step
    if 'mc_MassFitResult' in file_name:
        step = 'roots_MC/Mass/'
        fit_store = 'fitresult_pdfMASS_Tot_dsAB'
    elif 'Mass_FixedFitResult' in file_name:
        step = 'roots/2DFit_No_Weight/Mass/'
        fit_store = 'fitresult_pdfMASS_Tot_dsAB'
    elif 'CtauErrResult' in file_name:
        step = 'roots/2DFit_No_Weight/CtauErr/'
        fit_store = 'sPlot doesn\'t make a fit'
    elif 'CtauResResult' in file_name:
        step = 'roots/2DFit_No_Weight/CtauRes/'
        fit_store = 'fitresult_GaussModel_Tot_ctauResCutDS'
    elif 'CtauBkgResult' in file_name:
        step = 'roots/2DFit_No_Weight/CtauBkg/'
        fit_store = 'fitresult_pdfTot_Bkg_dataw_Bkg'
    elif 'CtauTrueResult' in file_name:
        step = 'roots/2DFit_No_Weight/CtauTrue/'
        fit_store = 'fitresult_TrueModel_Tot_reducedDS_MC'
    elif '2DFitResult' in file_name:
        step = 'roots/2DFit_No_Weight/Final/'
        fit_store = 'fitresult_pdfCTAUMASS_Tot_dsToFit'
    else:
        step = 'No_file'
        #print('Error: File missed')


    # read fit result
    file_path  = dir_path + step + file_name
    is_file_exist = os.path.isfile(file_path)
    if not is_file_exist:
        return ('No', 'File')
    elif step == 'roots/2DFit_No_Weight/CtauErr/':
        return ('sPlot', 'no fit')
    else:
        infile = TFile(file_path)
        fit_result = infile.Get(fit_store)
        hesse_code = fit_result.status() # Hesse code
        edm_temp = format(fit_result.edm(), '.5E') # Define edm_temp for formatting
        edm = f'{edm_temp}\t'
        mll = fit_result.minNll()
        fit_return = [hesse_code, edm, format(mll,'.5E')]
        
        floating_list = fit_result.floatParsFinal()
        for idx in range(len(floating_list)):
            fit_return.append(str(floating_list[idx]).strip())    
        return fit_return


def main():
    # list to make bin list
    forward_pt = ['3.0', '6.5', '12.0', '50.0']
    forward_cent = ['0', '40', '80', '180']
    mid_pt = ['6.5', '9.0', '12.0', '15.0', '20.0', '50.0']
    mid_cent = ['0', '20', '40', '60', '80', '100', '180']

    # list to save file names
    forward_pt_files = []
    forward_cent_files =[]
    mid_pt_files = []
    mid_cent_files = []


    # get input file names - forward
    make_pt_file_list('1.6', '2.4', '0', '180', forward_pt, forward_pt_files)
    make_cent_file_list('1.6', '2.4', '3.0', '50.0', forward_cent, forward_cent_files)

    # mid rapidity
    make_pt_file_list('0.0', '1.6', '0', '180', mid_pt, mid_pt_files)
    make_cent_file_list('0.0', '1.6', '6.5', '50.0', mid_cent, mid_cent_files)


    # get and print fit results
    fit_result_list = []
    #print("file\tHesse code\tedm")
    for file_name in forward_pt_files:
        fit_result_list.append(check_fit_status(file_name))
    for file, fit_result, in zip(forward_pt_files, fit_result_list):
        printable_result = f'{file}'
        for parameter in fit_result:
            printable_result += '\t' + f'{parameter}'
        print(printable_result)
    fit_result_list.clear()
    print("")

    for file_name in forward_cent_files:
        fit_result_list.append(check_fit_status(file_name))
    for file, fit_result, in zip(forward_cent_files, fit_result_list):
        printable_result = f'{file}'
        for parameter in fit_result:
            printable_result += '\t' + f'{parameter}'
        print(printable_result)
    fit_result_list.clear()
    print("")

    for file_name in mid_pt_files:
        fit_result_list.append(check_fit_status(file_name))
    for file, fit_result, in zip(mid_pt_files, fit_result_list):
        printable_result = f'{file}'
        for parameter in fit_result:
            printable_result += '\t' + f'{parameter}'
        print(printable_result)
    fit_result_list.clear()
    print("")

    for file_name in mid_cent_files:
        fit_result_list.append(check_fit_status(file_name))
    for file, fit_result, in zip(mid_cent_files, fit_result_list):
        printable_result = f'{file}'
        for parameter in fit_result:
            printable_result += '\t' + f'{parameter}'
        print(printable_result)


if __name__ == '__main__':
    main()