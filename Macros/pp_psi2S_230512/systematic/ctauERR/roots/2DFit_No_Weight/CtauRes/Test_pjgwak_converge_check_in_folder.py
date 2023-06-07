import os
import argparse
import subprocess as sub
import ROOT
from ROOT import TFile, RooFit
parser = argparse.ArgumentParser(description='')
parser.add_argument('--target', required=True, help='인풋 파일')
args = parser.parse_args()

dir_path = str(args.target)
infile = TFile(dir_path)
fit_result = infile.Get('fitresult_GaussModel_Tot_ctauResCutDS')
hesse_code = fit_result.status() # Hesse code
edm = fit_result.edm()
mll = fit_result.minNll()
fit_return = [hesse_code, format(edm,'.5E'), format(mll,'.5E')]
floating_list = fit_result.floatParsFinal()
for idx in range(len(floating_list)):
    if (floating_list[idx].getVal() - floating_list[idx].getError() > floating_list[idx].getRange()[0]) and (floating_list[idx].getVal() + floating_list[idx].getError() < floating_list[idx].getRange()[1]) :
        a = '[Okay]'
    else:
        a = '\033[31m' + '[Stuck]'
    print(f'{a} // {floating_list[idx].GetName()} \033[0m // Final Value : {floating_list[idx].getVal(): .4f} // Lower Limit : {floating_list[idx].getRange()[0]} // Upper Limit : {floating_list[idx].getRange()[1]}')
    # "\033[0m" is system command for cosmetic (No color in terminal)


exit(1)
def make_pt_file_list(rapi_low, rapi_high, bin_tags, file_list):
    for pt_low, pt_high in zip(bin_tags, bin_tags[1:]):
        mc_mass = mass = f'mc_MassFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        mass = f'Mass_FixedFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_err = f'CtauErrResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_res= f'CtauResResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_bkg= f'CtauBkgResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_true= f'CtauTrueResult_Inclusive_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0.root'
        final= f'2DFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        file_list.append(mc_mass)
        file_list.append(mass)
        file_list.append(ctau_err)
        file_list.append(ctau_res)
        file_list.append(ctau_bkg)
        file_list.append(ctau_true)
        file_list.append(final)

def make_cent_file_list(rapi_low, rapi_high, pt_low, pt_high, bin_tags, file_list):
    for cent_low, cent_high in zip(bin_tags, bin_tags[1:]):
        mc_mass = f'mc_MassFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        mass = f'Mass_FixedFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_err = f'CtauErrResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_res= f'CtauResResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_bkg= f'CtauBkgResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
        ctau_true= f'CtauTrueResult_Inclusive_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0.root'
        final= f'2DFitResult_pt{pt_low}-{pt_high}_y{rapi_low}-{rapi_high}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root'
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
        #print('Error: No step matched')


    # read fit result
    dir_path  = dir_path + step + file_name
    is_file_exist = os.path.isfile(file_path)
    if not is_file_exist:
        return ('No', 'File')
    elif step == 'roots/2DFit_No_Weight/CtauErr/':
        return ('sPlot', 'no fit')
    else:
        infile = TFile(file_path)
        fit_result = infile.Get(fit_store)
        hesse_code = fit_result.status() # Hesse code
        edm = fit_result.edm()
        mll = fit_result.minNll()
        fit_return = [hesse_code, format(edm,'.5E'), format(mll,'.5E')]
        
        floating_list = fit_result.floatParsFinal()
        for idx in range(len(floating_list)):
            fit_return.append(str(floating_list[idx]).strip())    
        return fit_return


def main():
    # list to make bin list
    forward_pt = ['3.0', '6.5', '12.0', '50.0']
    forward_cent = ['', ''] # pp doesn't have centrality. Leave two empty string to run this code.
    mid_pt = ['6.5', '9.0', '12.0', '15.0', '50.0']
    mid_cent = ['', '']

    # list to save file names
    forward_pt_files = []
    forward_cent_files =[]
    mid_pt_files = []
    mid_cent_files = []


    # get input file names - forward
    make_pt_file_list('1.6', '2.4', forward_pt, forward_pt_files)
    make_cent_file_list('1.6', '2.4', '3.0', '50.0', forward_cent, forward_cent_files)

    # mid rapidity
    make_pt_file_list('0.0', '1.6', mid_pt, mid_pt_files)
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