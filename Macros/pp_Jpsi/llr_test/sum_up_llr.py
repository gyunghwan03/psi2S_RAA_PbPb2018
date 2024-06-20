# ===== Skeleton ===== #
# 빈 별로
# 결과물 모양은 같게
# Forward
# 3.5 - 6.5
# 6.5 - 9
# 9 -12
# 12 - 50
# 3.5 - 50

# Mid
# 6.5 - 9
# 9 - 12
# 12 - 15
# 15 - 20
# 20 - 25
# 25 - 50
# 6.5 - 50

# 빈 설정 -> basic_path + list + loop
# 1,2,3,4 차 계산 해서 기억 -> fcn, root code 실행 + loop
# 표 만들고 저장

# ==================== #

import subprocess as sub

def cal_p_value(nll_low_dimension, nll_high_dimension, diff_dimension):
    # Call .C file
    fit_command = f"root -l -b -q cal_p_value.C\'({nll_low_dimension},{nll_high_dimension},{diff_dimension})\'"
    process = sub.Popen(fit_command, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
    output, error = process.communicate()
    p_value = float(output.decode().strip().split()[-1])

    # Below one line is added to decrease length of p-value
    # Long numbers make the table of LaText get cut.
    p_value_str = '{:.2f}'.format(p_value)

    return p_value_str


# Set the kinematics bins
#fwd_pt = [(6.5, 9.0)] # For run test
fwd_pt = [(3.5, 6.5), (6.5, 9.0), (9.0, 12.0), (12.0, 50.0), (3.5, 50.0)]
mid_pt = [(6.5, 9.0), (9.0, 12.0), (12.0, 15.0), (15.0, 20.0), (20.0, 25.0), (25.0, 50.0), (6.5, 50.0)]
rapidity = [(0, 1.6), (1.6, 2.4)]
bkg_names = ['Expo', '1stCheby', '2ndCheby', '3rdCheby', '4thCheby', '5thCheby', '6thCheby']

in_dir_path = './figs/2DFit_LLR/Mass'  # Input file's directory


# Forward
# 빈 별로 계산해서 저장하기
for pt_low, pt_high in fwd_pt:
    nlls = {}
    # Fill the nll values
    for bkg_name in bkg_names:
        in_file_name = f'Mass_Inclusive_pt{pt_low}-{pt_high}_y1.6-2.4_muPt0.0_{bkg_name}.txt'
        input_path = in_dir_path + '/' + in_file_name

        # Read NLL from input
        with open(input_path, 'r') as in_file:
            content = in_file.read()
        #print(content)
        nlls[f'{bkg_name}'] = content
    #print(nlls)

    # Calculate p-values
    # Call .C file
    p_values = {'0_0':'0', }
    
    # 1st -> vs 1 -> 0
    p_values['0_1'] = cal_p_value(nlls['Expo'], nlls['1stCheby'], 1)
    
    # 2nd -> vs 0, 1
    p_values['0_2'] = cal_p_value(nlls['Expo'], nlls['2ndCheby'], 2)
    p_values['1_2'] = cal_p_value(nlls['1stCheby'], nlls['2ndCheby'], 1)

    # 3rd -> vs 0, 1, 2
    p_values['0_3'] = cal_p_value(nlls['Expo'], nlls['3rdCheby'], 3)
    p_values['1_3'] = cal_p_value(nlls['1stCheby'], nlls['3rdCheby'], 2)
    p_values['2_3'] = cal_p_value(nlls['2ndCheby'], nlls['3rdCheby'], 1)
    
    # 4th -> vs 0, 1, 2, 3
    p_values['0_4'] = cal_p_value(nlls['Expo'], nlls['4thCheby'], 4)
    p_values['1_4'] = cal_p_value(nlls['1stCheby'], nlls['4thCheby'], 3)
    p_values['2_4'] = cal_p_value(nlls['2ndCheby'], nlls['4thCheby'], 2)
    p_values['3_4'] = cal_p_value(nlls['3rdCheby'], nlls['4thCheby'], 1)
    
    # 5th -> vs 0, 1, 2, 3. 4
    p_values['0_5'] = cal_p_value(nlls['Expo'], nlls['5thCheby'], 5)
    p_values['1_5'] = cal_p_value(nlls['1stCheby'], nlls['5thCheby'], 4)
    p_values['2_5'] = cal_p_value(nlls['2ndCheby'], nlls['5thCheby'], 3)
    p_values['3_5'] = cal_p_value(nlls['3rdCheby'], nlls['5thCheby'], 2)
    p_values['4_5'] = cal_p_value(nlls['4thCheby'], nlls['5thCheby'], 1)

    # 6th -> vs 0, 1, 2, 3. 4, 5
    p_values['0_6'] = cal_p_value(nlls['Expo'], nlls['6thCheby'], 6)
    p_values['1_6'] = cal_p_value(nlls['1stCheby'], nlls['6thCheby'], 5)
    p_values['2_6'] = cal_p_value(nlls['2ndCheby'], nlls['6thCheby'], 4)
    p_values['3_6'] = cal_p_value(nlls['3rdCheby'], nlls['6thCheby'], 3)
    p_values['4_6'] = cal_p_value(nlls['4thCheby'], nlls['6thCheby'], 2)
    p_values['5_6'] = cal_p_value(nlls['5thCheby'], nlls['6thCheby'], 1)
    #print(p_values)

    # To avoid LaTex table shortening


    # Make a table and save to txt file
    out_path = f'tables/llr_table_pt{pt_low}_{pt_high}_y1.6_2.4.txt'
    with open(out_path, 'w') as out_file:
        #out_file.write('\\documentclass[10pt]{article}\n')
        #out_file.write('\\usepackage[usenames]{color} %used for font color\n')
        #out_file.write('\\usepackage{amssymb} %maths\n')
        #out_file.write('\\usepackage{amsmath} %maths\n')
        #out_file.write('\\usepackage[utf8]{inputenc} %useful to type directly diacritic characters\n')
        #out_file.write('\\begin{document}\n')
        #out_file.write('\n\n\n')
        out_file.write('\\begin{table}[h]\n')
        out_file.write('\\caption{$p_T=' + f'{pt_low}-{pt_high}' +'$ GeV, $y=1.6-2.4 $}\n')
        out_file.write('\\centering\n')
        out_file.write('\\begin{tabular}{c|c|c|c|c|c|c|c}\n')
        out_file.write('\\hline\n')
        out_file.write(' N & NLL & p(H0: & p(H0: & p(H0: & p(H0: & p(H0: & p(H0: \\\\ \n')
        out_file.write('   &     &  N=0) &  N=1) & N=2)  & N=3)  & N=4)  & N=5)  \\\\ \n')
        out_file.write('\\hline\n')
        out_file.write('\\hline\n')
        out_file.write(f' 0 & {nlls["Expo"]} & & & & & & \\\\ \n')
        out_file.write(f' 1 & {nlls["1stCheby"]} & {p_values["0_1"]} \\% & & & &    \\\\ \n')
        out_file.write(f' 2 & {nlls["2ndCheby"]} & {p_values["0_2"]} \\% & {p_values["1_2"]} \\% & & & & \\\\ \n')
        out_file.write(f' 3 & {nlls["3rdCheby"]} & {p_values["0_3"]} \\% & {p_values["1_3"]} \\% & {p_values["2_3"]} \\% & & &  \\\\ \n')
        out_file.write(f' 4 & {nlls["4thCheby"]} & {p_values["0_4"]} \\% & {p_values["1_4"]} \\% & {p_values["2_4"]} \\% & {p_values["3_4"]} \\% & & \\\\ \n')
        out_file.write(f' 5 & {nlls["5thCheby"]} & {p_values["0_5"]} \\% & {p_values["1_5"]} \\% & {p_values["2_5"]} \\% & {p_values["3_5"]} \\% & {p_values["4_5"]} \\% & \\\\ \n')
        out_file.write(f' 6 & {nlls["6thCheby"]} & {p_values["0_6"]} \\% & {p_values["1_6"]} \\% & {p_values["2_6"]} \\% & {p_values["3_6"]} \\% & {p_values["4_6"]} \\% & {p_values["5_6"]} \\% \\\\ \n')
        out_file.write('\\hline\n')
        out_file.write('\\end{tabular}\n')
        out_file.write('\\end{table}\n')
        #out_file.write('\n\n\n')
        #out_file.write('\\end{document}\n')
        



# Mid
# 빈 별로 계산해서 저장하기
for pt_low, pt_high in mid_pt:
    nlls = {}
    # Fill the nll values
    for bkg_name in bkg_names:
        in_file_name = f'Mass_Inclusive_pt{pt_low}-{pt_high}_y0.0-1.6_muPt0.0_{bkg_name}.txt'
        input_path = in_dir_path + '/' + in_file_name

        # Read NLL from input
        with open(input_path, 'r') as in_file:
            content = in_file.read()
        #print(content)
        nlls[f'{bkg_name}'] = content
    #print(nlls)

    # Calculate p-values
    # Call .C file
    p_values = {'0_0':'0', }
    
    # 1st -> vs 1 -> 0
    p_values['0_1'] = cal_p_value(nlls['Expo'], nlls['1stCheby'], 1)
    
    # 2nd -> vs 0, 1
    p_values['0_2'] = cal_p_value(nlls['Expo'], nlls['2ndCheby'], 2)
    p_values['1_2'] = cal_p_value(nlls['1stCheby'], nlls['2ndCheby'], 1)

    # 3rd -> vs 0, 1, 2
    p_values['0_3'] = cal_p_value(nlls['Expo'], nlls['3rdCheby'], 3)
    p_values['1_3'] = cal_p_value(nlls['1stCheby'], nlls['3rdCheby'], 2)
    p_values['2_3'] = cal_p_value(nlls['2ndCheby'], nlls['3rdCheby'], 1)
    
    # 4th -> vs 0, 1, 2, 3
    p_values['0_4'] = cal_p_value(nlls['Expo'], nlls['4thCheby'], 4)
    p_values['1_4'] = cal_p_value(nlls['1stCheby'], nlls['4thCheby'], 3)
    p_values['2_4'] = cal_p_value(nlls['2ndCheby'], nlls['4thCheby'], 2)
    p_values['3_4'] = cal_p_value(nlls['3rdCheby'], nlls['4thCheby'], 1)
    
    # 5th -> vs 0, 1, 2, 3. 4
    p_values['0_5'] = cal_p_value(nlls['Expo'], nlls['5thCheby'], 5)
    p_values['1_5'] = cal_p_value(nlls['1stCheby'], nlls['5thCheby'], 4)
    p_values['2_5'] = cal_p_value(nlls['2ndCheby'], nlls['5thCheby'], 3)
    p_values['3_5'] = cal_p_value(nlls['3rdCheby'], nlls['5thCheby'], 2)
    p_values['4_5'] = cal_p_value(nlls['4thCheby'], nlls['5thCheby'], 1)

    # 6th -> vs 0, 1, 2, 3. 4, 5
    p_values['0_6'] = cal_p_value(nlls['Expo'], nlls['6thCheby'], 6)
    p_values['1_6'] = cal_p_value(nlls['1stCheby'], nlls['6thCheby'], 5)
    p_values['2_6'] = cal_p_value(nlls['2ndCheby'], nlls['6thCheby'], 4)
    p_values['3_6'] = cal_p_value(nlls['3rdCheby'], nlls['6thCheby'], 3)
    p_values['4_6'] = cal_p_value(nlls['4thCheby'], nlls['6thCheby'], 2)
    p_values['5_6'] = cal_p_value(nlls['5thCheby'], nlls['6thCheby'], 1)
    #print(p_values)

    # Make a table and save to txt file
    out_path = f'tables/llr_table_pt{pt_low}_{pt_high}_y0.0_1.6.txt'
    with open(out_path, 'w') as out_file:
        #out_file.write('\\documentclass[10pt]{article}\n')
        #out_file.write('\\usepackage[usenames]{color} %used for font color\n')
        #out_file.write('\\usepackage{amssymb} %maths\n')
        #out_file.write('\\usepackage{amsmath} %maths\n')
        #out_file.write('\\usepackage[utf8]{inputenc} %useful to type directly diacritic characters\n')
        #out_file.write('\\begin{document}\n')
        #out_file.write('\n\n\n')
        out_file.write('\\begin{table}[h]\n')
        out_file.write('\\caption{$p_T=' + f'{pt_low}-{pt_high}' +'$ GeV, $y=0.0-1.6 $}\n')
        out_file.write('\\centering\n')
        out_file.write('\\begin{tabular}{c|c|c|c|c|c|c|c}\n')
        out_file.write('\\hline\n')
        out_file.write(' N & NLL & p(H0: &  p(H0:   &  p(H0:  &   p(H0: & p(H0: & p(H0: \\\\ \n')
        out_file.write('   &     &   N=0) &   N=1)  &   N=2)  &     N=3) &  N=4) & N=5)  \\\\ \n')
        out_file.write('\\hline\n')
        out_file.write('\\hline\n')
        out_file.write(f' 0 & {nlls["Expo"]} & & & & & & \\\\ \n')
        out_file.write(f' 1 & {nlls["1stCheby"]} & {p_values["0_1"]} \\% & & & &    \\\\ \n')
        out_file.write(f' 2 & {nlls["2ndCheby"]} & {p_values["0_2"]} \\% & {p_values["1_2"]} \\% & & & & \\\\ \n')
        out_file.write(f' 3 & {nlls["3rdCheby"]} & {p_values["0_3"]} \\% & {p_values["1_3"]} \\% & {p_values["2_3"]} \\% & & &  \\\\ \n')
        out_file.write(f' 4 & {nlls["4thCheby"]} & {p_values["0_4"]} \\% & {p_values["1_4"]} \\% & {p_values["2_4"]} \\% & {p_values["3_4"]} \\% & & \\\\ \n')
        out_file.write(f' 5 & {nlls["5thCheby"]} & {p_values["0_5"]} \\% & {p_values["1_5"]} \\% & {p_values["2_5"]} \\% & {p_values["3_5"]} \\% & {p_values["4_5"]} \\% & \\\\ \n')
        out_file.write(f' 6 & {nlls["6thCheby"]} & {p_values["0_6"]} \\% & {p_values["1_6"]} \\% & {p_values["2_6"]} \\% & {p_values["3_6"]} \\% & {p_values["4_6"]} \\% & {p_values["5_6"]} \\% \\\\ \n')
        out_file.write('\\hline\n')
        out_file.write('\\end{tabular}\n')
        out_file.write('\\end{table}\n')
        #out_file.write('\n\n\n')
        #out_file.write('\\end{document}\n')