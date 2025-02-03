import subprocess

subprocess.call("root -l -b -q MassFit_FixPar_Data_pp.C'(6.5, 9, 0, 1.6, 2)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data_pp.C'(9, 12, 0, 1.6, 2)'", shell=True)