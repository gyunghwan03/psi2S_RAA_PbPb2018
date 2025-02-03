import subprocess

subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 2)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 0, 180, 2)'", shell=True)