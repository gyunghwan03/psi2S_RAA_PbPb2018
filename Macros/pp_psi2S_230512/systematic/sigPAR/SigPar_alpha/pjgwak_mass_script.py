import subprocess


subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(9, 12, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(12, 50, 1.6, 2.4)'", shell=True)

subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(3.5, 50, 1.6, 2.4)'", shell=True)


subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6)'", shell=True)

subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(25, 50, 0, 1.6)'", shell=True)

subprocess.call("root -l -b -q MassFit_FixPar_Data.C'(6.5, 50, 0, 1.6)'", shell=True)


# ======================================================================= #


subprocess.call("root -l -b -q Final2DFit.C'(3.5, 6.5, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(6.5, 9, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(9, 12, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(12, 50, 1.6, 2.4)'", shell=True)

subprocess.call("root -l -b -q Final2DFit.C'(3.5, 50, 1.6, 2.4)'", shell=True)


subprocess.call("root -l -b -q Final2DFit.C'(6.5, 9, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(9, 12, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(12, 15, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(15, 20, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(20, 25, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(25, 50, 0, 1.6)'", shell=True)

subprocess.call("root -l -b -q Final2DFit.C'(6.5, 50, 0, 1.6)'", shell=True)