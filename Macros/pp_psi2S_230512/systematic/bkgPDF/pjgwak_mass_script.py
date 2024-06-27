import subprocess


#subprocess.call("root -l MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4)'", shell=True)
#subprocess.call("root -l MassFit_FixPar_Data.C'(6.5, 9, 1.6, 2.4)'", shell=True)
subprocess.call("root -l MassFit_FixPar_Data.C'(12, 50, 1.6, 2.4)'", shell=True)
subprocess.call("root -l MassFit_FixPar_Data.C'(3.5, 50, 1.6, 2.4)'", shell=True)


subprocess.call("root -l MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6)'", shell=True)
subprocess.call("root -l MassFit_FixPar_Data.C'(9, 12, 0, 1.6)'", shell=True)
subprocess.call("root -l MassFit_FixPar_Data.C'(12, 15, 0, 1.6)'", shell=True)

subprocess.call("root -l MassFit_FixPar_Data.C'(6.5, 50, 0, 1.6)'", shell=True)


#======================================================================= 

#subprocess.call("root -l Final2DFit.C'(3.5, 6.5, 1.6, 2.4)'", shell=True)
#subprocess.call("root -l Final2DFit.C'(6.5, 9, 1.6, 2.4)'", shell=True)
subprocess.call("root -l Final2DFit.C'(12, 50, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q Final2DFit.C'(3.5, 50, 1.6, 2.4)'", shell=True)

subprocess.call("root -l -b Final2DFit.C'(6.5, 9, 0, 1.6)'", shell=True)
subprocess.call("root -l -b Final2DFit.C'(9, 12, 0, 1.6)'", shell=True)
subprocess.call("root -l -b Final2DFit.C'(12, 15, 0, 1.6)'", shell=True)

subprocess.call("root -l -b Final2DFit.C'(6.5, 50, 0, 1.6)'", shell=True)