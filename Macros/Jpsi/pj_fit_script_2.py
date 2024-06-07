import subprocess as sub
import argparse

# Usage example
# python3 pj_fit_script.py --step=mass 

parser = argparse.ArgumentParser()
parser.add_argument('--step', required=True)
args = parser.parse_args()

#6.5 - 12
#12  - 50


#6.5 - 9
#9 - 12
#12 - 15
#15 - 20
#6.5 - 50
#6.5-7.5, 7.5-8.5, 8.5-9.5, 9.5-11, 11-13, 13-15 , 15-17.5, 17.5-20, 20-25, 25-30, 30-50

#1st value
#pt = [6.5, 50]
#y = [0, 2.4]

#pt = [7.5 , 8.5]
#pt = [8.5 , 9.5]
#pt = [9.5 , 11.]

#####################
#y = [1.6, 2.4]
#c = [0, 180]
#pt = [3.5 , 6.5]
#pt = [6.5 , 9.]
#pt = [9. , 12.]
#pt = [12. , 50.]
######################

#y = [1.6, 2.4]
#pt = [3.5 , 50.]
#c = [0, 20]
#c = [20, 40]
#c = [40, 60]
#c = [60, 80]
#c = [80, 100]
#c = [100, 180]

y = [0., 1.6]
pt = [6.5 , 50.]
#c = [0, 20]
#c = [20, 40]
#c = [40, 60]
c = [60, 80]
#c = [80, 100]
#c = [100, 180]

######################
#y = [0, 1.6]
#c = [0, 180]
#pt = [6.5 , 9.]
#pt = [9. , 12.]
#pt = [12. , 15.]
#pt = [15. , 20.]
#pt = [20. , 25.]
#pt = [25. , 50.]
######################

steps = {}
steps['mc'] = f"mc_MassFit_HighpT.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
steps['mass'] = f"MassFit_FixPar_Data.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
steps['err'] = f"CtauErr.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
steps['res'] = f"CtauRes.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
steps['true'] = f"CtauTrue.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
#steps['bkg'] = f"CtauBkg.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},c[1]})\'"
steps['bkg'] = f"CtauBkg_2exp.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
#if pt[0] < 6.5: # Use LowPt.C code for low pT region
#    steps['bkg'] = f"CtauBkg_LowPt.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
#else:
#    steps['bkg'] = f"CtauBkg.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"
steps['fin'] = f"Final2DFit.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{c[0]},{c[1]})\'"


fit_command = steps[args.step]
sub.call('root -l ' + fit_command, shell=True)
#sub.call('root -l -b -q ' + fit_command, shell=True)
