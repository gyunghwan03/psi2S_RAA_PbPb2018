import subprocess as sub


### rapidity 1.6 - 2.4 ###
y = [1.6, 2.4]

pt = [3.5, 6.5]
# Fit: Expo ~ cheby6
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [6.5, 9]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [3.5, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)



### rapidity 0 - 1.6 ###
y = [0, 1.6]

pt = [6.5, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [6.5, 9]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 15]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [15, 20]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)
#for i in range(0, 7):
#    fit_command = f"MassFit_LLR_mc_params.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
#    sub.call('root -l -b -q ' + fit_command, shell=True)
#
pt = [20, 25]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [25, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)




##### Fit specific degree of bkg function #####
### rapidity 1.6 - 2.4 ###
'''
y = [1.6, 2.4]

pt = [3.5, 6.5]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [6.5, 9]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 50]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [3.5, 50]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)



### rapidity 0 - 1.6 ###
y = [0, 1.6]

pt = [6.5, 50]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [6.5, 9]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 15]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [15, 20]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)
#for i in [0,5,6]:
#    fit_command = f"MassFit_LLR_mc_params.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
#    sub.call('root -l -b -q ' + fit_command, shell=True)
#
pt = [20, 25]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [25, 50]
for i in [0,5,6]:
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)
'''