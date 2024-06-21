import subprocess as sub

cents = [(0,20), (20,40), (40,60), (60,80), (80,100), (100,180)]

### rapidity 1.6 - 2.4 ###
y = [1.6, 2.4]

# Cent 0 - 90 %
pt = [3.5, 6.5]
# Fit: Expo ~ cheby6
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [6.5, 9]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)


# Cent dependency
pt = [3.5, 50]
for i in range(0, 7):
    for cent in cents:
        fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{cent[0]},{cent[1]},{i})\'"
        sub.call('root -l -b -q ' + fit_command, shell=True)



### rapidity 0 - 1.6 ###
y = [0, 1.6]


pt = [6.5, 9]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [9, 12]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [12, 15]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [15, 20]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [20, 25]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)

pt = [25, 50]
for i in range(0, 7):
    fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},0,180,{i})\'"
    sub.call('root -l -b -q ' + fit_command, shell=True)


# Cent dependency
pt = [6.5, 50]
for i in range(0, 7):
    for cent in cents:
        fit_command = f"MassFit_LLR.C\'({pt[0]},{pt[1]},{y[0]},{y[1]},{cent[0]},{cent[1]},{i})\'"
        sub.call('root -l -b -q ' + fit_command, shell=True)
