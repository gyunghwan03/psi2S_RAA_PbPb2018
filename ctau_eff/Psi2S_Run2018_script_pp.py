import subprocess

# NPMC
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(6.5, 9, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(9, 12, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)