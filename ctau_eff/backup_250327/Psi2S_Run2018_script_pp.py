import subprocess

# PRMC
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(6.5, 9, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(9, 12, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(12, 15, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(15, 20, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(20, 25, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(25, 40, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)

subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(6.5, 40, 0, 1.6, 0, 3.3, 4.1, 0)'", shell=True)

# NPMC
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(6.5, 9, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(9, 12, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(12, 15, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(15, 20, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(20, 25, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(25, 40, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)

subprocess.call("root -l -b -q psuedo_proper_decay_length_pp.C'(6.5, 40, 0, 1.6, 0, 3.3, 4.1, 1)'", shell=True)
