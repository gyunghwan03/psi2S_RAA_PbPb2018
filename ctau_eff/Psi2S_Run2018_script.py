import subprocess

subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 3, 6.5, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 4, 6.5, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 12, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 12, 50, 1.6, 2.4)'", shell=True)

subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 40, 3, 50, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(20, 80, 3, 50, 1.6, 2.4)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(80, 180, 3, 50, 1.6, 2.4)'", shell=True)

subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 9, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 9, 12, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 12, 15, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 15, 20, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 180, 20, 50, 0, 1.6)'", shell=True)

subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(0, 20, 6.5, 50, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(20, 40, 6.5, 50, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(40, 60, 6.5, 50, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(60, 80, 6.5, 50, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(80, 100, 6.5, 50, 0, 1.6)'", shell=True)
subprocess.call("root -l -b -q psuedo_proper_decay_length.C'(100, 180, 6.5, 50, 0, 1.6)'", shell=True)