import subprocess
from concurrent.futures import ThreadPoolExecutor

#### 1S ####
# PRMC
commands_1S_PRMC = [
    "root -l -b -q syst_1S_pp.C'(6.5, 9, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(9, 12, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(12, 15, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(15, 20, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(20, 25, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(25, 40, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(6.5, 40, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(3.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(6.5, 9, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(9, 12, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(12, 40, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q syst_1S_pp.C'(3.5, 40, 1.6, 2.4, 0, 2.6, 3.5, 0)'"
]

# NPMC
commands_1S_NPMC = [
    "root -l -b -q syst_1S_pp.C'(6.5, 9, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(9, 12, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(12, 15, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(15, 20, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(20, 25, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(25, 40, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(6.5, 40, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(3.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(6.5, 9, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(9, 12, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(12, 40, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q syst_1S_pp.C'(3.5, 40, 1.6, 2.4, 0, 2.6, 3.5, 1)'"
]

#### 2S ####
# PRMC
commands_2S_PRMC = [
    "root -l -b -q syst_2S_pp.C'(6.5, 9, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(9, 12, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(12, 15, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(15, 20, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(20, 25, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(25, 40, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(6.5, 40, 0, 1.6, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(3.5, 6.5, 1.6, 2.4, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(6.5, 9, 1.6, 2.4, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(9, 12, 1.6, 2.4, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(12, 40, 1.6, 2.4, 0, 3.3, 4.1, 0)'",
    "root -l -b -q syst_2S_pp.C'(3.5, 40, 1.6, 2.4, 0, 3.3, 4.1, 0)'"
]

## NPMC
commands_2S_NPMC = [
    "root -l -b -q syst_2S_pp.C'(6.5, 9, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(9, 12, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(12, 15, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(15, 20, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(20, 25, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(25, 40, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(6.5, 40, 0, 1.6, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(3.5, 6.5, 1.6, 2.4, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(6.5, 9, 1.6, 2.4, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(9, 12, 1.6, 2.4, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(12, 40, 1.6, 2.4, 0, 3.3, 4.1, 1)'",
    "root -l -b -q syst_2S_pp.C'(3.5, 40, 1.6, 2.4, 0, 3.3, 4.1, 1)'"
]

# 병렬 실행 함수
def run_command(command):
    subprocess.call(command, shell=True)

# ThreadPoolExecutor를 사용하여 병렬 실행
with ThreadPoolExecutor(max_workers=12) as executor:  # 최대 4개의 작업을 동시에 실행
    #executor.map(run_command, commands_1S_PRMC + commands_1S_NPMC + commands_2S_PRMC + commands_2S_NPMC)
    executor.map(run_command, commands_2S_PRMC + commands_2S_NPMC)
