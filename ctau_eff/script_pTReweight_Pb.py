import subprocess
from concurrent.futures import ThreadPoolExecutor

# 실행할 명령어 리스트
## 1S ##
commands_1S = [
    # PRMC
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 6.5, 7.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 7.5, 8.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 8.5, 9.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 9.5,  11, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 11,   13, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 13,   15, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 15, 17.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 17.5, 20, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 20, 22.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 22.5, 25, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 25, 27.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 30,   40, 0, 1.6, 0, 2.6, 3.5, 0)'",

    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 3.5, 4.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 4.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 6.5, 8.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 8.5,  12, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  12,  15, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  15,  20, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  20,  40, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    # NPMC
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 6.5, 7.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 7.5, 8.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 8.5, 9.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 9.5,  11, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 11,   13, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 13,   15, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 15, 17.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 17.5, 20, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 20, 22.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 22.5, 25, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 25, 27.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 30,   40, 0, 1.6, 0, 2.6, 3.5, 1)'",

    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 3.5, 4.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 4.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 6.5, 8.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180, 8.5,  12, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  12,  15, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  15,  20, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length_1S.C'(0, 180,  20,  40, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
]

## 2S ##
commands_2S = [
    # PRMC
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 7.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 7.5, 8.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 8.5, 9.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 9.5,  11, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 11,   13, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 13,   15, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 15, 17.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 17.5, 20, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 20, 22.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 22.5, 25, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 25, 27.5, 0, 1.6, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 30,   40, 0, 1.6, 0, 2.6, 3.5, 0)'",

    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 3.5, 4.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 4.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 8.5, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 8.5,  12, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  12,  15, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  15,  20, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  20,  40, 1.6, 2.4, 0, 2.6, 3.5, 0)'",
    # NPMC
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 7.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 7.5, 8.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 8.5, 9.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 9.5,  11, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 11,   13, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 13,   15, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 15, 17.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 17.5, 20, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 20, 22.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 22.5, 25, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 25, 27.5, 0, 1.6, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 30,   40, 0, 1.6, 0, 2.6, 3.5, 1)'",

    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 3.5, 4.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 4.5, 6.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 6.5, 8.5, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180, 8.5,  12, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  12,  15, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  15,  20, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
    "root -l -b -q psuedo_proper_decay_length.C'(0, 180,  20,  40, 1.6, 2.4, 0, 2.6, 3.5, 1)'",
]


# 병렬 실행 함수
def run_command(command):
    subprocess.call(command, shell=True)

# ThreadPoolExecutor를 사용하여 병렬 실행
with ThreadPoolExecutor(max_workers=12) as executor:  # 최대 4개의 작업을 동시에 실행
    executor.map(run_command, commands_1S + commands_2S)