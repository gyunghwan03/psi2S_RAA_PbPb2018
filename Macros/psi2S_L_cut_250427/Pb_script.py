import subprocess
from concurrent.futures import ThreadPoolExecutor

commands_nominal = [
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 0, 180, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 0, 180, 2)'", 


"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 0, 20, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 20, 40, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 40, 60, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 60, 80, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 80, 100, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 100, 180, 1)'", 


"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 0, 20, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 20, 40, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 40, 60, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 60, 80, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 80, 100, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 100, 180, 2)'", 


"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4,  0, 180, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4,  0, 180, 2)'", 

"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 0, 20, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 20, 60, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 60, 100, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 100, 180, 1)'", 


"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 0, 20, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 20, 60, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 60, 100, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 100, 180, 2)'", 
]

commands_pTreweight = [
    "root -l -b -q MassFit_FixPar_Data.C'(6.5, 7.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(7.5, 8.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(8.5, 9.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(9.5,  11, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(11,   13, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(13,   15, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(15, 17.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(17.5, 20, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(20, 22.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(22.5, 25, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(25, 27.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(27.5, 30, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(30,   40, 0, 1.6, 0, 180, 1)'", 

    "root -l -b -q MassFit_FixPar_Data.C'(3.5, 4.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(4.5, 6.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(6.5, 8.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(8.5,  12, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(12,   15, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(15,   20, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q MassFit_FixPar_Data.C'(20,   40, 1.6, 2.4, 0, 180, 1)'",

    "root -l -b -q MassFit_FixPar_Data.C'(6.5, 7.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(7.5, 8.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(8.5, 9.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(9.5,  11, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(11,   13, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(13,   15, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(15, 17.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(17.5, 20, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(20, 22.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(22.5, 25, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(25, 27.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(27.5, 30, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q MassFit_FixPar_Data.C'(30,   40, 0, 1.6, 0, 180, 2)'", 

    "root -l -b -q MassFit_FixPar_Data.C'(3.5, 4.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(4.5, 6.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(6.5, 8.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(8.5,  12, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(12,   15, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(15,   20, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q MassFit_FixPar_Data.C'(20,   40, 1.6, 2.4, 0, 180, 2)'",
]



# 병렬 실행 함수
def run_command(command_pTreweight):
    subprocess.call(command_pTreweight, shell=True)

# ThreadPoolExecutor를 사용하여 병렬 실행
with ThreadPoolExecutor(max_workers=4) as executor:  # 최대 4개의 작업을 동시에 실행
    executor.map(run_command, commands_pTreweight)
