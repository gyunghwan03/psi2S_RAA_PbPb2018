import subprocess
from concurrent.futures import ThreadPoolExecutor
commands = [
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 2)'", 

"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 2)'", 


"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5,  40, 1.6, 2.4, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5,  40, 1.6, 2.4, 2)'", 
]

# 병렬 실행 함수
def run_command(command):
    subprocess.call(command, shell=True)

# ThreadPoolExecutor를 사용하여 병렬 실행
with ThreadPoolExecutor(max_workers=1) as executor:  # 최대 4개의 작업을 동시에 실행
    executor.map(run_command, commands)