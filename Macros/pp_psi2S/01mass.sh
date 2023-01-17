#!/bin/bash

root -l -b <<EOF
.L MassFit_FixPar_Data.C 
.q
EOF

for pt in  '3,6.5' 
do
	root -l -b -q 'MassFit_FixPar_Data_LowPt.C('$pt',1.6,2.4,0,200)'
done
for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
do
	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,1.6,0,200)'
done

for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,200'
do
    root -l -b -q 'MassFit_FixPar_Data.C(6.5,50,0,1.6,'$cent')'
done
for cent in '0,40' '40,80' '80,200'
do
    root -l -b -q 'MassFit_FixPar_Data_LowPt.C(3,50,1.6,2.4,'$cent')'
done


#for pt in  '3,6.5' 
#do
#	root -l -b -q 'MassFit_FixPar_Data_LowPt.C('$pt',1.6,2.4,0,180)'
#done
#for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
#do
#	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,2.4,0,180)'
#done
#
#for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
#do
#    root -l -b -q 'MassFit_FixPar_Data.C(6.5,50,0,2.4,'$cent')'
#done
#for cent in '0,40' '40,80' '80,180'
#do
#    root -l -b -q 'MassFit_FixPar_Data.C(3,50,1.6,2.4,'$cent')'
#done
