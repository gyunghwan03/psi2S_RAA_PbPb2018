#!/bin/bash

root -l -b <<EOF
.L Final2DFit.C 
.q
EOF



for pt in  '3.5,6.5' '6.5,9' '9,12' '12,40' '3.5,40'
do
	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4,0,180)'
done

for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,40' '6.5,40' #'20,40'
do
	root -l -b -q 'Final2DFit.C('$pt',0,1.6,0,180)'
done

for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
do
   root -l -b -q 'Final2DFit.C(6.5,40,0,1.6,'$cent')'
done

for cent in  '0,20' '20,60' '60,100' '100,180'
do
   root -l -b -q 'Final2DFit.C(3.5,40,1.6,2.4,'$cent')'
done
#

#for pt in  '3,6.5' 
#do
#	root -l -b -q 'MassFit_FixPar_Data_LowPt.C('$pt',1.6,2.4,0,180)'
#done
#for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,40'
#do
#	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,2.4,0,180)'
#done
#
#for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
#do
#    root -l -b -q 'MassFit_FixPar_Data.C(6.5,40,0,2.4,'$cent')'
#done
#for cent in '0,40' '40,80' '80,180'
#do
#    root -l -b -q 'MassFit_FixPar_Data.C(3,40,1.6,2.4,'$cent')'
#done
