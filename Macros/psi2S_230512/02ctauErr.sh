#!/bin/bash

root -l -b <<EOF
.L CtauErr.C
.q
EOF

for y in '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4'
do
	root -l -b -q 'CtauErr.C(6.5,50,'$y',0,180)'
done

#for pt in  '3,6.5'
#do
#   root -l -b -q 'CtauErr.C('$pt',1.6,2.4,0,180)'
#done
#for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
#for pt in '20,30' '30,50'
#do
#   root -l -b -q 'CtauErr.C('$pt',0,1.6,0,180)'
#done
#for pt in '3,4.5' '4.5,6.5' '6.5,12' '12,50'
#do
#    root -l -b -q 'CtauErr.C('$pt',1.6,2.4,0,180)'
#done

#for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
#do
#    root -l -b -q 'CtauErr.C(6.5,50,0,1.6,'$cent')'
#done
#for cent in '0,40' '40,80' '80,180'
#do
#    root -l -b -q 'CtauErr.C(3,50,1.6,2.4,'$cent')'
#done

#for pt in  '3,6.5'
#do
#   root -l -b -q 'CtauErr.C('$pt',1.6,2.4,0,180)'
#done
#for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
#do
#   root -l -b -q 'CtauErr.C('$pt',0,2.4,0,180)'
#done
#
#for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
#do
#    root -l -b -q 'CtauErr.C(6.5,50,0,2.4,'$cent')'
#done
#for cent in '0,40' '40,80' '80,180'
#do
#    root -l -b -q 'CtauErr.C(3,50,1.6,2.4,'$cent')'
#done