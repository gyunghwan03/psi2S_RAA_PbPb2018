#!/bin/bash
root -l -b -q 'draw_Bfraction_y0_1p6.C'
root -l -b -q 'draw_Bfraction_y1p6_2p4.C'
for sys in '0' '1'
do
root -l -b -q 'draw_Raa_psi2S_y0_1p6_pT.C('$sys')'
root -l -b -q 'draw_Raa_psi2S_y1p6_2p4_pT.C('$sys')'
root -l -b -q 'draw_Raa_psi2S_y0_1p6_Cent.C('$sys')'
root -l -b -q 'draw_Raa_psi2S_y1p6_2p4_Cent.C('$sys')'
root -l -b -q 'compare_pT_PR.C('$sys')'
root -l -b -q 'compare_pT_NP.C('$sys')'
root -l -b -q 'compare_Npart.C('$sys')'
done