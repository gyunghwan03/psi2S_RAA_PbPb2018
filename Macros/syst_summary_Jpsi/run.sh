#!/bin/bash

root -l -b -q syst_HFdown.C
root -l -b -q syst_HFup.C
root -l -b -q merge_HF.C
root -l -b -q syst_bkgPDF.C
root -l -b -q syst_ctauBkg.C
root -l -b -q syst_ctauErr.C
root -l -b -q syst_ctauRes.C
root -l -b -q syst_ctauTrue.C
root -l -b -q syst_bFrac.C
root -l -b -q syst_sigPAR.C
root -l -b -q syst_sigPDF.C
root -l -b -q syst_eff.C
root -l -b -q syst_acc.C
root -l -b -q compute_total_syst.C
root -l -b -q draw_bFrac.C
root -l -b -q draw_syst_ver2.C

