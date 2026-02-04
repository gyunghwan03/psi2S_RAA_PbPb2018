#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_3_3_patch1
eval `scramv1 runtime -sh`
cd /data/hwan/psi2S_RAA_PbPb2018/Eff_Acc
root -l -b -q 'decay_length_OniaTree_v2_pp_2S.C(1, 1, true, 0, true, 3.0, 40.0, 0.0, 2.4)'
