#!/bin/bash

#root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C( ptLow, ptHigh, yLow, yHigh, cLow, cHigh, isTnP, isPtWeight)'


# Prompt
# Forward Rapidity
# # Pt
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 3.5, 1.6, 2.4, 0, 181, true, true)' &
# exit ### TEST EXIT ###

root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 3.5, 5.0, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 5.0, 6.5, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 6.5, 12.0, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 12.0, 50., 1.6, 2.4, 0, 181, true, true)'

# # Centrality
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 0, 20, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 20, 40, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 40, 60, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 60, 80, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 80, 100, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 1.6, 2.4, 100, 181, true, true)'

# Mid Rapidity
# # Pt
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 6.5, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 6.5, 9.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 9.0, 12.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 12.0, 15.0, 0.0, 1.6, 0, 181, true, true)' & 
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 15.0, 20.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 20.0, 25.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 25.0, 30.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 30.0, 50.0, 0.0, 1.6, 0, 181, true, true)'

# # Centrality
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 0, 20, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 20, 40, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 40, 60, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 60, 80, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 80, 100, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(true, 0.0, 50., 0.0, 1.6, 100, 181, true, true)'

# Non-prompt
# # Pt
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 3.5, 1.6, 2.4, 0, 181, true, true)' &

root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 3.5, 5.0, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 5.0, 6.5, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 6.5, 12.0, 1.6, 2.4, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 12.0, 50., 1.6, 2.4, 0, 181, true, true)'

# # Centrality
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 0, 20, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 20, 40, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 40, 60, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 60, 80, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 80, 100, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 1.6, 2.4, 100, 181, true, true)' 

# Mid Rapidity
# # Pt
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 6.5, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 6.5, 9.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 9.0, 12.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 12.0, 15.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 15.0, 20.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 20.0, 25.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 25.0, 30.0, 0.0, 1.6, 0, 181, true, true)' &
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 30.0, 50.0, 0.0, 1.6, 0, 181, true, true)'

# # Centrality
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 0, 20, true, true)'
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 20, 40, true, true)'
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 40, 60, true, true)'
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 60, 80, true, true)'
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 80, 100, true, true)'
root -l -b -q 'getEfficiency_psi_pbpb_SYSTNP.C(false, 0.0, 50., 0.0, 1.6, 100, 181, true, true)'
