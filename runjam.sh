#!/bin/bash

# Basic options
runjam_cascade_mode=cascade
runjam_nevent=1
runjam_ievent_begin=0
runjam_output_directory=jam
runjam_seed=18371
runjam_jamseed=

# Options for output
runjam_phasespace_enabled=1
runjam_phasespace_fname=phasespace.dat
runjam_phasespace_fname0=phasespace0.dat
runjam_output_phbin=false
runjam_output_phbin0=false

# Options for Cooper-Frye formula
runjam_turnsOffViscousEffect=0
runjam_oversampling_factor=1.0
runjam_switching_temperature=-1
runjam_resodata=ResonanceJam.dat

# Note: currently this should be 0
runjam_switch_weak_decay=0

hydrojet_directory=test
hydrojet_eospce=6
hydrojet_kintmp=5
hydrojet_baryonfree=1
hydrojet_deltat=0.3
hydrojet_deltah=0.3
hydrojet_deltax=0.3
hydrojet_deltay=0.3

# Options for JAM
runjam_phi_decays=1
runjam_decay_only=0

export "${!runjam_@}" "${!hydrojet_@}"
./runjam.exe
