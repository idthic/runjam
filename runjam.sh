#!/bin/bash

# Basic options
runjam_cascade_mode=cascade
runjam_nevent=1
runjam_ievent_begin=0
runjam_output_directory=jam
runjam_seed=18371
runjam_jamseed=

#------------------------------------------------------------------------------
# Model options

# JAM options
runjam_phi_decays=1
runjam_decay_only=0

# Note: currently this should be 0
runjam_switch_weak_decay=0

# Options for Cooper-Frye formula
runjam_turnsOffViscousEffect=0
runjam_oversampling_factor=1.0
runjam_switching_temperature=-1
runjam_resodata=ResonanceJam.dat

#------------------------------------------------------------------------------
# Output options

# Options for phasespace data output
runjam_output_phdat=true
runjam_output_phdat0=true
runjam_fname_phdat=$runjam_output_directory/phasespace.dat
runjam_fname_phdat0=$runjam_output_directory/phasespace0.dat

# Options for phasespace binary output
runjam_output_phbin=false
runjam_output_phbin0=false
runjam_fname_phbin=$runjam_output_directory/phasespace.bin
runjam_fname_phbin0=$runjam_output_directory/phasespace0.bin

# Options for single-event phasespace data output
runjam_output_phbin_indexed=false
runjam_output_phbin0_indexed=false
runjam_output_index_start=0

#------------------------------------------------------------------------------
# Input options

hydrojet_directory=test
hydrojet_eospce=6
hydrojet_kintmp=5
hydrojet_baryonfree=1
hydrojet_deltat=0.3
hydrojet_deltah=0.3
hydrojet_deltax=0.3
hydrojet_deltay=0.3

export "${!runjam_@}" "${!hydrojet_@}"
./runjam.exe
