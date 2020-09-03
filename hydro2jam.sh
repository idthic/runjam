#!/bin/bash

# Basic options
hydro2jam_cascade_mode=cascade
hydro2jam_nevent=1
hydro2jam_ievent_begin=0
hydro2jam_output_directory=jam
hydro2jam_seed=18371
hydro2jam_jamseed=

# Options for output
hydro2jam_phasespace_enabled=1
hydro2jam_phasespace_fname=phasespace.dat
hydro2jam_phasespace_fname0=phasespace0.dat
hydro2jam_output_phbin=false
hydro2jam_output_phbin0=false

# Options for Cooper-Frye formula
hydro2jam_turnsOffViscousEffect=0
hydro2jam_oversampling_factor=1.0
hydro2jam_switching_temperature=-1
hydro2jam_resodata=ResonanceJam.dat

# Note: currently this should be 0
hydro2jam_switch_weak_decay=0

hydrojet_directory=test
hydrojet_eospce=6
hydrojet_kintmp=5
hydrojet_baryonfree=1
hydrojet_deltat=0.3
hydrojet_deltah=0.3
hydrojet_deltax=0.3
hydrojet_deltay=0.3

# Options for JAM
hydro2jam_phi_decays=1
hydro2jam_decay_only=0

export "${!hydro2jam_@}" "${!hydrojet_@}"
./hydro2jam.exe
