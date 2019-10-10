#!/bin/bash

# Basic options
hydro2jam_seed=18371
hydro2jam_jamseed=
hydro2jam_output_directory=jam
hydro2jam_nevent=1
hydro2jam_ievent_begin=0

# Options for output
hydro2jam_phasespace_enabled=1
hydro2jam_phasespace_fname=phasespace.dat
hydro2jam_phasespace_fname0=phasespace0.dat

# Options for Cooper-Frye formula
hydro2jam_turnsOffViscousEffect=0

hydrojet_directory=test
hydrojet_eospce=6
hydrojet_kintmp=5
hydrojet_resodata=dict/ResonanceJam.dat

hydro2jam_deltat=0.3
hydro2jam_deltah=0.3
hydro2jam_deltax=0.3
hydro2jam_deltay=0.3

# Options for JAM
hydro2jam_phi_decays=1

export "${!hydro2jam_@}" "${!hydrojet_@}"
./hydro2jam.exe
