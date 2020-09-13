# runjam - Cooper-Frye sampler & CLI for JAM1 cascades

An efficient hadron sampler based on the Cooper-Frye formula. This
also provides a command-line interface for JAM cascades useful for
event-by-event calculations of heavy-ion collisions.

- The Monte-Carlo sampling by the combination of the inversion method
  & rejection method.  Isotropic components are sampled by the
  inversion method based on the pre-calculated interpolation table of
  one-dimentional integrations.  Anisotropy is taken into account by
  the rejection sampling.
- The isotropic part is reduced to one-dimensional integration by
  considering the symmery in the local rest frame.  High-precision
  interpolation of integral tables by the combination of the Lagrange
  interpolation of the Chebyshef nodes and the cubic splines.
- Baryon density and diffusion are currently not considered.

## Compile

### Prerequisites

This package requires the library `libjam` from JAM version 1.820 or
before.  The source code of JAM can be obtained from the following URL:

  http://www.aiu.ac.jp/~ynara/jam/

In the JAM package, the parameter `mxv` defined in `src/jam1.inc`
should be rewritten to 200000 before the compile of JAM.
For illustration, the JAM library can be installed by the following commands.

```bash
$ wget http://www.aiu.ac.jp/~ynara/jam/old/jam-1.820.tar.bz2
$ tar xf jam-1.820.tar.bz
$ cd jam-1.820
$ EDIT src/jam1.inc
$ export F77=gfortran
$ ./configure --prefix=$HOME/opt/jam/1.820
$ make -j
$ make install
```

### Compile

`runjam.exe` can be compiled by the following commands.

```bash
$ git clone https://github.com/idthic/runjam.git
$ cd runjam
$ ./configure --prefix=$HOME/opt/idt --with-jam=$HOME/opt/jam/1.820
$ make
$ make install
```

### License

This programs is provided under [GPLv2](LICENSE).

```
runjam (idt) - Sample hadrons by Cooper-Frye formula / Run JAM cascade.
Copyright (C) 2011-2020, Koichi Murase @akinomyoga

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
```

## Usage

```console
$ runjam.exe --help
usage: runjam [SUBCOMMAND] [OPTIONS|VAR=VALUE]

SUBCOMMAND
  cascade (default)
  decay
  sample

OPTIONS and VARIABLES
  Variables can be specified by the environment variables or the command-line
  argument of the form `VAR=VALUE'.  In command-line arguments, the prefix
  `runjam_' can be omitted.

      runjam_mode=SUBCOMMAND
  -n, runjam_nevent=INT [1]              number of events to process
  -s, runjam_seed=INT [18371]            seed for random numbers in runjam
      runjam_jamseed=INT [runjam_seed]   seed for random numbers in JAM
  -t, runjam_oversampling_factor=NUM [1] number of test particles
  -w, runjam_switch_weak_decay=BOOL [false] enable weak decays
      runjam_phi_decays=BOOL [true]

 Output options
  -o,        runjam_output_directory=DIR [out]   directory of output files
  --fphase,  runjam_fname_phdat=FILE []          output filename
  --fphase0, runjam_fname_phdat0=FILE []         output filename
             runjam_output_phdat=BOOL [true]     output phasespace data
             runjam_output_phdat0=BOOL [true]    output phasespace0 data
             runjam_output_phbin=BOOL [false]    output binary phasespace
             runjam_output_phbin0=BOOL [false]   output binary phasespace0
             runjam_output_phdat_indexed=BOOL [false]
             runjam_output_phdat0_indexed=BOOL [false]
             runjam_output_index_start=INT [0]
  -d INT
    0        Disable all output format
    1        Enable only 'phdat' and 'phdat0'
    2        Enable only 'phbin' and 'phbin0'
    3        Enable only 'phdat_indexed' and 'phdat0_indexed'

 Initialization options
  -i ICSPEC       specify initial condition
    c0lrf:FILE    sample particles from the hypersurface data from
                  rfh c0lrf format "hypersurface_v1.txt"
    hydrojet:DIR  sample particles using the hypersurface data from
                  hydrojet (DIR/freezeout.dat, DIR/position.dat)
    phase:FILE    load particle lists from the text format "phasespace.dat".
    phase1:FILE   load a particle list from FILE. This performs an additional
                  check to require that FILE contains only a single event.
    phbin:FILE    load particle lists from the binary format "ph000k.bin".
    psample:FILE  read a particle list from the file in the format of
                  runjam "particlesample_pos.dat"

 Resonance list
  -r, --resodata, runjam_resodata=FILE   resonance data
  -p, runjam_eospce=INT [6]              eospce
  -k, runjam_kintmp=INT [5]              freezeout temperature type

 Options for viscous sampler
  --switching-temperature, runjam_switching_temperature=TEMP [155]
                    an advice to switching temperature in MeV
  runjam_turnsOffViscousEffect=INT

 Options for hydrojet hypersurface
  --hydrojet-bfree, hydrojet_baryonfree=INT [1]   baryonfree
  --hydrojet-dir,   hydrojet_directory=DIR [test] directory of freezeout.dat
  --hydrojet-dt,    hydrojet_deltat=NUM [0.3]     delta tau
  --hydrojet-dx,    hydrojet_deltax=NUM [0.3]     delta x
  --hydrojet-dy,    hydrojet_deltay=NUM [0.3]     delta y
  --hydrojet-dh,    hydrojet_deltah=NUM [0.3]     delta eta
  hydrojet_reverse_particles=BOOL [false]         (debug) perform z-reflection
  hydrojet_shuffle_particles=BOOL [false]         (debug) shuffle particles
  hydrojet_rotate_freezeout=BOOL [false]          (debug) rotate freezeout data

 Other options
  --help          show this help
  --version       show version information

EXAMPLE

$ ./runjam cascade -s 12345 -o jam -i c0lrf:hypersurface_v1.txt
$ ./runjam decay   -s 12345 -o jam -i phase:phasespace0.in
$ ./runjam sample  -s 12345 -o jam -i c0lrf:hypersurface_v1.txt -n 1000 -d 2
```

### Input/output data format of phasespace data (`phdat`)

This contains the sequence of event data in text format.

Each event data starts with a line containing two fields.  The first field
contains the number of particles in the event.  When this number is `-999`, it
means that this is the end of the file and that no more events are available.
The second column contains the oversampling factor, i.e., the factor for the
number of particles in the test-particle method for the Boltzmann equation.
This is normally identity which corresponds to the real number of hadrons.

Then as many lines as the number of the particles are followed.  Each line
carries one particle information with the following format:

```
ks    kf    px    py    pz    m    x    y    z    t
```


`ks` is `1` for stable qparticle and `2` for unstable particle.  `kf` is PDG
Monte-Carlo code.  `px py pz m` are the momentum and the mass in the unit
[GeV].  `x y z t` are the spacetime point of the last interaction vertex or the
formation vertex.

### Input/output data format of phasespace binary (`phbin`)

This is the binary format.  This contains the sequence of event frames.  Each
event frame starts with an eight-byte header followed by a sequence of particle
entries.  The structure is summarized in the following pseudo code.

```c
struct PARTICLE {
  // 1: stable, 2: unstable particle
  UINT32 ks;

  // PDG Monte-Carlo code
  UINT32 kf;
  
  // Momentum and mass in [GeV]
  SINGLE px, py, pz, m;
  
  // The formation point or the last interaction point in [fm]
  SINGLE x, y, z, t;
};

struct PHBIN_EVENT_FRAME {
  FOURCC magic; // 4-bytes. Fixed to be a string "EvPh"
  UINT32 numberOfParticles;
  PARTICLE data[numberOfParticles];
};

struct PHBIN {
  PHBIN_EVENT_FRAME fileContents[numberOfEvents];
};
```

### Input data format of particle sample data (`psamp`)

This data format contains just a sequence of particle information without any
header information.  Each particle is represented in a line of the following
format:

```
px    py    pz    e    m    ireso    tau    x    y    eta
```

`e` is the energy in the unit [GeV].  `ireso` is the particle ID in the
resonance list described below.  `tau` and `eta` are spacetime coordinates
specified by the Bjorken proper time and the spacetime rapidity.


### Input data format of rfh hypersurface data (`c0lrf`)

The first line contains two fields.  The first field indicating the file format
is always `c0lrf`.  The second field indicates the grid coordinates used in the
hydrodynamic calculation.  This needs to be `taueta` or `taueta-tilde`.
`taueta` is the natural tensor basis defined by the differentiation with
respect to the coordinate variables.  `taueta-tilde` corresponds to the
orthonormalized tensor basis, i.e., Galilei frame comoving with the fixed
spatial coodinates.  The subsequent lines are the information of hypersurface
elements.  Each line corresponds to a hypersurface element of the following
format:

```
dir ds0 ds1 ds2 ds3 x0 x1 x2 x3 dx0 dx1 dx2 dx3 T u0 u1 u2 u3 e p pi11 pi12 pi13 pi22 pi23 pi33
```

In this format, the Lorentz indices (0, 1, 2, 3) correspond to (tau, eta, x,
y), respectively.  `dir` is the direction of the hypersurface.  `ds0 ds1 ds2
ds3` contains the hypersurface element dσ<sub>μ</sub> (direction and
size). `x0 x1 x2 x3` contains the center position of the hypersurface
element. `dx0 dx1 dx2 dx3` contains the extent of the hypersurface element
along the grid coordinates. `T e p` are temperature, energy density, and
pressure, respectively. `u0 u1 u2 u3` is the flow velocity at the point in the
grid coordinates.  `pi11 pi12 pi13 pi22 pi23 pi33` are the shear stress tensor
in the local rest frame (singly boosted from the grid cordinates).

### Input data format of hydrojet hypersurface data (`hydrojet`)

The hypersurface data from hydrojet is saved in two separate files
`freezeout.dat` and `position.dat`.  Both files contain the same numbers of the
lines and are just a sequence of lines each corresponding to a hypersurface
element.  The line in `freezeout.dat` has the following format:

```
isbulk ds dx dy T eta mu vxp vyp yvp
```

`ds dx dy` contain `ds0 dx dy`, respectively if `isbulk` is 1, or otherwise
contain `ds1 dx dy`, respectively.  All the tensor basis are in `taueta-tilde`
coodinates. `T mu` are temperature and chemical potential, respectively.  Only
`mu` = 0.0 is supported by `runjam`.  `eta` is the spacetime rapidity of the
position of the hypersurface element.  `vxp` and `vyp` are the flow velocity
(dx/dt, dy/dt).  `yvp` is the rapidity of the flow velocity in the z direction.
The line in `position.dat` has the following format:

```
tau x y eta
```

`tau x y eta` is the position of the hypersurface element.  Note that these
fields in both files actually do not have to form a line for each hypersurface
element, i.e., each field can be placed in an independent line as in the output
of the original `hydrojet`.


### Change resonance list

The resonance list can be specified by the option `-r, --resodata,
runjam_resodata=FILE`.  This option selects the file that contains the
information of particle species to be sampled.  The default resonance
list file is [`data/ResonanceJam.dat`](data/ResonanceJam.dat).  When
you run cascades the default should be used.

Here the file format of the resonance list file is described.  Empty
lines and lines starting with `#` are ignored.  Each line contains the
information of a particle species with the following format:

- Column 1: Mass
- Column 2: Degeneracy
- Column 3: Effective degeneracy (average number of pions after decays)
- Column 4: Chemical potential
- Column 5: Statistics (1: boson, 2: fermion)
- Column 6: Is antiparticle
- Column 7: KEY.  The name used to determine the filename for the resonance.
- Column 8: Description
- Column 9+: corresponding PDG Monte-Carlo codes.  Note that there can be multiple codes corresponding to the isospin degenracy.

When this option is not specified, the resonance list file is determined
based on the other options `--hydrojet-pce, hydrojet_eospce=INT` and `--hydrojet-ftemp, hydrojet_kintmp=INT`.
The default is `eospce=6` and `kintmp=5` so that `ResonanceJam.dat` is used.

- `eospce=0`: `ResonancePCE.dat` will be used as the list of particles. 21 resonances are contained.
- `eospce=1` The same as `eospce=0` but with the chemical potentials (PCE)
  - `kintmp=1`: `ResonancePCE.T080.dat`. Freezeout temperature 80 MeV
  - `kintmp=2`: `ResonancePCE.T100.dat`. Freezeout temperature 100 MeV
  - `kintmp=3`: `ResonancePCE.T120.dat`. Freezeout temperature 120 MeV
  - `kintmp=4`: `ResonancePCE.T140.dat`. Freezeout temperature 140 MeV
  - `kintmp=5`: `ResonancePCE.T160.dat`. Freezeout temperature 160 MeV
- `eospce=4`: `ResonanceEosqJam.dat`. EOS-Q
- `eospce=5`: `ResonancePCE.New.dat`. `eospce=1` with updated masses
- `eospce=10`: `ResonanceCharged.Massless.dat`.
- `eospce=11`: `ResonanceCharged.dat`.
- `eospce=12` The same as `eospce=11` but with the chemical potentials (PCE)
  - `kintmp=1`: `ResonanceCharged.T080.dat`. Freezeout temperature 80 MeV
  - `kintmp=2`: `ResonanceCharged.T100.dat`. Freezeout temperature 100 MeV
  - `kintmp=3`: `ResonanceCharged.T120.dat`. Freezeout temperature 120 MeV
  - `kintmp=4`: `ResonanceCharged.T140.dat`. Freezeout temperature 140 MeV
  - `kintmp=5`: `ResonanceCharged.T160.dat`. Freezeout temperature 160 MeV
- `eospce=13`: `ResonancePhi.dat`. phi and J/psi mesons
- `eospce=14`: `ResonancePhi.T100.dat`. phi and J/psi mesons with PCE T = 100 MeV.
- Otherwise: `ResonanceJam.dat`
