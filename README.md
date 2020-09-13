# runjam

## Compile

### Requirements

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

`runjam.exe` will be created with the following commands.

```bash
$ ./configure --prefix=$HOME/opt/idt --with-jam=$HOME/opt/jam/1.820
$ make
$ make install
```


## Usage

Please check the output of `runjam.exe --help`.

### Option `--resodata, runjam_resodata=FILE`

This option can be used to select the file that contains the list of sampled particles.
Empty lines and lines starting with `#` are ignored.

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
- `eospce=10`: `ResonanceCharged.dat`.
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

## Changes 0.2..0.3

- change option prefix: `hydro2jam_*` -> `runjam_*`
- reorganize mode options
  - rename option `runjam_cascade_mode` -> `runjam_mode`
  - remove option `runjam_decay_only`
- fix bugs that some options did not work because of option name typos:
  - `runjam_switch_weak_decay`
  - `--hydrojet-dt, hydrojet_deltat`
  - `--hydrojet-dx, hydrojet_deltax`
  - `--hydrojet-dy, hydrojet_deltay`
  - `--hydrojet-dh, hydrojet_deltah`
- split the option for multi-event phasespace files
  - remove option `runjam_phasespace_enabled`
  - new option `runjam_output_phdat`
  - new option `runjam_output_phdat0`
- rename options for multi-event phasespace filenames
  - rename option `runjam_phasespace_fname` -> `runjam_fname_phdat`
  - rename option `runjam_phasespace_fname0` -> `runjam_fname_phdat0`
- add options for phasespace binary filenames
  - new option `runjam_fname_phbin`
  - new option `runjam_fname_phbin0`
- add options to save single-event phasespace files
  - new option `runjam_output_phbin_indexed`
  - new option `runjam_output_phbin0_indexed`
  - rename option `runjam_ievent_begin` -> `runjam_output_index_start`
- ParticleSampleHydrojet: rename debug options
  - rename option `ParticleSample_ReverseParticleList` -> `hydrojet_reverse_particles`
  - rename option `ParticleSample_ShuffleParticleList` -> `hydrojet_shuffle_particles`
  - rename option `HydroSpectrum_RotateFreezeoutData` -> `hydrojet_rotate_freezeout`

## Changes 0.1..0.2

- General options:
  - Option `-dirJAM PATH`   -> `-o PATH`
  - Option `-f PATH`        -> `--fphase=PATH`
  - Option `-f0 PATH`       -> `--fphase0=PATH`
  - Option `-resodata PATH` -> `--resodata=PATH`
- Options related to `hydrojet` inputs
  - Option `-dir PATH`  -> `--hydrojet-dir=PATH`
  - Option `-ftemp INT` -> `--hydrojet-ftemp=INT`
  - Option `-pce INT`   -> `--hydrojet-pce=INT`
  - Option `-bfree INT` -> `--hydrojet-bfree=INT`
  - Option `-dt NUM`    -> `--hydrojet-dt=NUM`
  - Option `-dx NUM`    -> `--hydrojet-dx=NUM`
  - Option `-dy NUM`    -> `--hydrojet-dy=NUM`
  - Option `-dh NUM`    -> `--hydrojet-dh=NUM`
