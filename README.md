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

The list of particles

- Column 1: Mass
- Column 2: Degeneracy
- Column 3: Effective degeneracy (average number of pions after decays)
- Column 4: Chemical potential
- Column 5: Statistics (1: boson, 2: fermion)
- Column 6: Is antiparticle

### Option `--hydrojet-ftemp, hydrojet_kintmp=INT`

- 1 = 80 MeV
- 2 = 100 MeV
- 3 = 120 MeV
- 4 = 140 MeV
- 5 = 160 MeV [default]

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
  - rename option `ParticleSample_ReverseParticleList` -> `hydrojet_shuffle_particles`
  - rename option `HydroSpectrum__RotateFreezeoutData` -> `hydrojet_rotate_freezeout`

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
