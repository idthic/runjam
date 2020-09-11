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

- change variable prefix: `hydro2jam_*` -> `runjam_*`
- rename environment: `ParticleSample_ReverseParticleList` -> `hydrojet_reverse_particles`
- rename environment: `ParticleSample_ReverseParticleList` -> `hydrojet_shuffle_particles`
- fix bugs that some options did not work because of option name typos:
  - `runjam_switch_weak_decay`
  - `--hydrojet-dt, hydrojet_deltat`
  - `--hydrojet-dx, hydrojet_deltax`
  - `--hydrojet-dy, hydrojet_deltay`
  - `--hydrojet-dh, hydrojet_deltah`

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
