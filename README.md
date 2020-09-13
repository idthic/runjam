# runjam

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

`runjam.exe` can be created with the following commands.

```bash
$ git clone https://github.com/idthic/runjam.git
$ cd runjam
$ ./configure --prefix=$HOME/opt/idt --with-jam=$HOME/opt/jam/1.820
$ make
$ make install
```

This programs is provided under [GPLv2](LICENSE).

```
runjam (idt) - Sample hadrons by Cooper-Frye formula / Run JAM cascade.
Copyright (C) 2013-2020, Koichi Murase @akinomyoga

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

Please check the output of `runjam.exe --help`.

### Resonance list

The resonance list can be specified by the option `-r, --resodata, runjam_resodata=FILE`.
This option selects the file that contains the list of particle species to be sampled.
Empty lines and lines starting with `#` are ignored.
Each line contains the information of a particle species with the following format:

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
