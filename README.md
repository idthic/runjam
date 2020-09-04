# hydro2jam

### Option `--resodata, hydro2jam_resodata=FILE`

The list of particles

- Column 1: Mass
- Column 2: Degeneracy
- Column 3: Effective degeneracy (average number of pions after decays)
- Column 4: Chemical potential
- Column 5: Statistics (1: boson, 2: fermion)
- Column 6: Is antiparticle

## Change 0.1..0.2

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
