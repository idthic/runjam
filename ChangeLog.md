# ChangeLog

## Interface changes 0.2..0.3

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
- resonance list options
  - new option `-r PATH`
  - rename option `--hydrojet-pce, hydrojet_eospce` -> `-p, runjam_eospce`
  - rename option `--hydrojet-ftemp, hydrojet_kintmp` -> `-k, runjam_kintmp`

## Interface changes 0.1..0.2

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
