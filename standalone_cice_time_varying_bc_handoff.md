# Standalone CICE time-varying regional boundary conditions: handoff and development plan

Last updated: 31 July 2026

## Purpose

This document is the handoff for moving the tested PSKRIPS time-varying CICE boundary-condition implementation into a fresh standalone CICE sandbox, then developing it toward the more complete regional-boundary design introduced by [CICE PR #1110](https://github.com/CICE-Consortium/CICE/pull/1110).

The work has two deliberately separate phases:

1. Reproduce the current coupled-model implementation in standalone CICE without changing its behavior.
2. Refactor and extend that standalone implementation, adopting useful upstream architecture while retaining time interpolation, per-edge selection, efficient I/O, and the existing thermodynamic handling.

Do not combine these phases in the first commit. A behavior-preserving standalone baseline is needed to distinguish porting errors from later design changes.

## Scope and terminology

The current PSKRIPS feature reads time-varying sea-ice target data and relaxes the CICE state toward those targets on selected regional boundary edges. It does not change MITgcm ocean OBCS or PWRF coupling.

This document uses three distinct terms:

- **Target acquisition**: read, validate, cache, and time-interpolate the external NetCDF target fields.
- **Exact halo boundary condition**: overwrite exterior ghost cells with prescribed targets at every point in the timestep where CICE halo exchanges would otherwise replace them.
- **Interior restoring/sponge**: relax active-domain cells toward targets over a configurable width and timescale.

The current implementation combines target acquisition with a legacy edge-restoring operation called once near the beginning of each outer CICE timestep. PR #1110 separates these concerns. The standalone development should eventually do the same, but the first port must preserve the current behavior.

## Recommended clean baseline

### CICE and Icepack revisions

Use the official `CICE6.6.3` release for the first standalone port:

| Component | Revision |
| --- | --- |
| Official CICE release | `CICE6.6.3` |
| CICE tag commit | `c734c12bfb19632cd700cb06b22c3ea62b4adeae` |
| Icepack submodule recorded by that release | `daa41638c6cef298583f35ee4dec5ab1bb077aea` |

The exact coupled checkout is commit `52dc4b111c61a3698a814e73e28735a82c52cb4a`, which describes as `CICE6.6.2-14-g52dc4b1`. It is not an official release and its former `regional_coupled` remote branch is no longer present in the CICE-Consortium repository. It consists of the `CICE6.6.2` release, upstream post-release fixes, and local PSKRIPS configuration/driver commits.

`CICE6.6.3` is the preferred clean base because:

- It is the newest official CICE release as of this handoff.
- Its unmodified `ice_restoring.F90` and standalone `CICE_FinalMod.F90` are byte-for-byte identical to those files at the coupled-model commit before the time-varying BC modifications.
- It contains the upstream post-6.6.2 fixes on which the coupled branch was based, plus subsequent released fixes.
- Its newer Icepack still provides the mushy enthalpy interfaces used by `bc_fields_v1`.
- It does not yet contain PR #1110, so the current behavior-preserving port and later upstream-inspired redesign remain separable.

Do not use `CICE6.6.2` merely to match the nearest tag: doing so would require separately identifying and replaying the upstream fixes already incorporated in 6.6.3. Do not start from current `main`, because that would combine this port with newer restoring architecture and unrelated development.

Create the clean sandbox as follows, substituting a new destination for `CICE_WORK`:

```bash
PSKRIPS_ROOT=/nfs/home/paulinan/models/coupled_model_apptainer/pskrips
CICE_WORK=/path/to/fresh/CICE-time-varying-bc
CICE_RELEASE=CICE6.6.3

git clone --branch "$CICE_RELEASE" --recurse-submodules \
  https://github.com/CICE-Consortium/CICE.git "$CICE_WORK"
git -C "$CICE_WORK" submodule update --init --recursive
git -C "$CICE_WORK" switch -c feature/time-varying-regional-bc

git -C "$CICE_WORK" rev-parse HEAD
git -C "$CICE_WORK" submodule status icepack
git -C "$CICE_WORK" status --short
```

The last three commands must report CICE `c734c12bfb19632cd700cb06b22c3ea62b4adeae`, Icepack `daa41638c6cef298583f35ee4dec5ab1bb077aea`, and no local changes before the port begins.

### Compatibility audit against the coupled checkout

The release choice was checked against the coupled CICE source before writing these instructions:

| File/API | CICE 6.6.3 relationship to the coupled baseline | Port consequence |
| --- | --- | --- |
| `ice_restoring.F90` | Byte-for-byte identical before the time-varying changes | The canonical current restorer can be used as a focused replacement |
| Standalone `CICE_FinalMod.F90` | Byte-for-byte identical | Add only `restore_forcing_finalize` import and call |
| Standalone `CICE_InitMod.F90` and `CICE_RunMod.F90` | Release versions retain the existing restore initialization/timestep calls | Preserve the release files; do not copy coupled drivers |
| `ice_forcing.F90` | Coupled baseline has a two-line PSKRIPS-specific JRA grid-name addition; current canonical file also has unrelated ERA5 work | Add only the new BC declarations to the release file |
| `ice_init.F90` | CICE 6.6.3 contains later official initialization changes | Apply BC integration surgically; never replace the release file wholesale |
| Icepack mushy enthalpy API | `icepack_enthalpy_snow` and `icepack_enthalpy_mush` are present in the release's Icepack | `bc_fields_v1` is source-compatible, subject to build and runtime tests |

This audit supports using 6.6.3, but it does not replace the standalone compile, disabled regression, schema tests, or restart-exactness gates below.

### Source of truth for the present implementation

The canonical coupled source is:

```text
models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice
```

Do not use the build-staged copies under the CICE submodule. Some canonical files contain working-tree changes not represented by the coupled CICE commit, so use the current files themselves rather than `git show HEAD:<path>`.

At the time of this handoff, the relevant canonical checksums are:

| File | SHA-256 |
| --- | --- |
| `ice_restore_forcing.F90` | `e42fac6766b655369c78d1b9736d9089bc3dfc61b3a60293ed5d492a2110816b` |
| `ice_restoring.F90` | `bc7a02f6fb41b5068c8d12fb30e5aedce4d907ddf6d0f806c10387729d2a1bc7` |
| `ice_init.F90` | `0a461259c72885112200eeb6e2d054ce7a9b28673bd6911b7d9b809e4ce1edc8` |
| `ice_forcing.F90` | `0994a4e432b462de35d2a51edf1b778133960e65d650de6100d9f628c9431189` |
| `CICE_FinalMod.F90` | `762270ba924855f312b4ff6443fc20eda6c8820f218d40b3d5ef70af4a9ab980` |

Verify them before starting:

```bash
sha256sum \
  "$PSKRIPS_ROOT"/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice/ice_restore_forcing.F90 \
  "$PSKRIPS_ROOT"/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice/ice_restoring.F90 \
  "$PSKRIPS_ROOT"/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice/ice_init.F90 \
  "$PSKRIPS_ROOT"/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice/ice_forcing.F90 \
  "$PSKRIPS_ROOT"/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice/CICE_FinalMod.F90
```

If a checksum has changed, review the newer diff and update this handoff rather than assuming the old file map remains valid.

## What the current implementation does

The behavior to preserve in the initial standalone port is:

- Daily or monthly linear interpolation using record indices.
- Cyclic climatologies or finite non-cyclic time series.
- Independent west, east, south, and north enable flags.
- Master-rank NetCDF reads, scatter to local blocks, and two-record caching.
- Cached block geometry and Icepack metadata to avoid repeated queries and block discovery.
- Three forcing representations:
  - `full_trcrn`: `aicen`, `vicen`, `vsnon`, and `trcrn`.
  - `bc_fields_v1`: `aicen`, `vicen`, `vsnon`, `Tsfc`, `Tinz`, and `Sinz`; CICE reconstructs thermodynamic tracers using the configured mushy-layer physics.
  - `structure_only`: `aicen`, `vicen`, and `vsnon`; no thermodynamic target is read from the file. The present implementation snapshots tracer values on the physical edge and extends them outward, so physical-edge tracer values are not directly restored while ghost tracers follow that zero-gradient target.
- Exact schema detection: partial `Tsfc/Tinz/Sinz` groups and ambiguous files containing both `trcrn` and the complete thermodynamic group are rejected.
- Exactly 365 records for cyclic daily forcing and exactly 12 for cyclic monthly forcing.
- Runtime coverage checks for non-cyclic forcing.
- The index-time convention: daily record 0 is 1 January of `fyear_init`; monthly record 0 is January of `fyear_init`.
- Land masking and extension from physical edge cells into the exterior ghost cells used by the legacy restorer.
- Finalization that closes the NetCDF handle and deallocates caches before CICE releases file units.

Important limitations of the baseline are:

- It applies the legacy restoring operation once near the beginning of the outer timestep. Later `bound_state` calls can replace halo values.
- Its halo-restoring loops explicitly note that they need modification for `nghost > 1`.
- The restoring work is charged to `timer_bound`, making boundary exchange and restoring cost hard to distinguish.
- It does not prescribe ice velocity.
- It does not provide a configurable interior sponge mask or width.
- The current category-aware thermodynamic update is more careful than upstream's initial PR implementation, but edge cases involving categories appearing or disappearing still need targeted tests.

## Part I: port the current feature into standalone CICE

### Port map

| Coupled canonical source | Standalone destination | Action |
| --- | --- | --- |
| `ice_restore_forcing.F90` | `cicecore/cicedyn/infrastructure/ice_restore_forcing.F90` | Add exact current file |
| `ice_restoring.F90` | `cicecore/cicedyn/infrastructure/ice_restoring.F90` | Replace the 6.6.3 file with the exact current implementation |
| `ice_forcing.F90` | `cicecore/cicedyn/general/ice_forcing.F90` | Apply only the BC declarations; omit ERA5 changes |
| `ice_init.F90` | `cicecore/cicedyn/general/ice_init.F90` | Apply the BC namelist/default/broadcast/diagnostic hunks |
| `CICE_InitMod.F90` | `cicecore/drivers/standalone/cice/CICE_InitMod.F90` | No change; verify its existing initialization call |
| `CICE_RunMod.F90` | `cicecore/drivers/standalone/cice/CICE_RunMod.F90` | No change; verify its existing timestep call |
| `CICE_FinalMod.F90` | `cicecore/drivers/standalone/cice/CICE_FinalMod.F90` | Add only reader finalization |
| Machine `ice_in` templates | `configuration/scripts/ice_in` | Add disabled defaults, then enable them only in test cases |
| PSKRIPS converters/validator/tests | Standalone tools/tests location | Copy or pin as a versioned external dependency |

### 1. Make a baseline commit before changing behavior

Record the CICE/Icepack revisions and create a clean branch. After each logical step below, commit separately. The recommended initial commit sequence is:

1. `Add time-varying restore forcing reader`
2. `Expose restore forcing namelist options`
3. `Integrate restore reader with standalone lifecycle`
4. `Add forcing conversion and validation tools`
5. `Add standalone regional BC tests`

This makes later upstream-inspired refactoring reviewable and reversible.

### 2. Copy the two BC implementation files

In a CICE 6.6.3 checkout, the complete `ice_restoring.F90` replacement is boundary-specific and can be copied as a unit. The unmodified 6.6.3 restorer is byte-identical to the unmodified restorer at the coupled-model commit, so this diff isolates the current time-varying implementation. Add the new reader beside it:

```bash
BC_SOURCE="$PSKRIPS_ROOT/models/PSKRIPS/PSKRIPSv3_CICE/ciceCode/mitgcm_cice"

cp "$BC_SOURCE/ice_restore_forcing.F90" \
  "$CICE_WORK/cicecore/cicedyn/infrastructure/ice_restore_forcing.F90"
cp "$BC_SOURCE/ice_restoring.F90" \
  "$CICE_WORK/cicecore/cicedyn/infrastructure/ice_restoring.F90"
```

The standard CICE `Filepath` includes `cicecore/cicedyn/infrastructure`, and the dependency generator discovers all `.F90` files in that path. No hard-coded source list should be necessary. Confirm this during the build rather than editing the Makefile pre-emptively.

Review the replacement before committing:

```bash
git -C "$CICE_WORK" diff --stat
git -C "$CICE_WORK" diff -- \
  cicecore/cicedyn/infrastructure/ice_restoring.F90 \
  cicecore/cicedyn/infrastructure/ice_restore_forcing.F90
```

### 3. Add only the BC declarations to `ice_forcing.F90`

Do **not** copy the coupled `ice_forcing.F90` wholesale. It also contains approximately 160 lines of unrelated, machine-specific ERA5 binary-forcing support. The standalone BC port only needs these declarations in module `ice_forcing`:

```fortran
      character(char_len_long), public :: &
         restore_ice_data_file  ! file name for time-varying ice boundary data

      character(char_len), public :: &
         restore_ice_data_type   ! 'legacy', 'daily_netcdf', or 'monthly_netcdf'

      logical (kind=log_kind), public :: &
         restore_ice_cycle_year, &
         restore_ice_use_west  , &
         restore_ice_use_east  , &
         restore_ice_use_south , &
         restore_ice_use_north
```

Place the character variables next to the existing forcing file names and the logical variables next to `restore_ocn`; retain the existing `trestore` declaration. Verify that the resulting diff contains no `ERA5_files`, `ERA5_data`, or `file_year` changes.

If standalone regional experiments later require the custom ERA5 reader, port it in a separate commit and test it independently of the boundary work.

### 4. Expose and distribute the new namelist values in `ice_init.F90`

Do not copy the canonical coupled `ice_init.F90` over CICE 6.6.3. The release contains later official initialization, history-restart, wave, grid-format, and residual-ice changes that must be preserved. Use the canonical file only to locate the BC additions:

```bash
git diff --no-index -- \
  "$CICE_WORK/cicecore/cicedyn/general/ice_init.F90" \
  "$BC_SOURCE/ice_init.F90"
```

That comparison will contain both the desired BC integration and unrelated release differences. Apply only the six BC integration points below to the 6.6.3 file.

The completed change must include all six integration points:

1. Import from `ice_forcing`:
   - `restore_ice_data_type`
   - `restore_ice_data_file`
   - `restore_ice_cycle_year`
   - `restore_ice_use_west/east/south/north`
2. Add all seven values to `forcing_nml`.
3. Set safe defaults:

   ```fortran
   restore_ice_data_type = 'legacy'
   restore_ice_data_file = 'unknown_restore_ice_file'
   restore_ice_cycle_year = .true.
   restore_ice_use_west  = .true.
   restore_ice_use_east  = .true.
   restore_ice_use_south = .true.
   restore_ice_use_north = .true.
   ```

4. Broadcast each value from the master rank.
5. Print the selected values to the diagnostic log when `restore_ice=.true.`.
6. Preserve the existing disabled default `restore_ice=.false.` and `trestore=90`.

Also add the new options to `configuration/scripts/ice_in` beside the existing `restore_ice` setting so newly generated standalone cases expose the interface:

```fortran
    restore_ice             = .false.
    restore_ice_data_type   = 'legacy'
    restore_ice_data_file   = 'unknown_restore_ice_file'
    restore_ice_cycle_year  = .true.
    restore_ice_use_west    = .true.
    restore_ice_use_east    = .true.
    restore_ice_use_south   = .true.
    restore_ice_use_north   = .true.
```

Keeping all edge flags true by default preserves the old all-edge behavior when the feature is enabled, while the master `restore_ice=.false.` default preserves normal CICE behavior.

### 5. Integrate with the standalone driver lifecycle

Do not copy the coupled `CICE_InitMod.F90` or `CICE_RunMod.F90`. They contain coupled-driver behavior that does not belong in standalone CICE. CICE 6.6.3 already has the required calls:

```bash
rg -n "ice_HaloRestore_init|ice_HaloRestore" \
  "$CICE_WORK/cicecore/drivers/standalone/cice/CICE_InitMod.F90" \
  "$CICE_WORK/cicecore/drivers/standalone/cice/CICE_RunMod.F90"
```

Expected behavior:

- `CICE_InitMod.F90` calls `ice_HaloRestore_init` once.
- `CICE_RunMod.F90` calls `ice_HaloRestore` once per outer timestep when `restore_ice` is true.

Make only the following standalone finalization change in `cicecore/drivers/standalone/cice/CICE_FinalMod.F90`:

```fortran
      use ice_restore_forcing, only: restore_forcing_finalize
```

and, before `release_all_fileunits`:

```fortran
      call restore_forcing_finalize()
      call release_all_fileunits
```

Do not copy the canonical coupled finalizer wholesale; it contains coupled preprocessor guards and unrelated formatting.

### 6. Port the external forcing tools

The present tools live outside the CICE submodule:

```text
environment.cice-setup.yml
tools/cice_setup_python.sh
tools/cice_thermo_utils.py
tools/convert_cice_history_to_restore_forcing.py
tools/convert_climate_ice_to_restore_forcing.py
tools/validate_cice_restore_forcing.py
tests/test_cice_restore_tools.py
```

Copy these into an appropriate `tools/restore_forcing` and test location in the standalone development repository, or keep them in a sibling utilities repository and record that revision in the standalone case. Do not duplicate them silently: there should be one clearly versioned converter/validator implementation associated with the Fortran schema.

Create the Python environment before running the tests:

```bash
cd "$PSKRIPS_ROOT"
mamba env create -f environment.cice-setup.yml
./tools/cice_setup_python.sh -m unittest discover \
  -s tests -p 'test_cice_restore_tools.py' -v
```

The environment provides Python 3.11, NumPy, SciPy, pandas, xarray, dask, netCDF4, cftime, cf_xarray, xESMF, and ESMPy. The test command currently cannot run on a host where the `cice_setup` mamba environment has not been created; that is an environment prerequisite, not a test skip.

The tools must retain these capabilities:

- Daily history input: `full_trcrn`, `bc_fields_v1`, or `structure_only`.
- Daily restart input: `full_trcrn` or `structure_only`.
- Monthly climate concentration, ice volume, and snow input: regrid with xESMF and write `structure_only` using the established five-category thickness template.
- Resolve the target grid from the active `ice_in`, with an explicit target-grid override available.
- Require or derive all paths; do not add machine-specific defaults.
- Write zero-based record-index time coordinates anchored to January of `fyear_init`.
- Validate dimensions, schema completeness, finite and nonnegative values, summed concentration, record counts, time anchoring, and non-cyclic simulation coverage.

See [cice_time_varying_boundaries.md](cice_time_varying_boundaries.md) for the current converter and validator command lines.

### 7. Configure a standalone regional case

The first regional standalone case should match the coupled CICE configuration as closely as possible:

```fortran
&grid_nml
    grid_format = 'nc'
    grid_type   = 'regional'
    grid_ice    = 'C'
    grid_atm    = 'A'
    grid_ocn    = 'C'
    grid_file   = '/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/grid_pskrips.nc'
    kmt_type    = 'file'
    kmt_file    = '/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/kmt_pskrips.nc'
    ncat        = 5
    nilyr       = 7
    nslyr       = 1
/

&domain_nml
    nprocs            = 12
    nx_global         = 210
    ny_global         = 240
    block_size_x      = 35
    block_size_y      = 30
    processor_shape   = 'slenderX2'
    distribution_type = 'cartesian'
    ew_boundary_type  = 'open'
    ns_boundary_type  = 'open'
    maskhalo_dyn      = .true.
    maskhalo_remap    = .true.
    maskhalo_bound    = .true.
/
```

Use the existing regional grid, mask, and initial condition as external test data. Do not commit large forcing, regridding-weight, grid, or restart files to the source repository.

For the known `bc_fields_v1` dataset, use:

```fortran
&forcing_nml
    restore_ice             = .true.
    restore_ice_data_type   = 'daily_netcdf'
    restore_ice_data_file   = '/nfs/scratch/paulinan/cice-dirs/input/CICE_data/bdry/pskrips/cice_ERA5_ECCO_2017_bc_fields_v1_2017.nc'
    restore_ice_cycle_year  = .true.
    restore_ice_use_west    = .true.
    restore_ice_use_east    = .true.
    restore_ice_use_south   = .true.
    restore_ice_use_north   = .true.
    trestore                = 90
/
```

The known coupled test used:

- Start: `2017-02-01 00:00:00`
- `dt=60 s`
- `restart_ext=.false.`
- `ktherm=2`
- `advection='remap'`
- five ice categories and seven ice layers
- the regional initial condition `iced_pskrips_2017-02-01.onice_hi0p5_hs0p01.thin_tail_itd_minarea0p05.nc`

A NetCDF time-varying forcing file is a valid source with `restart_ext=.false.`. The legacy extended-restart source still requires `restart_ext=.true.` for open boundaries.

Do not copy the coupled case namelist blindly. In standalone mode, PWRF and MITgcm no longer supply atmospheric and ocean fields. For mechanical boundary tests, CICE's deterministic `atm_data_type='default'` and `ocn_data_type='default'` are suitable if their default fields meet the test objective. For science comparisons, construct or select standalone atmospheric/ocean forcing explicitly and hold it identical between control and boundary-enabled runs.

### 8. Build checks

Use a NetCDF-enabled CICE build. On Raapoi, add or adapt machine files under `configuration/scripts/machines/`, create a standalone case with `cice.setup`, and build from the generated case directory. A representative sequence is:

```bash
cd "$CICE_WORK"
./cice.setup --case /path/to/cases/pskrips_bc_smoke \
  --mach raapoi --env native --pes 12x1 --grid gx1
cd /path/to/cases/pskrips_bc_smoke
./cice.build
```

The current local `Macros.raapoi_native` and `env.raapoi_native` are useful references, but they are untracked build-port files in the coupled CICE checkout. Put reviewed copies in `configuration/scripts/machines/` in the standalone branch, ensure `PSKRIPS_NATIVE_PREFIX` points to the native NetCDF installation, and avoid copying the coupled `cice.settings` because it selects the NUOPC `mitgcm_cice` driver and `libcice` target.

After building, verify:

```bash
find /path/to/cice/run -name 'ice_restore_forcing.o' -size +0 -print
rg -n "ice_restore_forcing|ice_restoring" /path/to/cases/pskrips_bc_smoke/logs/cice.bldlog.*
```

Also confirm that the final linked target is the standalone `cice` executable, not only `libcice.a`, and that the build uses the standalone driver path.

Perform one debug build with bounds, uninitialized-variable, and floating-point exception checks before performance testing. Perform production comparisons with the same compiler and optimization flags used for the control.

### 9. Initial standalone validation sequence

Run these tests in order. Do not begin the upstream-inspired refactor until all acceptance gates pass.

#### Gate A: reader and converter tests

- Python converter/validator unit tests pass.
- The known forcing file passes the read-only validator against the active standalone `ice_in`.
- Synthetic files for all three schemas are accepted.
- Partial thermodynamic groups, ambiguous schemas, wrong dimensions, invalid concentrations, non-finite values, and bad time coverage are rejected.

Example validation:

```bash
cd "$PSKRIPS_ROOT"
./tools/cice_setup_python.sh tools/validate_cice_restore_forcing.py \
  --forcing-file /nfs/scratch/paulinan/cice-dirs/input/CICE_data/bdry/pskrips/cice_ERA5_ECCO_2017_bc_fields_v1_2017.nc \
  --ice-in /path/to/cases/pskrips_bc_smoke/ice_in \
  --data-type daily_netcdf --cycle-year
```

#### Gate B: disabled regression

- Build the modified executable.
- Run `restore_ice=.false.` and `restore_ice_data_type='legacy'`.
- Compare against an unmodified standalone executable using identical namelists and forcing.
- Require bit-for-bit restart equivalence and exact history values after ignoring creation metadata.

This proves that adding the module and finalizer does not perturb normal CICE runs.

#### Gate C: legacy time-varying behavior

- Run multi-rank `bc_fields_v1` forcing across at least one daily record transition.
- Confirm the reader loads the expected adjacent record pair and linearly interpolates the midpoint.
- Confirm only selected edges change when each edge flag is tested independently.
- Confirm land cells remain zero and padded/ghost cells are populated from the correct physical edge.
- Repeat with `full_trcrn` and monthly `structure_only` forcing.
- For `structure_only`, assert thermodynamic tracer values in the restored physical-edge cells are unchanged by the restoring operation, and separately verify the documented zero-gradient tracer extension into ghost cells.

#### Gate D: restart exactness

- Run a continuous case across a forcing-record transition.
- Split an otherwise identical run before that transition and continue from restart.
- Require byte-identical final CICE restart fields and exact history values after ignoring metadata.

The coupled implementation already passed this pattern: the final CICE restarts were byte-identical, CICE history values matched exactly, and MITgcm ice-state outputs and mediator restarts matched. Standalone CICE must reproduce the CICE portion independently.

#### Gate E: decomposition and performance

- Repeat at no fewer than two MPI decompositions and require identical global results.
- Record `timer_bound`, total timestep time, record-load count, and resident memory.
- Run long enough to include cached timesteps and at least one new-record read.
- Compare enabled and disabled runs using the same executable.

The coupled optimization reduced the incremental boundary cost from approximately 58.10 s to 16.73 s in the measured run, about 71.2%, while retaining exact outputs. Treat that as evidence that metadata caching and boundary-only interpolation matter, not as a standalone performance threshold. Establish a fresh standalone baseline.

### 10. Commit the reproducible standalone baseline

The baseline milestone is complete only when the branch contains:

- The CICE 6.6.3 and corresponding Icepack revisions in its provenance record.
- The new forcing reader and the current legacy restorer.
- Minimal `ice_forcing`, `ice_init`, and standalone finalizer integration.
- Disabled defaults in the standard `ice_in` template.
- Versioned converter and validator tools, or an explicit pinned external-tools dependency.
- A small synthetic forcing fixture or fixture generator suitable for CI.
- A documented regional case manifest listing external files and their checksums.
- Passing Gates A-E.

Tag or otherwise identify this point before changing the boundary semantics. A suggested tag is `time-varying-bc-standalone-legacy-v1`.

## Part II: what to adopt from CICE PR #1110

PR #1110 should be treated as an architectural reference, not a patch to cherry-pick into this older CICE revision. It targets a substantially newer codebase and changes many files. Its data-source implementation also does not replace the present time-varying reader. The assessment used PR head `1032b33f6f427b8dcd882afddc8a2ce3825999a9`; GitHub records merge commit `52cb686ad3e0801f24ddabdc8ba2413efe24fa04`.

### Adopt these ideas

1. **Separate data acquisition, exact halo assignment, and interior restoring.**

   The upstream design exposes the conceptual equivalents of `ice_restoring_getdata`, `ice_restoring_halo`, and `ice_restoring_interior`. This separation makes call timing, field selection, and testing much clearer.

2. **Reapply prescribed halo values after CICE overwrites halos.**

   The most important lesson is that a regional boundary condition cannot be guaranteed by setting the ghost cells only once at the start of the outer timestep. `bound_state` and velocity halo exchanges occur again during thermodynamics, dynamics, and transport. Exact prescribed halos must be reinstated immediately after the relevant exchanges.

3. **Separate exact exterior boundaries from an optional interior sponge.**

   Upstream supports boundary-field selection plus independently configured interior restoring fields, masks, widths, and timescales. This is a cleaner physical model than using one relaxation operation for both purposes.

4. **Use field selection.**

   Upstream can independently select state, category area, ice volume, snow volume, tracers, and velocity. We should support a compatible field-selection concept so experiments can distinguish structural, thermodynamic, and dynamic boundary effects.

5. **Use a dedicated restoring timer.**

   Restoring should not be hidden in `timer_bound`. At minimum, add a `timer_restore`; preferably report target update/I/O, thermodynamic reconstruction, halo application, and interior-restoring time separately during development.

6. **Make halo loops generic in `nghost`.**

   Upstream avoids the current one-ghost-cell assumption. All new exact-halo logic should derive ranges from CICE block geometry and `nghost`, then test more than one ghost width.

7. **Add full-to-nest tests and test the advection choice.**

   Upstream testing found much larger errors with remapping advection than with upwind in one nested-boundary experiment (roughly 0.12 versus 0.003 in the reported comparison). The discussion suspected technical/advection halo treatment rather than only the restoring formula. PSKRIPS currently uses `advection='remap'`, so this is a priority science test.

### Retain or improve these PSKRIPS capabilities

Do not lose the following when adopting the upstream structure:

- Daily and monthly time interpolation.
- Cyclic and non-cyclic series.
- A single multi-record forcing file rather than reopening a separately named restart file every timestep.
- Two-record caching and master-rank reads.
- Per-edge enable flags.
- Three forcing schemas, including `bc_fields_v1` and intentional `structure_only` support.
- Strict file/schema validation.
- Category-thickness preservation in the direct thermodynamic-restoring path.
- Clean finalization and disabled-mode regression behavior.
- The optimization that avoids repeated metadata queries and full-domain interpolation.

### Do not adopt these PR details directly

- Do not require `restart_ext=.true.` for NetCDF time-varying forcing. That upstream condition follows from its extended-restart data source, not from open boundaries in general.
- Do not adopt the hard-coded, per-timestep extended-restart filename scheme as the primary data source.
- Do not remove per-edge selection.
- Do not replace the present category-aware thermodynamic logic with unconditional linear blending.
- Do not make exact halo values depend on the legacy `trestore` timescale. Exact halo assignment and interior relaxation have different semantics.
- Do not add velocity restoration until the forcing schema, grid staggering, rotation, and validation rules are explicitly defined.

## Detailed staged implementation plan

Each stage below has its own acceptance gate. Preserve the legacy standalone tag as the comparison baseline throughout.

### Stage 0: freeze behavior and terminology

Before refactoring:

- Save small expected-output fixtures for each schema.
- Save a continuous/restart-split reference pair.
- Save per-edge and interpolation diagnostics.
- Record the standalone performance profile.
- Document whether target values represent cell means at record times and confirm that the existing linear interpolation convention is scientifically intended.

Rename nothing yet. The purpose of Stage 0 is a stable oracle.

**Acceptance:** all Part I gates pass, and the result is tagged.

### Stage 1: add observability without changing results

Add a dedicated `timer_restore` to the CICE timer infrastructure. During development, split or count:

- target-record selection;
- NetCDF record loads;
- scatter/cache updates;
- target interpolation;
- `bc_fields_v1` thermodynamic reconstruction;
- boundary application;
- interior application, once it exists.

Add counters printed at finalization or diagnostic intervals:

- number of model updates;
- number of record loads;
- cached-pair hits;
- number of boundary cells updated by edge and field group.

Counters should be cheap and disabled or low-overhead in production.

**Acceptance:** restart and history outputs remain bit-for-bit identical to the Stage 0 baseline. The record-load count changes only at forcing transitions.

### Stage 2: refactor into three APIs without changing legacy behavior

Refactor `ice_restoring` around three responsibilities:

```fortran
call ice_restore_targets_update(model_time)
call ice_restore_halo(field_group)
call ice_restore_interior(field_group, dt)
```

The names may follow current CICE style, but the responsibilities should be explicit:

- `targets_update` selects/interpolates records once per outer timestep and exposes immutable target arrays for that timestep.
- `halo` assigns exterior ghost cells only.
- `interior` modifies physical-domain cells only.

Initially, retain a compatibility wrapper:

```fortran
call ice_HaloRestore
```

that invokes the refactored code in exactly the old order and reproduces legacy results. Move the reader and cache implementation behind a small query API rather than letting dynamics code know about NetCDF.

Keep target arrays distinct from model state. Mark intent carefully so exact-halo routines cannot accidentally mutate the cached targets.

**Acceptance:** all Stage 0 results remain bit-for-bit identical in `legacy_nudged` mode, with no performance regression larger than measurement noise.

### Stage 3: introduce an explicit boundary mode

Add a mode switch with backward-compatible default:

```fortran
restore_ice_mode = 'legacy_nudged'
```

Proposed values:

- `legacy_nudged`: current behavior; `trestore` applies to the legacy edge operation.
- `prescribed_halo`: exact exterior halo values; no interior relaxation unless separately enabled.
- `prescribed_halo_sponge`: exact exterior halo values plus independently configured interior restoring.

Reject unknown modes and incompatible combinations at startup. Print a clear summary of the resolved behavior.

For future alignment with upstream naming, introduce field selectors modeled on `set_boundary_flds` and `restore_flds`, while retaining explicit PSKRIPS schema semantics. One workable interface is:

```fortran
set_boundary_flds = 'state'
restore_flds      = 'none'
```

Support explicit tokens rather than ambiguous substrings. Define `state` in documentation as `aicen`, `vicen`, `vsnon`, plus thermodynamic tracers only when the selected forcing schema supplies them. A `structure_only` file must never silently become full thermodynamic state.

**Acceptance:** default namelists still select legacy behavior; selecting prescribed mode changes only the intended halo cells in a one-step test.

### Stage 4: implement exact prescribed state halos

Implement state-halo assignment for `aicen`, `vicen`, `vsnon`, and schema-authorized tracers. Use cached block metadata, but replace all one-cell assumptions with ranges derived from `nghost`, `ilo/ihi`, `jlo/jhi`, and actual east/north padding.

Apply exact state halos immediately after every `bound_state` that can affect the subsequent computation. In CICE 6.6.3, audit at least:

- initialization in `ice_init.F90`;
- `update_state` in `ice_step_mod.F90`;
- upwind and remap transport paths in `ice_transport_driver.F90`;
- any ridging, wave, or restart path that exchanges state before consuming boundary values.

Avoid scattering NetCDF data on every subcycle. Update/interpolate targets once per outer timestep, then reapply already-local target arrays at each required halo point.

Corner policy must be explicit. Recommended behavior:

- A corner is active if either adjoining enabled edge is active.
- If both edges prescribe the same cell from the same gridded target, assign once.
- Never apply relaxation twice at corners.

After exact assignment, recalculate aggregates only where downstream code requires them. Do not perform unnecessary whole-domain aggregation for a boundary strip.

**Acceptance:** instrumented assertions immediately before consumers show exact target values on all enabled exterior halo cells after every relevant exchange. Disabled edges retain native CICE halo behavior. Tests pass for `nghost=1`, `2`, and `3` if supported by the configuration.

### Stage 5: add optional interior sponge restoring

Implement interior restoring separately from halo assignment, following the useful PR #1110 concepts:

```fortran
restore_flds     = 'state'
restore_mask     = 'linear'   ! 'none', 'all', 'constant', or 'linear'
restore_width    = N          ! physical-domain cells inward from enabled edges
restore_timescale = seconds   ! use an unambiguous unit
```

Define the weight mathematically. For example, for a linear mask of width `W`, let the boundary-adjacent physical cell have weight 1 and the last sponge cell have weight `1/W`, then apply a stable relaxation factor derived from `dt/timescale`. Specify the instantaneous limit and reject negative timescales.

Build the sponge mask once from global edge distance, land mask, enabled-edge flags, and decomposition metadata. At corners, combine edge weights with `max`, not a sum, to prevent double restoring. Apply each field group once per outer timestep, not once per dynamics subcycle.

Keep `trestore` only for legacy compatibility. Use a new explicitly documented timescale for the new sponge API, preferably seconds to avoid hidden day conversions.

**Acceptance:** analytic small-grid tests reproduce the expected mask and exponential/linear update; MPI decompositions produce identical global masks and results; `restore_width=0` or `restore_flds='none'` leaves the active domain untouched.

### Stage 6: formalize thermodynamic category semantics

Use the discussion around PR #1110 as a warning that tracer restoration is not just component-wise blending. Implement and test these rules explicitly:

- `full_trcrn`: targets all supplied tracer values, subject to CICE consistency checks.
- `bc_fields_v1`: reconstructs the thermodynamic tracer subset from `Tsfc/Tinz/Sinz` with the active Icepack thermodynamics, while preserving non-thermodynamic tracer behavior unless separately selected.
- `structure_only`: changes structure but preserves tracer densities/properties on pre-existing ice.
- When a category gains ice from effectively zero concentration, initialize its thermodynamics from an explicit source: supplied target thermo, a documented neighboring/source category rule, or a physically valid default. Never divide by a near-zero old concentration.
- When a category loses all ice, clear dependent extensive quantities consistently.
- Enforce bounds and call the relevant Icepack consistency logic after category creation/removal.

Add cell-level tests for:

- unchanged category concentration with changed volume;
- concentration increase/decrease;
- category creation and deletion;
- snow appearing/disappearing;
- near-zero `aicen`;
- enthalpy and salinity bounds;
- conservation behavior expected from each mode.

Preserve the current target category thickness behavior in legacy mode. Any new semantics should be selected explicitly and documented.

**Acceptance:** all three schemas pass category-transition tests without NaNs, invalid Icepack states, or unintended changes to non-selected tracers.

### Stage 7: decide whether to prescribe velocity

PR #1110 supports velocity boundary selection, but the current forcing schemas do not contain CICE `uvel/vvel`. Treat velocity as a separate feature.

Before implementation, specify:

- variable names and NetCDF dimensions;
- whether values are on the CICE U grid or a source grid;
- rotation convention and units;
- how corner and land points are handled;
- whether velocities are instantaneous targets or time means;
- the locations of every velocity halo update in the dynamics subcycling.

Extend the converter and validator only after this schema is fixed. If velocity is not supplied, retain native CICE velocity halo behavior and log that choice.

**Acceptance:** a solid-body or uniform-flow synthetic test verifies staggering, sign, rotation, edge selection, and exact reapplication after velocity halo exchanges.

### Stage 8: full-to-nest regional science tests

Create a parent standalone CICE run and a nested/regional run whose target data are extracted from the parent. Prefer writing the parent data to the established multi-record forcing format rather than adopting per-timestep restart filenames.

Run this matrix:

| Dimension | Required cases |
| --- | --- |
| Boundary semantics | legacy nudged, exact halo, exact halo plus sponge |
| Advection | remap, upwind |
| Thermodynamics | `full_trcrn`, `bc_fields_v1`, `structure_only` where scientifically meaningful |
| Edges | all, each single edge, representative adjacent pair |
| MPI layout | at least two decompositions |
| Restart | continuous and split across a record transition |
| Ghost width | every supported `nghost`, at least 1 and 2 |

Diagnostics should include:

- halo-minus-target norms before and after exact reapplication;
- active-domain error by distance from each boundary;
- category area/volume/snow error;
- thermodynamic tracer error;
- ice velocity error if enabled;
- total-area and total-volume budgets;
- errors before and after transport;
- remap technical/metric halo diagnostics near open boundaries;
- performance by restoring sub-timer.

The remap/upwind discrepancy reported during PR #1110 review makes the transport-stage diagnostics essential. If exact state halos do not reduce remap errors, inspect metric, face, and technical halos used by remap rather than increasing the sponge strength without diagnosis.

**Acceptance:** define quantitative tolerances from the parent/nest experiment before selecting the production boundary mode. Exact halos must be exact at consumption points; active-domain error tolerances may differ by advection method but must be documented and scientifically justified.

### Stage 9: harden the data interface

Before returning the work to the coupled model:

- Add a schema/version global attribute and preserve strict auto-detection as a cross-check.
- Record source CICE/Icepack revisions, source files, grid checksum, converter revision, and generation command in NetCDF attributes or a sidecar manifest.
- Decide how leap years are handled. The current cyclic daily contract is exactly 365 records; a future 366-day mode must be an explicit new convention, not inferred silently.
- Add chunking/compression performance tests for realistic regional files.
- Consider parallel or rank-local I/O only if profiling shows master-read/scatter is a bottleneck; retain the simple cached implementation otherwise.
- Keep validator checks synchronized with every Fortran schema change.
- Ensure all error messages identify the offending file, variable, expected dimensions, and model configuration.

**Acceptance:** malformed and provenance-incomplete files fail early with actionable messages; valid production-scale files meet the standalone I/O target.

### Stage 10: bring the proven design back to PSKRIPS

Only after standalone acceptance:

1. Identify the standalone commits containing reader changes, pure refactors, exact-halo call sites, sponge support, and science fixes.
2. Port those commits into the canonical PSKRIPS CICE source directory, not the staged submodule build copy.
3. Reconcile standalone driver call sites with the one-step-per-call coupled driver without changing target update frequency.
4. Rebuild `libcice` and the complete coupled executable.
5. Repeat the coupled disabled control, all three schema tests, multi-rank tests, performance comparison, and continuous/restart-split comparison.
6. Check PWRF and MITgcm outputs to ensure CICE-only changes do not introduce coupled restart divergence.

Keep the legacy mode available for one full comparison cycle. Change the production default only after the exact-halo/sponge configuration has passed the nested science tests.

## Proposed test inventory

The standalone branch should eventually automate the following:

| Test | Purpose | Required result |
| --- | --- | --- |
| Modified executable, BC off | Non-regression | Bit-for-bit with unmodified CICE |
| Daily interpolation endpoints | Record mapping | Exact record values |
| Daily interpolation midpoint | Linear interpolation | Analytic midpoint |
| Monthly interpolation and year wrap | Monthly/cyclic logic | Analytic values |
| Non-cyclic final coverage | Bounds checking | Valid endpoint accepted; out-of-range time rejected |
| Each forcing schema | Schema semantics | Correct state/tracer selection |
| Partial thermo group | Schema safety | Startup failure |
| Each single edge | Edge selection | Only requested edge modified |
| Adjacent edges and corners | Corner policy | One consistent update, no double relaxation |
| Land/padding | Geometry safety | No ice on land; no padding corruption |
| Multiple MPI layouts | Decomposition independence | Identical global fields |
| Multiple ghost widths | Halo generality | Correct prescribed strips |
| Continuous vs split restart | Restart exactness | Bit-for-bit final state |
| Category creation/removal | Thermodynamic safety | Valid bounded state |
| Full-to-nest upwind/remap | Regional science | Quantified, acceptable errors |
| Enabled/disabled same executable | Performance | Measured incremental cost |
| Record transition profile | Cache efficiency | One new record load, no redundant full-file reads |

## External data and known coupled evidence

Known external test data:

```text
/nfs/scratch/paulinan/cice-dirs/input/CICE_data/bdry/pskrips/cice_ERA5_ECCO_2017_bc_fields_v1_2017.nc
/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/grid_pskrips.nc
/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/kmt_pskrips.nc
/nfs/scratch/paulinan/cice-dirs/input/CICE_data/ic/pskrips/iced_pskrips_2017-02-01.onice_hi0p5_hs0p01.thin_tail_itd_minarea0p05.nc
```

Known coupled run directories, useful for provenance and comparison but not substitutes for standalone tests:

```text
/nfs/scratch/paulinan/pskrips/bctest
/nfs/scratch/paulinan/pskrips/bcofftest
/nfs/scratch/paulinan/pskrips/perftest
/nfs/scratch/paulinan/pskrips/bcoffperftest
/nfs/scratch/paulinan/pskrips/bcresttest
/nfs/scratch/paulinan/pskrips/bcconttest
```

The completed coupled tests established that:

- the reader works with the regional five-category, seven-layer `bc_fields_v1` data on multiple ranks;
- a boundary-disabled control completes normally;
- cached metadata and boundary-only interpolation substantially reduce incremental cost while retaining exact output;
- a restart split across a forcing transition reproduces the continuous CICE state exactly.

These results reduce porting risk, but standalone integration, exact-halo semantics, sponge behavior, multiple ghost widths, and full-to-nest remap behavior remain to be demonstrated.

## Decision points that require explicit scientific choices

Do not hide these choices in implementation details:

1. Should production use exact prescribed halos only, or exact halos plus an interior sponge?
2. What sponge width and timescale are appropriate for the PSKRIPS regional grid?
3. Which fields are prescribed for each forcing schema?
4. What initializes thermodynamics when a previously empty category gains ice under `structure_only` forcing?
5. Is velocity prescription required, and what is its source/staggering convention?
6. Is remap advection acceptable once exact and technical halos are handled, or should regional runs use upwind?
7. How should leap-day forcing be represented?

The standalone parent/nest experiments should provide the evidence for these decisions.

## Recommended first development issue

The first issue in the standalone repository should be narrowly scoped:

> Reproduce PSKRIPS time-varying regional ice restoring on the official `CICE6.6.3` release (`c734c12b`) with no semantic changes. Port the reader, minimal namelist integration, lifecycle finalization, converters/validator, and regional tests. Demonstrate disabled bit-for-bit behavior, all three forcing schemas, multi-rank interpolation across a record transition, and continuous/restart-split exactness. Do not implement PR #1110 exact-halo or sponge semantics in this issue.

The next issue should implement Stages 1-2 only: observability and API separation with bit-for-bit legacy results. Exact halo behavior then begins in a third, independently reviewable issue.

## References

- [CICE PR #1110: regional restoring and boundary work](https://github.com/CICE-Consortium/CICE/pull/1110)
- [PR #1110 restoring implementation at the reviewed head](https://github.com/CICE-Consortium/CICE/blob/1032b33f6f427b8dcd882afddc8a2ce3825999a9/cicecore/cicedyn/infrastructure/ice_restoring.F90)
- [PR #1110 user-guide changes](https://github.com/CICE-Consortium/CICE/blob/1032b33f6f427b8dcd882afddc8a2ce3825999a9/doc/source/user_guide/ug_implementation.rst)
- [PR #1110 full-to-nest test script](https://github.com/CICE-Consortium/CICE/blob/1032b33f6f427b8dcd882afddc8a2ce3825999a9/configuration/scripts/tests/full2nest.sh)
- [PR #1110 review discussion of remap boundary error](https://github.com/CICE-Consortium/CICE/pull/1110#discussion_r3306081942)
- [PR #1110 thermodynamic-restoring discussion](https://github.com/CICE-Consortium/CICE/pull/1110#issuecomment-4594615018)
- [Current PSKRIPS forcing and runtime guide](cice_time_varying_boundaries.md)
