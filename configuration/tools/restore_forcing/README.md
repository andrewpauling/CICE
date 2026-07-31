# CICE restore-forcing tools

These tools convert and validate the NetCDF inputs used by the standalone
time-varying regional ice-restoring implementation.

They were ported on 31 July 2026 from:

```text
/nfs/home/paulinan/models/coupled_model_apptainer/pskrips
```

The converters and environment retain the canonical implementation. The unit
test has a path-only import change. The validator's zero-area consistency check
was corrected during this port to distinguish a truly zero category from a
small positive thin-tail category; a regression test records that behavior.

Canonical SHA-256 checksums:

| Source | SHA-256 |
| --- | --- |
| `environment.cice-setup.yml` | `50b6cdda3b460d542f8cbb789b46ddce88e1b9cfd89d5adcfbe0d357437306b0` |
| `tools/cice_setup_python.sh` | `d03cabadc602750b83fb45a7f26c44f5d8ada0d4e244c69aafa5c324ce84ec16` |
| `tools/cice_thermo_utils.py` | `c681af5f1f8c3faa1a75f70d093915ca62dbe978ed1f8dab52e6e0106185dd70` |
| `tools/convert_cice_history_to_restore_forcing.py` | `47cd3afe5f83b00adc5dac98790676d7ee9669d8ed6c7af9883b02896365a079` |
| `tools/convert_climate_ice_to_restore_forcing.py` | `bf54173370646f1207a179f66c8667a7bc81382dbdb01f73688d5952a68957bc` |
| `tools/validate_cice_restore_forcing.py` | `1d895f35e630b2a2459d6c6a1aaded8bdd10f0d19a62e1f394861192f0911ea8` |
| `tests/test_cice_restore_tools.py` | `74f1832efb6ff2313489a6cc3f6a37a438a7728d51278903062e3d5f3e3ff724` |

Create the environment and run the tests from this directory:

```bash
mamba env create -f environment.yml
./cice_setup_python.sh -m unittest discover -s tests -p 'test_cice_restore_tools.py' -v
```

Validate a configured case:

```bash
./cice_setup_python.sh validate_cice_restore_forcing.py \
  --forcing-file /path/to/restore_forcing.nc \
  --ice-in /path/to/case/ice_in \
  --data-type daily_netcdf \
  --cycle-year
```

The existing `cice_ERA5_ECCO_2017_bc_fields_v1_2017.nc` file predates the
zero-based converter convention and stores `time=1..365`. Validate that file
with the additional `--allow-one-based-time` compatibility switch. The CICE
reader uses record position rather than the coordinate values.
