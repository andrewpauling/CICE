# P-SKRIPS ERA5 standalone baseline

This manifest records the successful boundary-disabled regional run that
preceded the standalone time-varying boundary-condition port.

## Provenance

- Run date: 31 July 2026
- Slurm job: `3319512`
- Case: `pskrips_era5_test`
- CICE version: `CICE_6.6.3`
- CICE release commit: `c734c12bfb19632cd700cb06b22c3ea62b4adeae`
- Icepack commit: `daa41638c6cef298583f35ee4dec5ab1bb077aea`
- Machine port commit: `2176b3c596c5188573fd7646382a220246b5ec6a`
- Build-relevant source state: `8654227` (`Add ERA5 forcing support for P-SKRIPS`)
- Simulation interval: 2005-01-01 00:00:00 through 2005-01-02 00:00:00
- Timestep: 3600 seconds
- MPI tasks: 30
- Time-varying ice boundary restoring: not yet implemented

The case was created with:

```text
./cice.setup -c /nfs/home/paulinan/cice-dirs/cases/pskrips_era5_test \
  -s regional -g gx1 -m raapoi -e intel
```

The run log ends with `CICE COMPLETED SUCCESSFULLY`.

## SHA-256 checksums

### Case and executable

```text
6fe4ffe3d6c4e095efa9d7faa431156fdbfbbfbdda82b76505b5de7325e39564  /nfs/home/paulinan/cice-dirs/cases/pskrips_era5_test/ice_in
20d621cc4ef58e3eea1c5c09f7adc643c9e0c26c5985e5c8f243c4e27a23b477  /nfs/home/paulinan/cice-dirs/cases/pskrips_era5_test/cice.settings
cdc83f4e5b913879f8171e121cdb902bf9acfc58a0ea1692b78c9fbfd8d50d5b  /nfs/home/paulinan/cice-dirs/cases/pskrips_era5_test/env.raapoi_intel
db7612436f0e5b600cb555da011aa7e0992abf69e0b8fab6934a5bf513ec634d  /nfs/home/paulinan/cice-dirs/cases/pskrips_era5_test/logs/cice.runlog.260731-142749
f4799c1eaa130a56476cb0ed8e1ff1e7564240e6f127d8a30bbd63c0caf3eb3a  /nfs/scratch/paulinan/cice-dirs/runs/pskrips_era5_test/cice
```

### Output

```text
1d4ed13a5e7160dc94a2b328f7b07ebc28cdcf017672994f73edf49c98b33c60  /nfs/scratch/paulinan/cice-dirs/runs/pskrips_era5_test/history/iceh_ic.2005-01-01-00000.nc
dcc1ec2042c16ddd0a2838f1302266d38257445b1f3fe96cda177d26ef59f680  /nfs/scratch/paulinan/cice-dirs/runs/pskrips_era5_test/restart/iced.2005-01-02-00000.nc
28e81e394f16663f7b77556fe0d37febaf696615c8b54f67e53e0e3a82233845  /nfs/scratch/paulinan/cice-dirs/runs/pskrips_era5_test/restart/iceh_rd.2005-01-02-00000.nc
60016973766dc46f1d41372d729632dc1a4e7294a971875283be5bd02104a182  /nfs/scratch/paulinan/cice-dirs/runs/pskrips_era5_test/restart/iceh_rm.2005-01-02-00000.nc
```

### External input

```text
a5019e17723c66612bb4193331438ff77d51277d9a49d61b4c5016ace013abe7  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/grid_pskrips.nc
a145cbce791d374e869e9a7561ff4d18fb22f8cf8ae69e3aba4579f72e760b90  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/kmt_pskrips.nc
5a519d70440fbfe8aaa0267770f4a590de4c0503024e80271701a66722a076ec  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/bath_pskrips.nc
10fee89f0e78e7e44828469a3f929bc2d86707abfe0fc97b445c39bfc5a0791b  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/ic/pskrips/iced_pskrips_v6.2005-01-01.nc
8d621449bb2bab4828abc98dc946c77222f1cb4bed44f7900219b9c791853d22  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/forcing/pskrips/ERA5/8XDAILY/ERA5_03hr_forcing_2005.nc
c83521dbeeabe74fa4d645468d39ec00da2d78c8ba844536a037a26cfe28bbc5  /nfs/scratch/paulinan/cice-dirs/input/CICE_data/forcing/pskrips/CESM/MONTHLY/ocean_forcing_clim_2D_pskrips.nc
```

Binary NetCDF hashes preserve the exact original artifacts. For later
disabled-mode regression testing, also compare all NetCDF variable values
while ignoring creation-time metadata.
