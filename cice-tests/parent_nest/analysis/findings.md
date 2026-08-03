# Parent-to-nest remap/upwind findings

All four controlled five-day transport runs completed successfully.  The
parent and nest use the same P-SKRIPS grid geometry, an exact-index 84 x 80
subdomain, identical fixed C-grid velocities, and scheme-specific daily
boundary states taken from the matching parent run.

The final corrected comparison is:

| Scheme | AICE max error | AICE RMS error | Relative AICE integral error | AICE >= 0.15 mask disagreement |
| --- | ---: | ---: | ---: | ---: |
| remap | 6.82496e-1 | 6.59596e-2 | -2.95574e-2 | 5.20833e-3 |
| upwind | 3.13295e-2 | 1.23370e-3 | -1.51592e-4 | 0 |

The remap AICE RMS error is 53.46 times the upwind error.  Remap loses about
2.96% of integrated ice area over five days, versus 0.0152% for upwind.
Because category thickness and snow/ice ratios are fixed in this transport
test, VICE and VSNO have the same relative integral errors.

The imposed velocity is eastward and northward, so west and south are the
inflow edges.  Away from the four corners, final outer-ring AICE RMS errors
for west/east/south/north are:

- remap: 0.45946 / 0.01007 / 0.27811 / 0.00349
- upwind: 0.00794 / 0 / 0.00580 / 0

Remap error grows monotonically each day (AICE RMS 0.01772, 0.03243,
0.04490, 0.05595, 0.06596).  It remains detectable through ring 9 from the
nearest edge and is exactly zero from ring 10 inward after this five-step
test.  Upwind grows much more slowly (0.000325 to 0.001234), is zero from
ring 5 inward, and does not move the 0.15 concentration contour.

The prescribed-state halo machinery itself passes its assertions in both
nested runs: `physical_max_delta`, `strict_interior_max_delta`,
`boundary_max_delta`, and `extension_max_error` are all zero.  The extended
nested restart also makes `uvel`, `vvel`, `uvelE`, and `vvelN` identical to
the matching parent stencil, including the exterior velocity cells needed
when standalone dynamics is disabled.

This isolates the remaining problem to the remap path rather than the forcing
reader or state-halo reapplication.  Inside `horizontal_remap`, CICE derives
departure points and reconstruction moments from the prescribed state, then
performs native halo exchanges on `dpx/dpy`, `mc/mx/my`, and `tc/tx/ty`.
Those technical arrays are not covered by `ice_HaloRestore_apply_halo`.
Consequently, an open-boundary exchange can replace the derived exterior
values before remap consumes them.  Upwind consumes the prescribed state and
face velocities directly and therefore largely reproduces the parent.

The next implementation test should preserve or reconstruct the remap
technical halos on enabled prescribed boundaries immediately after those
internal halo exchanges, then rerun this same four-case experiment.  Success
criteria should be an upwind-scale remap error, no 0.15-mask disagreement,
and no degradation of the exact-halo assertions.
