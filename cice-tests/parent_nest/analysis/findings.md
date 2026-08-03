# Parent-to-nest remap/upwind findings

All four controlled five-day transport runs completed successfully.  The
parent and nest use the same P-SKRIPS grid geometry, an exact-index 84 x 80
subdomain, identical fixed C-grid velocities, and scheme-specific daily
boundary states taken from the matching parent run.

Before reconstructing remap's scheme-internal open-boundary halos, the final
comparison was:

| Scheme | AICE max error | AICE RMS error | Relative AICE integral error | AICE >= 0.15 mask disagreement |
| --- | ---: | ---: | ---: | ---: |
| remap | 6.82496e-1 | 6.59596e-2 | -2.95574e-2 | 5.20833e-3 |
| upwind | 3.13295e-2 | 1.23370e-3 | -1.51592e-4 | 0 |

The remap AICE RMS error was 53.46 times the upwind error.  Remap lost about
2.96% of integrated ice area over five days, versus 0.0152% for upwind.
Because category thickness and snow/ice ratios are fixed in this transport
test, VICE and VSNO have the same relative integral errors.

The imposed velocity is eastward and northward, so west and south are the
inflow edges.  Away from the four corners, final outer-ring AICE RMS errors
for west/east/south/north are:

- remap: 0.45946 / 0.01007 / 0.27811 / 0.00349
- upwind: 0.00794 / 0 / 0.00580 / 0

The original remap error grew monotonically each day (AICE RMS 0.01772, 0.03243,
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

This isolated the problem to the remap path rather than the forcing reader or
state-halo reapplication.  With CICE's compile-time `nghost=1`,
`make_masks`, `construct_fields`, and `departure_points` deliberately operate
only on physical cells.  The resulting `dpx/dpy`, `mc/mx/my`, and `tc/tx/ty`
arrays are initialized to zero outside the physical domain.  Native halo
exchange then fills inter-block halos but correctly leaves global open-edge
halos unchanged.  Remap therefore consumed zero-valued technical halos at
the regional exterior even though the corresponding physical state and
velocity halos were prescribed correctly.  Upwind consumes those physical
fields directly and did not have this failure mode.

The implementation now reconstructs those remap-only arrays in memory after
their native halo exchange.  It is active only for enabled prescribed edges
on open boundaries.  The forcing interface remains physical: cell means come
from the prescribed state, departure points use the prescribed velocity, and
subcell gradients in the single available exterior layer are zero.  No
remap-specific fields are added to the boundary NetCDF format.

The complete four-case test was rebuilt and rerun.  Both parent controls and
the nested upwind sequence are bit-for-bit unchanged.  The post-fix final
comparison is:

| Scheme | AICE max error | AICE RMS error | Relative AICE integral error | AICE >= 0.15 mask disagreement |
| --- | ---: | ---: | ---: | ---: |
| remap | 2.52638e-2 | 2.43060e-3 | -5.03820e-4 | 2.97619e-4 |
| upwind | 3.13295e-2 | 1.23370e-3 | -1.51592e-4 | 0 |

Remap AICE RMS is reduced by 96.3%, the magnitude of its integral error by
98.3%, and the remap/upwind RMS ratio from 53.46 to 1.97.  Its daily AICE RMS
is now 0.000537, 0.001035, 0.001511, 0.001974, and 0.002431.  Final remap
outer-ring RMS for west/east/south/north is
0.01428 / 0.01007 / 0.00745 / 0.00349.  Thus the catastrophic inflow-edge
error is removed; the remaining boundary-local difference reflects the
first-order reconstruction forced by having one physical exterior layer,
whereas the parent can form its native limited gradients and midpoint
velocity from a wider stencil.

This result supports keeping the upstream-style physical boundary interface
and reconstructing scheme-specific quantities internally.  Exact
parent-to-nest remap equivalence would require a wider prescribed physical
stencil (or a general increase in `nghost`), not saved remap intermediates.
