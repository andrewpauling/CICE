# Standalone P-SKRIPS parent-to-nest transport test

This experiment compares a nested regional CICE run with an exact-index
subset of a larger parent run.  It is intentionally transport-only:

- the parent uses the full 210 x 240 P-SKRIPS grid;
- the nest is the fully wet 84 x 80 region `x=100:184, y=140:220` (Python
  zero-based, upper bound exclusive);
- thermodynamics, momentum dynamics, and ridging are disabled;
- a fixed eastward/northward velocity transports a smooth synthetic ice
  field;
- the nest reads a one-cell extended restart so its fixed C-grid face
  velocities match the parent outside the west and south edges (there is no
  dynamics step to reconstruct velocity halos in this transport-only test);
- parent and nest are run separately with `remap` and `upwind`;
- nested boundary records come directly from the matching parent restart
  files with no spatial interpolation; because `structure_only` projects a
  file-edge value into the adjacent ghost cell, each file edge is populated
  from the corresponding parent cell immediately outside the nest.

The workflow is:

1. `generate_inputs.py` creates the synthetic parent restart and exact nest
   grid/restart files.
2. `configure_cases.py` modifies four cases created by `cice.setup`.
3. Run both parent cases.
4. `make_boundary_forcing.py` creates scheme-specific six-record daily
   boundary files from the initial state and five parent restart outputs.
5. Run both nested cases.
6. `analyze.py` compares each final nested restart with the exact parent
   subset and writes `analysis/summary.json` and `analysis/summary.txt`.

The Python environment used on Raapoi is
`/nfs/home/paulinan/miniforge3/envs/cice_setup/bin/python`.
