```@meta
CurrentModule = PhaseSpaceDTFE
```

# Migrating from version 1

Version 2 keeps the point-query interface of version 1 (`SimBox`, `PS_DTFE_periodic`, `density`, `numberOfStreams`, `velocity`, `velocitySum`, `inSimplices`). The estimator itself is new, and some functions have been removed.

## Unchanged calls

The constructors of version 1 with the depth of the search tree are accepted; the depth is not used:

```julia
ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, vels, m, depth, sim_box)
ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, m, depth, sim_box)
```

These are equal to `PS_DTFE_periodic(coords_q, coords_x, vels, m, sim_box)` and `PS_DTFE(coords_q, coords_x, vels, m, sim_box)`.

## Removed functions

| Version 1 | Version 2 |
|---|---|
| `DTFE_periodic(coords_x, vels, m, depth, sim_box)` | `PS_DTFE(coords_x, coords_x, vels, m, sim_box)` |
| `DTFE_periodic(coords_x, m, depth, sim_box)` | `PS_DTFE(coords_x, coords_x, m, sim_box)` |
| `ps_dtfe_subbox(coords_q, coords_x, [vels,] m, depth, sim_box)` | `PS_DTFE(coords_q, coords_x, [vels,] m, sim_box)` for simulations whose estimator fits in memory; `ps_dtfe_grid` otherwise |
| `density_subbox(coords_arr, ps_dtfe_sb)` | `density.(coords_arr, Ref(ps_dtfe))`, or `density_grid(xs, ys, zs, ps_dtfe)` on a regular grid |
| `numberOfStreams_subbox(coords_arr, ps_dtfe_sb)` | `numberOfStreams.(coords_arr, Ref(ps_dtfe))`, or `numberOfStreams_grid(xs, ys, zs, ps_dtfe)` |
| `velocity_subbox`, `velocitySum_subbox` | `velocity.(coords_arr, Ref(ps_dtfe))`, `velocitySum.(coords_arr, Ref(ps_dtfe))` |

The estimator no longer writes files to disk (the `ps_dtfe` directory of the sub-box estimators of version 1).

## Changed results

* **Velocities.** Within a tetrahedron, the velocity is interpolated linearly with the barycentric coordinates of the point. Version 1 applied the transposed velocity gradient, which agrees only for irrotational, single-stream flow. Single-stream velocities of the 64³ test simulation differ by a median of 7% of the rms velocity.
* **Densities.** With the default `tessellation = :delaunay, vertex_density = :mass`, the estimator uses the same tessellation and vertex densities as version 1. With `jitter = 0`, the densities at 20,000 random points of the 64³ test simulation agree to 10⁻¹² at 99.4% of the points, and the numbers of streams agree everywhere. The points that differ lie within half a particle spacing of a face of the box, where version 1 tessellated translated copies of the particles instead of computing a periodic tessellation. With the default `jitter = 1e-5` the agreement is 98.3%; the jitter decides how nearly degenerate (almost cubic) cells of an initial lattice are split.
* **Saved estimators.** The fields of `PS_DTFE` have changed. Estimators of version 1 saved to file cannot be loaded; construct them again.
