# Hybrid 2D FEM for Optical-Fiber Modes

This is a project for the course "_Physical modeling and simulation_". Topic is "2D Eigenvalue Analysis of Optical Fiber (Vector FEM)". The goal is first derive the theoretical description of the 2D FEM electromagnetic vector eigenvalue problem, implementing the hybrid vector-nodal FEM solver, computing and visualization of the first six modes. Potential extension is the analysis of triangular-core hollow optical fiber geometry (lossy case).

## Geometris

The provided mesh in `Mesh_Files` is read as zero-based triangle
and `bcs.txt` edge connectivity, then converted to MATLAB's one-based indices.
The edges in `bcs.txt` are mapped to mesh edges without assuming that all are
on the exterior boundary; duplicate records are preserved.

## Files

- `run_optical_fiber_modes.m` is the main workflow, organized into physical
  and numerical parameters, mesh/material assignment, PEC identification,
  assembly, eigenvalue solving, and post-processing sections.
- `load_fem_mesh.m` loads the mesh and constructs canonical global edges,
  element-to-edge maps, edge orientation signs, and topological exterior
  edges. Exterior edges are found by triangle adjacency, not by assuming
  every record in `bcs.txt` is an exterior edge.
- `build_fem_dofs.m` numbers the unknowns: one transverse-vector DOF per
  global edge and one longitudinal-scalar DOF per node. It records local-to-
  global maps and the two index ranges in the combined vector `[E_t; E_z]`.
- `identify_pec_dofs.m` checks only the edge records in `mesh.bcEdges` and
  identifies transverse edge and longitudinal nodal DOFs for those whose
  endpoints lie on the outer circle. It also returns the unconstrained global
  indices needed to eliminate PEC unknowns from the system.
- `fem_basis.m` provides first-order Nedelec edge functions, their scalar
  curls, first-order Lagrange shape functions and gradients, and degree-1/2
  triangle quadrature rules.
- `calculate_integral.m` is the local element-matrix template. Given a block
  and triangle data, it returns the local contributions to `A` and `B`.
  Each block has its own TODO; the function errors until they are implemented.
- `assemble_hybrid_blocks.m` calls `calculate_integral` for each element and
  assembles the four independent blocks of both global matrices `A` and `B`.
  Edge orientation is accounted for in every block involving vector DOFs.
- `solve_beta_modes.m` solves a generalized eigenvalue pencil `A*x =
  lambda*B*x` near a shift using standard shift-invert. It does not derive
  or assemble `A` and `B`; those must come from the weak-form derivation.

## Remaining mathematical work

Implement the four cases in `calculate_integral.m`. Each case must return two
3-by-3 local matrices, the element contributions to `A` and `B` in the
derived pencil `A*x = lambda*B*x`, where `lambda = -beta^2`. `Ktz` and `Kzt` are intentionally independent.

## Generalized eigenvalue solver

Once the weak form is known to be linear in `beta^2`, identify matrices `A`
and `B` such that `A*x = lambda*B*x` with `lambda = -beta^2`. Eliminate PEC
DOFs from both matrices, then call:

```matlab
options.lossless = true;
modes = solve_beta_modes(Afree, Bfree, sigma, numberOfModes, options);
```

The matrices are assembled directly once with
`[A, B, blocks] = assemble_hybrid_blocks(mesh, dofs)`.

The solver applies the standard operator `(A - sigma*B)^(-1)*B` through
linear solves, asks `eigs` for its largest-magnitude eigenvalues `nu`, and
maps them back with `lambda = sigma + 1/nu`. It deliberately does not use
ARPACK's generalized mode with `B`, so `B` may be symmetric indefinite. Every
returned mode is checked using the relative residual of the original pencil.
For lossless materials, `options.lossless = true` additionally enforces that
`beta^2 = -lambda` is real to the configured tolerance. The suggested initial
shift is `sigma = -(k0*n_clad)^2`, with `n_clad = sqrt(epsilonR.cladding)`.

The mesh's `bcs.txt` records are loaded and mapped. PEC identification filters
only those listed edges by the outer-radius test, leaving room to identify
other boundary-condition types separately later. PEC unknowns are identified
but not yet applied to a matrix because assembly remains unimplemented.

## Global parameters

`k0` and the relative permittivities are set in the main script and must match
the optical-fiber specification. The mesh domain IDs are assigned as core =
3, cladding = 2, and air = 1 based on their radial locations. The PEC radius
is 0.5 mesh units, with a 1e-6 tolerance.