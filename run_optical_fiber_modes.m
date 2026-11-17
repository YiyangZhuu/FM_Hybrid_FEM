%% Physical and numerical parameters
k0 = 25;                       % Free-space wavenumber [1/m] between [18, 25]
epsilonR.core = 2.16;            % Core relative permittivity
epsilonR.cladding = 2.1;        % Cladding relative permittivity
epsilonR.air = 1.0;             % Air relative permittivity

% Domain IDs are inferred from the supplied mesh geometry:
% domain 3 is central, 2 is the surrounding ring, and 1 is the outer region.
materialDomain.core = 3;
materialDomain.cladding = 2;
materialDomain.air = 1;

outerRadius = 0.5;              % PEC boundary radius in mesh coordinate units
boundaryTolerance = 1e-6;      % Mesh nodes are within 7e-8 of r = 0.5

%% Mesh loading and material assignment

projectDirectory = fileparts(mfilename('fullpath'));
addpath(projectDirectory);

meshDirectory = fullfile(projectDirectory, 'Mesh_Files');
mesh = load_fem_mesh(meshDirectory);
dofs = build_fem_dofs(mesh);

mesh.k0 = k0;
mesh.relativePermittivity = nan(size(mesh.domains));
mesh.relativePermittivity(mesh.domains == materialDomain.core) = epsilonR.core;
mesh.relativePermittivity(mesh.domains == materialDomain.cladding) = epsilonR.cladding;
mesh.relativePermittivity(mesh.domains == materialDomain.air) = epsilonR.air;

%% PEC boundary-condition identification and application
% Identify only the circular exterior boundary. The identified global DOFs
% will be eliminated from the assembled system
pec = identify_pec_dofs(mesh, dofs, outerRadius, boundaryTolerance);
freeDofIds = pec.freeGlobalIds;

fprintf('Mesh: %d nodes, %d triangles, %d edges.\n', ...
    size(mesh.nodes, 1), size(mesh.elements, 1), size(mesh.edgeNodes, 1));
fprintf('Outer PEC: %d transverse edge DOFs, %d longitudinal nodal DOFs.\n', ...
    numel(pec.vectorEdgeIds), numel(pec.scalarNodeIds));

%% Matrix assembly
[A, B, blocks] = assemble_hybrid_blocks(mesh, dofs);

%% Generalized eigenvalue solver
% TODO: After deriving K(beta) = A - lambda*B with lambda = -beta^2,
% eliminate the PEC DOFs from both assembled coefficient matrices:
%
% Afree = A(freeDofIds, freeDofIds);
% Bfree = B(freeDofIds, freeDofIds);
% sigma = -(k0 * sqrt(epsilonR.cladding))^2;
% options.lossless = true;
% modes = solve_beta_modes(Afree, Bfree, sigma, numberOfModes, options);
%
% solve_beta_modes uses standard shift-invert on (A - sigma*B)\B, maps
% lambda = sigma + 1/nu, and checks each residual in the original pencil.
% This call must wait until the derivation identifies A and B.

%% Mode post-processing and visualization
% TODO: Reconstruct full fields by inserting zeros at PEC-constrained DOFs.
% TODO: Separate transverse edge coefficients and longitudinal nodal values.
% TODO: Plot mode profiles and report propagation constants.
