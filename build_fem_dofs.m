% BUILD_FEM_DOFS 
% Define the global unknown numbering for the hybrid system.
% Vector unknowns (Ex,Ey) are attached to global edges; scalar unknowns (Ez)
% are attached to global nodes. The two ranges define their positions in the
% concatenated global vector [E_t; E_z].

function dofs = build_fem_dofs(mesh)
    requiredFields = {'nodes', 'elements', 'edgeNodes', ...
        'elementEdges', 'elementEdgeSigns'};
    if ~isstruct(mesh) || ~all(isfield(mesh, requiredFields))
        error('build_fem_dofs:InvalidMesh', ...
            'mesh must be returned by load_fem_mesh.');
    end

    dofs.vectorCount = size(mesh.edgeNodes, 1);
    dofs.scalarCount = size(mesh.nodes, 1);
    dofs.totalCount = dofs.vectorCount + dofs.scalarCount;
    dofs.vectorElementIds = mesh.elementEdges;
    dofs.vectorElementSigns = mesh.elementEdgeSigns;
    dofs.scalarElementIds = mesh.elements;
    dofs.vectorRange = 1:dofs.vectorCount;
    dofs.scalarRange = dofs.vectorCount + (1:dofs.scalarCount);
end
