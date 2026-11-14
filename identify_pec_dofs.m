% IDENTIFY_PEC_DOFS
% Find hybrid field DOFs on the outer circular PEC boundary.
% Only the boundary edges listed in mesh.bcEdges are considered. A listed
% edge is PEC when both endpoints lie within tolerance of outerRadius.
% Scalar unknowns are selected at the endpoints of those PEC edges.

function pec = identify_pec_dofs(mesh, dofs, outerRadius, tolerance)
    nodeRadius = hypot(mesh.nodes(:, 1), mesh.nodes(:, 2));
    isOuterNode = abs(nodeRadius - outerRadius) <= tolerance;
    pec.bcRecordIds = find(all(isOuterNode(mesh.bcEdges), 2));
    pec.vectorEdgeIds = unique(mesh.bcEdgeIds(pec.bcRecordIds), 'stable');
    pec.scalarNodeIds = unique( ...
        mesh.bcEdges(pec.bcRecordIds, :).', 'stable');
    pec.vectorGlobalIds = pec.vectorEdgeIds;
    pec.scalarGlobalIds = dofs.vectorCount + pec.scalarNodeIds;
    pec.constrainedGlobalIds = [pec.vectorGlobalIds; pec.scalarGlobalIds];
    pec.freeGlobalIds = setdiff((1:dofs.totalCount).', ...
        pec.constrainedGlobalIds, 'stable');
    pec.outerRadius = outerRadius;
    pec.tolerance = tolerance;
end
