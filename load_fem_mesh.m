% LOAD_FEM_MESH
% Load triangle data and construct global edge topology.

function mesh = load_fem_mesh(meshDirectory)
    nodeFile = fullfile(meshDirectory, 'nodes.txt');
    elementFile = fullfile(meshDirectory, 'elements.txt');
    domainFile = fullfile(meshDirectory, 'domains.txt');
    boundaryFile = fullfile(meshDirectory, 'bcs.txt');

    nodes = readmatrix(nodeFile);
    elements = readmatrix(elementFile) + 1; % MATLAB is 1-based
    domains = readmatrix(domainFile);
    boundaryEdges = readmatrix(boundaryFile) + 1;

    % Some assertions
    if size(nodes, 2) ~= 2 || isempty(nodes) || any(~isfinite(nodes(:)))
        error('load_fem_mesh:InvalidNodes', ...
            'nodes.txt must contain finite x-y pairs.');
    end
    if size(elements, 2) ~= 3 || isempty(elements) || ...
            any(elements(:) < 1) || any(elements(:) > size(nodes, 1)) || ...
            any(elements(:) ~= fix(elements(:)))
        error('load_fem_mesh:InvalidElements', ...
            'elements.txt must contain valid three-node triangle connectivity.');
    end
    if numel(domains) ~= size(elements, 1)
        error('load_fem_mesh:InvalidDomains', ...
            'domains.txt must contain one material/domain ID per element.');
    end
    if size(boundaryEdges, 2) ~= 2 || ...
            any(boundaryEdges(:) < 1) || any(boundaryEdges(:) > size(nodes, 1)) || ...
            any(boundaryEdges(:) ~= fix(boundaryEdges(:)))
        error('load_fem_mesh:InvalidBoundaryEdges', ...
            'bcs.txt must contain valid two-node edge connectivity.');
    end


    [edgeNodes, elementEdges, elementEdgeSigns] = ...
        build_edge_topology(elements);
    sortedBCEdges = sort(boundaryEdges, 2);
    [isMeshEdge, bcEdgeIds] = ismember(sortedBCEdges, edgeNodes, 'rows');
    if any(~isMeshEdge)
        error('load_fem_mesh:InvalidBCEdges', ...
            'At least one edge in bcs.txt is not present in the triangle mesh.');
    end

    mesh.nodes = nodes;
    mesh.elements = elements;
    mesh.domains = domains(:);
    mesh.bcEdges = boundaryEdges;
    mesh.bcEdgeIds = bcEdgeIds;
    mesh.edgeNodes = edgeNodes;
    mesh.elementEdges = elementEdges;
    mesh.elementEdgeSigns = elementEdgeSigns;
end

function [edgeNodes, elementEdges, elementEdgeSigns] = ...
        build_edge_topology(elements)
    localEdgeNodes = [elements(:, [1, 2]); ...
                    elements(:, [2, 3]); ...
                    elements(:, [3, 1])];
    sortedLocalEdgesNodes = sort(localEdgeNodes, 2);
    [edgeNodes, ~, edgeIds] = unique(sortedLocalEdgesNodes, 'rows');
    numberOfElements = size(elements, 1);
    elementEdges = [edgeIds(1:numberOfElements), ...
                    edgeIds(numberOfElements + 1:2 * numberOfElements), ...
                    edgeIds(2 * numberOfElements + 1:3 * numberOfElements)];
    localSigns = 2 * (localEdgeNodes(:, 1) < localEdgeNodes(:, 2)) - 1;
    elementEdgeSigns = [localSigns(1:numberOfElements), ...
                        localSigns(numberOfElements + 1:2 * numberOfElements), ...
                        localSigns(2 * numberOfElements + 1:3 * numberOfElements)];
end
