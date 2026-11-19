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

    numberOfElements = size(elements, 1);
    elementJacobian = zeros(2, 2, numberOfElements);
    elementJacobianInverseTranspose = zeros(2, 2, numberOfElements);
    elementJacobianDeterminant = zeros(numberOfElements, 1);

    % Build the Jacobian first: J = [p2-p1, p3-p1].
    % Its determinant is positive for counter-clockwise triangles.
    for element = 1:numberOfElements
        elementNodes = nodes(elements(element, :), :);
        jacobian = [elementNodes(2, :) - elementNodes(1, :); ...
                    elementNodes(3, :) - elementNodes(1, :)].';
        determinant = det(jacobian);
        if ~isfinite(determinant) || determinant == 0
            error('load_fem_mesh:DegenerateElements', ...
                'All triangle elements must have finite, nonzero area.');
        end
        if determinant < 0
            % Swap the second and third nodes to reverse the triangle orientation.
            elements(element, [2, 3]) = elements(element, [3, 2]);
            elementNodes = nodes(elements(element, :), :);
            jacobian = [elementNodes(2, :) - elementNodes(1, :); ...
                        elementNodes(3, :) - elementNodes(1, :)].';
            determinant = det(jacobian);
        end
        elementJacobian(:, :, element) = jacobian;
        elementJacobianInverseTranspose(:, :, element) = inv(jacobian).';
        elementJacobianDeterminant(element) = determinant;
    end

    [edgeNodes, elementEdges, elementEdgeSigns] = ...
        build_edge_topology(elements);
    firstNodes = nodes(elements(:, 1), :);
    secondNodes = nodes(elements(:, 2), :);
    thirdNodes = nodes(elements(:, 3), :);
    elementEdgeLengths = [ ...
        hypot(secondNodes(:, 1) - firstNodes(:, 1), ...
              secondNodes(:, 2) - firstNodes(:, 2)), ...
        hypot(thirdNodes(:, 1) - secondNodes(:, 1), ...
              thirdNodes(:, 2) - secondNodes(:, 2)), ...
        hypot(firstNodes(:, 1) - thirdNodes(:, 1), ...
              firstNodes(:, 2) - thirdNodes(:, 2))];
    edgeDeltaX = [secondNodes(:, 1) - firstNodes(:, 1), ...
        thirdNodes(:, 1) - secondNodes(:, 1), ...
        firstNodes(:, 1) - thirdNodes(:, 1)];
    edgeDeltaY = [secondNodes(:, 2) - firstNodes(:, 2), ...
        thirdNodes(:, 2) - secondNodes(:, 2), ...
        firstNodes(:, 2) - thirdNodes(:, 2)];
    elementArea = abs(elementJacobianDeterminant) / 2;
    % For counter-clockwise local edges, the right-hand perpendicular is
    % the outward unit normal.
    elementEdgeNormals = cat(3, edgeDeltaY, -edgeDeltaX);
    elementEdgeNormals = elementEdgeNormals ./ ... % normalize
        reshape(elementEdgeLengths, [], 3, 1);
    if ~isequal(size(elementEdgeNormals), [numberOfElements, 3, 2])
        error('load_fem_mesh:InvalidElementNormals', ...
            'Element edge normals must have size [numberOfElements, 3, 2].');
    end
    if any(~isfinite(elementEdgeLengths(:))) || ...
            any(elementEdgeLengths(:) <= 0) || ...
            any(~isfinite(elementArea)) || any(elementArea <= 0)
        error('load_fem_mesh:DegenerateElements', ...
            'All triangle elements must have finite, positive edge lengths and area.');
    end

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
    mesh.elementEdgeLengths = elementEdgeLengths;
    mesh.elementArea = elementArea;
    mesh.elementJacobian = elementJacobian;
    mesh.elementJacobianInverseTranspose = elementJacobianInverseTranspose;
    mesh.elementJacobianDeterminant = elementJacobianDeterminant;
    mesh.elementEdgeNormals = elementEdgeNormals;
    mesh.elementNormals = elementEdgeNormals;
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
