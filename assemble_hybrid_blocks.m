% ASSEMBLE_HYBRID_BLOCKS
% Assemble the two matrices in A*x = lambda*B*x.
% blocks.A and blocks.B contain their respective Ktt, Ktz, Kzt, and Kzz
% sparse matrices. The two coupling blocks are assembled independently.

function [A, B, blocks] = assemble_hybrid_blocks(mesh, dofs)
    numberOfEdges = dofs.vectorCount;
    numberOfNodes = dofs.scalarCount;
    blocks.A = struct();
    blocks.B = struct();

    [blocks.A.Ktt, blocks.B.Ktt] = assemble_block( ...
        mesh, dofs, 'Ktt', 'vector', 'vector', numberOfEdges, numberOfEdges);
    [blocks.A.Ktz, blocks.B.Ktz] = assemble_block( ...
        mesh, dofs, 'Ktz', 'vector', 'scalar', numberOfEdges, numberOfNodes);
    [blocks.A.Kzt, blocks.B.Kzt] = assemble_block( ...
        mesh, dofs, 'Kzt', 'scalar', 'vector', numberOfNodes, numberOfEdges);
    [blocks.A.Kzz, blocks.B.Kzz] = assemble_block( ...
        mesh, dofs, 'Kzz', 'scalar', 'scalar', numberOfNodes, numberOfNodes);

    A = [blocks.A.Ktt, blocks.A.Ktz; ...
        blocks.A.Kzt, blocks.A.Kzz];
    B = [blocks.B.Ktt, blocks.B.Ktz; ...
        blocks.B.Kzt, blocks.B.Kzz];
end

function [matrixA, matrixB] = assemble_block(mesh, dofs, blockName, ...
        rowSpace, columnSpace, nRows, nColumns)
    numberOfElements = size(mesh.elements, 1);
    entriesPerElement = 9;
    % global indices of the local contributions to the global matrix
    rowIndices = zeros(numberOfElements * entriesPerElement, 1);
    columnIndices = zeros(numberOfElements * entriesPerElement, 1);
    % matrix entries of the local contributions to the global matrix
    valuesA = zeros(numberOfElements * entriesPerElement, 1);
    valuesB = zeros(numberOfElements * entriesPerElement, 1);
    entry = 0;

    for element = 1:numberOfElements
        elementData.nodes = mesh.elements(element, :);
        elementData.coordinates = mesh.nodes(elementData.nodes, :);
        elementData.domain = mesh.domains(element);
        elementData.edgeSigns = mesh.elementEdgeSigns(element, :);
        elementData.elementEdgeLengths = mesh.elementEdgeLengths(element, :);
        elementData.elementArea = mesh.elementArea(element);
        if isfield(mesh, 'relativePermittivity')
            elementData.relativePermittivity = mesh.relativePermittivity(element);
        end
        if isfield(mesh, 'k0')
            elementData.k0 = mesh.k0;
        end

        [localRows, rowSigns] = element_dofs(dofs, element, rowSpace);
        [localColumns, columnSigns] = ...
            element_dofs(dofs, element, columnSpace);
        [localA, localB] = calculate_integral(blockName, elementData);
        validate_local_matrix(localA, blockName, 'A');
        validate_local_matrix(localB, blockName, 'B');

        orientation = rowSigns(:) * columnSigns(:).';
        localA = localA .* orientation;
        localB = localB .* orientation;
        [rowGrid, columnGrid] = ndgrid(localRows, localColumns);
        indices = entry + (1:entriesPerElement);
        rowIndices(indices) = rowGrid(:);
        columnIndices(indices) = columnGrid(:);
        valuesA(indices) = localA(:);
        valuesB(indices) = localB(:);
        entry = entry + entriesPerElement;
    end

    matrixA = sparse(rowIndices, columnIndices, valuesA, nRows, nColumns);
    matrixB = sparse(rowIndices, columnIndices, valuesB, nRows, nColumns);
end

function validate_local_matrix(localMatrix, blockName, coefficientName)
    if ~isnumeric(localMatrix) || ~isequal(size(localMatrix), [3, 3]) || ...
            any(~isfinite(localMatrix(:)))
        error('assemble_hybrid_blocks:InvalidElementMatrix', ...
            ['The %s contribution for %s must be a finite numeric ', ...
            '3-by-3 matrix.'], coefficientName, blockName);
    end
end

function [ids, signs] = element_dofs(dofs, element, space)
    switch space
        case 'vector'
            ids = dofs.vectorElementIds(element, :);
            signs = dofs.vectorElementSigns(element, :);
        case 'scalar'
            ids = dofs.scalarElementIds(element, :);
            signs = ones(1, 3);
        otherwise
            error('assemble_hybrid_blocks:UnknownSpace', ...
                'Unknown local DOF space "%s".', space);
    end
end
