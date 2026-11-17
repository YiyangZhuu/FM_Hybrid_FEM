% CALCULATE_INTEGRAL 
% Return one block's local A and B coefficient matrices.
% Called once per triangle and block by assemble_hybrid_blocks. elementData
% contains coordinates, edge lengths, area, k0, and relativePermittivity.

function [localA, localB] = calculate_integral(blockName, elementData)
    [points, weights] = fem_basis('quadrature', 2);
    [Nt, curlNt] = fem_basis('vector', elementData.coordinates, points);
    [Nz, gradNz] = fem_basis('nodal', elementData.coordinates, points);
    localA = zeros(3, 3);
    localB = zeros(3, 3);
    if ~isnumeric(elementData.elementEdgeLengths) || ...
            ~isequal(size(elementData.elementEdgeLengths), [1, 3]) || ...
            any(~isfinite(elementData.elementEdgeLengths)) || ...
            any(elementData.elementEdgeLengths <= 0) || ...
            ~isnumeric(elementData.elementArea) || ...
            ~isscalar(elementData.elementArea) || ...
            ~isfinite(elementData.elementArea) || elementData.elementArea <= 0
        error('calculate_integral:InvalidElementGeometry', ...
            'Element edge lengths and area must be finite and positive.');
    end

    % The vector basis has unit edge-circulation normalization; scale each
    % local edge function by its edge length for tangential-field DOFs.
    edgeLengths = elementData.elementEdgeLengths;
    Nt = bsxfun(@times, Nt, reshape(edgeLengths, 1, 3, 1));
    curlNt = bsxfun(@times, curlNt, edgeLengths);
    jacobianDeterminant = 2 * elementData.elementArea;
    vectorMass = Nt(:, :, 1).' * (weights .* Nt(:, :, 1)) + ...
        Nt(:, :, 2).' * (weights .* Nt(:, :, 2));
    scalarMass = Nz.' * (weights .* Nz);
    integratedNtX = Nt(:, :, 1).' * weights;
    integratedNtY = Nt(:, :, 2).' * weights;

    switch upper(blockName)
        case 'KTT'
            localA = jacobianDeterminant * ( ...
                curlNt.' * (weights .* curlNt) - ...
                elementData.k0^2 * elementData.relativePermittivity * vectorMass);
            localB = jacobianDeterminant * vectorMass;
        case 'KTZ'
            localB = -jacobianDeterminant * ( ...
                integratedNtX * gradNz(:, 1).' + ...
                integratedNtY * gradNz(:, 2).');
        case 'KZT'
            % Where is grad Nz from?? there is no gradient. Check last page
            % I'm not sure abt this
            % localB = -jacobianDeterminant * ( ...
            %     gradNz(:, 1) * integratedNtX.' + ...
            %     gradNz(:, 2) * integratedNtY.');
        case 'KZZ'
            localB = jacobianDeterminant * ...
                elementData.k0^2 * elementData.relativePermittivity * scalarMass - ...
                elementData.elementArea * (gradNz * gradNz.');
        otherwise
            error('calculate_integral:UnknownBlock', ...
                'Unknown matrix block "%s".', blockName);
    end
end
