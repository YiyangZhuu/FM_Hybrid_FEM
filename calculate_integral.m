% CALCULATE_INTEGRAL 
% Return one block's local A and B coefficient matrices.
% Called once per triangle and block by assemble_hybrid_blocks. elementData
% contains nodes, coordinates, domain, edgeSigns, k0, and the element's
% relativePermittivity.

function [localA, localB] = calculate_integral(blockName, elementData)
    [points, weights] = fem_basis('quadrature', 2);
    [Nt, curlNt] = fem_basis('vector', elementData.coordinates, points);
    [Nz, gradNz] = fem_basis('nodal', elementData.coordinates, points);
    localA = zeros(3, 3);
    localB = zeros(3, 3);
    x = elementData.coordinates(:, 1);
    y = elementData.coordinates(:, 2);
    % det(J) = 2 * area, mapping to reference (0,0), (1,0), (0,1) triagle.
    doubleArea = (x(2) - x(1)) * (y(3) - y(1)) - (x(3) - x(1)) * (y(2) - y(1));

    % Atz, Azt, Azz are zeros
    switch upper(blockName)
        case 'KTT'
            localA = localA + doubleArea * (curlNt.' * (weights .* curlNt) - ...
                     elementData.k0^2 * elementData.relativePermittivity * ( ...
                     Nt(:,:,1).' * (weights .* Nt(:,:,1)) + ...
                     Nt(:,:,2).' * (weights .* Nt(:,:,2))));
            localB = localB +
        case 'KTZ'
            % TODO: Derive local A and B contributions for Ktz.
            % localB = localB +
        case 'KZT'
            % TODO: Derive local A and B contributions for Kzt independently.
            % localB = localB +
        case 'KZZ'
            % TODO: Derive local A and B contributions for Kzz.
            % localB = localB +
        otherwise
            error('calculate_integral:UnknownBlock', ...
                'Unknown matrix block "%s".', blockName);
    end
end
