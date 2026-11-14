% CALCULATE_INTEGRAL 
% Return one block's local A and B coefficient matrices.
% Called once per triangle and block by assemble_hybrid_blocks. elementData
% contains nodes, coordinates, domain, edgeSigns, k0, and the element's
% relativePermittivity.

function [localA, localB] = calculate_integral(blockName, elementData)
    switch upper(blockName)
        case 'KTT'
            % TODO: Derive local A and B contributions for Ktt.
        case 'KTZ'
            % TODO: Derive local A and B contributions for Ktz.
        case 'KZT'
            % TODO: Derive local A and B contributions for Kzt independently.
        case 'KZZ'
            % TODO: Derive local A and B contributions for Kzz.
        otherwise
            error('calculate_integral:UnknownBlock', ...
                'Unknown matrix block "%s".', blockName);
    end

    error('calculate_integral:NotImplemented', ...
        'The %s element coefficient matrices have not been implemented.', ...
        upper(blockName));
end
