% FEM_BASIS 
% First-order triangle bases and a basic reference-cell rule.
%   [N, curlN] = fem_basis('vector', coordinates, referencePoints)
%   [N, gradN] = fem_basis('nodal', coordinates, referencePoints)
%   [points, weights] = fem_basis('quadrature', polynomialDegree)
%
% Vector functions are ordered on local directed edges (1->2, 2->3, 3->1).
% Apply mesh.elementEdgeSigns when mapping local edge DOFs to global DOFs.
% Quadrature weights integrate on the reference triangle (area 1/2).

function varargout = fem_basis(kind, varargin)
    switch lower(kind)
        case 'vector'
            [varargout{1}, varargout{2}] = ...
                vector_basis(varargin{1}, varargin{2});
        case 'nodal'
            [varargout{1}, varargout{2}] = ...
                nodal_basis(varargin{1}, varargin{2});
        case 'quadrature'
            [varargout{1}, varargout{2}] = triangle_rule(varargin{1});
        otherwise
            error('fem_basis:UnknownKind', 'Unknown basis kind "%s".', kind);
    end
end

% returns basis: local Nij, curlBasis: curlNij
% basis stores each of the three local edge basis fucntion value at 
% each quadrature point, both in x and y.
% curlBasis stores the curl of each of the three local edge basis functions 
% at each quadrature point.
function [basis, curlBasis] = vector_basis(coordinates, referencePoints)
    validate_triangle_inputs(coordinates, referencePoints);
    x = coordinates(:, 1);
    y = coordinates(:, 2);
    doubleArea = (x(2) - x(1)) * (y(3) - y(1)) - ...
        (x(3) - x(1)) * (y(2) - y(1));
    if doubleArea == 0
        error('fem_basis:DegenerateTriangle', ...
            'Vector basis is undefined on a degenerate triangle.');
    end

    % gradients of the three barycentric coordinates lambda_(1,2,3)
    gradients = [(y(2) - y(3)), (x(3) - x(2)); ...
                (y(3) - y(1)), (x(1) - x(3)); ...
                (y(1) - y(2)), (x(2) - x(1))] / doubleArea;
    referenceBarycentric = [1 - sum(referencePoints, 2), referencePoints];
    physicalBarycentric = referenceBarycentric;
    edgePairs = [1, 2; 2, 3; 3, 1];
    numberOfPoints = size(referencePoints, 1);
    basis = zeros(numberOfPoints, 3, 2);
    curlBasis = zeros(numberOfPoints, 3);

    for edge = 1:3
        first = edgePairs(edge, 1);
        second = edgePairs(edge, 2);
        edgeLength = hypot(x(second) - x(first), y(second) - y(first));
        basis(:, edge, 1) = edgeLength * ...
            (physicalBarycentric(:, first) * gradients(second, 1) - ...
            physicalBarycentric(:, second) * gradients(first, 1));
        basis(:, edge, 2) = edgeLength * ...
            (physicalBarycentric(:, first) * gradients(second, 2) - ...
            physicalBarycentric(:, second) * gradients(first, 2));
        curlBasis(:, edge) = 2 * edgeLength * ...
            (gradients(first, 1) * gradients(second, 2) - ...
            gradients(first, 2) * gradients(second, 1));
    end
end

function [shapeValues, shapeGradients] = nodal_basis(coordinates, referencePoints)
    validate_triangle_inputs(coordinates, referencePoints);
    x = coordinates(:, 1);
    y = coordinates(:, 2);
    doubleArea = (x(2) - x(1)) * (y(3) - y(1)) - ...
        (x(3) - x(1)) * (y(2) - y(1));
    if doubleArea == 0
        error('fem_basis:DegenerateTriangle', ...
            'Nodal basis is undefined on a degenerate triangle.');
    end

    shapeValues = [1 - sum(referencePoints, 2), referencePoints];
    shapeGradients = [(y(2) - y(3)), (x(3) - x(2)); ...
                    (y(3) - y(1)), (x(1) - x(3)); ...
                    (y(1) - y(2)), (x(2) - x(1))] / doubleArea;
end

function [points, weights] = triangle_rule(polynomialDegree)
    if ~isnumeric(polynomialDegree) || ~isscalar(polynomialDegree) || ...
            ~any(polynomialDegree == [1, 2])
        error('fem_basis:UnsupportedQuadrature', ...
            'Supported triangle quadrature degrees are 1 and 2.');
    end

    if polynomialDegree == 1
        points = [1 / 3, 1 / 3];
        weights = 1 / 2;
    else
        points = [1 / 6, 1 / 6; 2 / 3, 1 / 6; 1 / 6, 2 / 3];
        weights = [1 / 6; 1 / 6; 1 / 6];
    end
end

function validate_triangle_inputs(coordinates, referencePoints)
    if ~isnumeric(coordinates) || ~isequal(size(coordinates), [3, 2]) || ...
            any(~isfinite(coordinates(:)))
        error('fem_basis:InvalidCoordinates', ...
            'Triangle coordinates must be a finite 3-by-2 array.');
    end
    if ~isnumeric(referencePoints) || size(referencePoints, 2) ~= 2 || ...
            any(~isfinite(referencePoints(:)))
        error('fem_basis:InvalidReferencePoints', ...
            'Reference points must be a finite N-by-2 array.');
    end
end
