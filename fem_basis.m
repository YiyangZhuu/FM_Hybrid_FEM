% FEM_BASIS
% First-order triangle bases and reference-cell/edge quadrature rules.
%   [N, curlN, basisNormal] = fem_basis('vector', element, rule)
%   [N, gradN] = fem_basis('nodal', element, rule)
%   [trianglePoints, triangleWeights, edgePoints, edgeWeights, edgeIds] = ...
%       fem_basis('quadrature', polynomialDegree)
%
% Input:
% element: triangle element data containing coordinates and precomputed
%          elementJacobianInverseTranspose
% rule: optional evaluation label 'triangle', 'edge', or 'line'; defaults to 'triangle'
% 'edge' returns physical vector values at edge quadrature points.
% 'line' returns physical vector values in N and dot(N, edge normal) in basisNormal.
% basisNormal is empty for 'triangle' and 'edge'.
% edgePoints: reference coordinates of the 2-point Gauss rule on each edge
% edgeWeights: parameter weights (sum 1 on every edge); multiply by physical
%              edge length to obtain physical line-integration weights
% edgeIds: local directed edge id (1=1->2, 2=2->3, 3=3->1)
% Output:
% N [numReferencePoints x 3]: shape function values evaluated at reference points with reference shape functions
% gradN [2 x 3]: physical gradient of each nodal shape function, with each
%               node stored as a column vector component (x/y by node).
% Vector functions are ordered on local directed edges (1->2, 2->3, 3->1).
% Apply mesh.elementEdgeSigns when mapping local edge DOFs to global DOFs.
% Quadrature weights integrate on the reference triangle (area 1/2).
% Reference nodal values, reference edge bases, and reference curls are cached.

% TODO:eliminate or the if case that catch that elementData does not contain that field.

% Each element transforms them once with its Jacobian inverse transpose.

function varargout = fem_basis(kind, varargin)
    switch lower(kind)
        case 'vector'
            [varargout{1:nargout}] = vector_basis(varargin{:});
        case 'nodal'
            [varargout{1:2}] = nodal_basis(varargin{:});
        case 'quadrature'
            [varargout{1:nargout}] = quadrature_rule(varargin{:});
        otherwise
            error('fem_basis:UnknownKind', 'Unknown basis kind "%s".', kind);
    end
end

% returns basis: local Nij, curlBasis: curlNij
% basis stores each of the three local edge basis fucntion value at 
% each quadrature point, both in x and y.
% curlBasis stores the curl of each of the three local edge basis functions 
% at each quadrature point.
function [basis, curlBasis, basisNormal] = vector_basis(element, rule)
    if nargin < 2
        rule = 'triangle';
    end
    referencePoints = basis_reference_points(rule);
    validate_basis_inputs(element, rule);

    [referenceBasis, referenceCurl] = reference_vector_basis(referencePoints);
    jacobianInverseTranspose = element.elementJacobianInverseTranspose;
    physicalBasis = zeros(size(referenceBasis));
    for point = 1:size(referenceBasis, 1)
        for edge = 1:size(referenceBasis, 2)
            referenceVector = squeeze(referenceBasis(point, edge, :));
            physicalBasis(point, edge, :) = jacobianInverseTranspose * referenceVector;
        end
    end

    basisNormal = [];
    if strcmp(lower(char(rule)), 'line')
        edgeNormals = element_edge_normals(element);
        [~, ~, ~, ~, edgeIds] = quadrature_rule(2);
        basisNormal = zeros(size(referenceBasis, 1), size(referenceBasis, 2));
        for point = 1:size(referenceBasis, 1)
            edge = edgeIds(point);
            normal = edgeNormals(edge, :);
            basisNormal(point, :) = squeeze(physicalBasis(point, :, :)) * normal(:);
        end
    end
    basis = physicalBasis;

    if isfield(element, 'elementJacobianDeterminant')
        jacobianDeterminant = element.elementJacobianDeterminant;
    elseif isfield(element, 'elementArea')
        jacobianDeterminant = 2 * element.elementArea;
    else
        jacobianDeterminant = abs(det(inv(jacobianInverseTranspose)));
    end
    curlBasis = referenceCurl / jacobianDeterminant;
end

function [shapeValues, shapeGradients] = nodal_basis(element, rule)
    if nargin < 2
        rule = 'triangle';
    end
    referencePoints = basis_reference_points(rule);
    validate_basis_inputs(element, rule);

    % Depending on what rules we are using
    shapeValues = reference_nodal_basis(referencePoints);
    % Gradient is stored as a column vector
    referenceGradients = [-1, 1, 0; 
                          -1, 0, 1]; % Fixed value, independent of reference values
    jacobianInverseTranspose = element.elementJacobianInverseTranspose;

    shapeGradients = jacobianInverseTranspose * referenceGradients;
end

function referencePoints = basis_reference_points(rule)
    if ~ischar(rule) && ~isstring(rule)
        error('fem_basis:InvalidRule', ...
            'Basis evaluation rule must be ''triangle'', ''edge'', or ''line''.');
    end
    rule = lower(char(rule));
    if ~any(strcmp(rule, {'triangle', 'edge', 'line'}))
        error('fem_basis:InvalidRule', ...
            'Basis evaluation rule must be ''triangle'', ''edge'', or ''line''.');
    end
    [trianglePoints, ~, edgePoints] = quadrature_rule(2);
    if strcmp(rule, 'triangle')
        referencePoints = trianglePoints;
    else
        referencePoints = edgePoints;
    end
end

function shapeValues = reference_nodal_basis(referencePoints)
    persistent cache
    if isempty(cache) || ~isequal(cache.referencePoints, referencePoints)
        cache.referencePoints = referencePoints;
        cache.shapeValues = [1 - sum(referencePoints, 2), referencePoints];
    end
    shapeValues = cache.shapeValues;
end

function [basis, curlBasis] = reference_vector_basis(referencePoints)
    persistent cache
    if isempty(cache) || ~isequal(cache.referencePoints, referencePoints)
        referenceGradients = [-1, -1; 1, 0; 0, 1];
        referenceBarycentric = [1 - sum(referencePoints, 2), referencePoints];
        edgePairs = [1, 2; 2, 3; 3, 1];
        numberOfPoints = size(referencePoints, 1);
        referenceBasis = zeros(numberOfPoints, 3, 2);
        referenceCurl = zeros(numberOfPoints, 3);
        for edge = 1:3
            first = edgePairs(edge, 1);
            second = edgePairs(edge, 2);
            referenceBasis(:, edge, 1) = referenceBarycentric(:, first) * ...
                referenceGradients(second, 1) - referenceBarycentric(:, second) * ...
                referenceGradients(first, 1);
            referenceBasis(:, edge, 2) = referenceBarycentric(:, first) * ...
                referenceGradients(second, 2) - referenceBarycentric(:, second) * ...
                referenceGradients(first, 2);
            referenceCurl(:, edge) = 2 * ( ...
                referenceGradients(first, 1) * referenceGradients(second, 2) - ...
                referenceGradients(first, 2) * referenceGradients(second, 1));
        end
        cache.referencePoints = referencePoints;
        cache.basis = referenceBasis;
        cache.curl = referenceCurl;
    end
    basis = cache.basis;
    curlBasis = cache.curl;
end

function varargout = quadrature_rule(polynomialDegree)
    if nargin < 1 || isempty(polynomialDegree)
        polynomialDegree = 2;
    end
    if ~isnumeric(polynomialDegree) || ~isscalar(polynomialDegree) || ...
            ~any(polynomialDegree == [1, 2])
        error('fem_basis:UnsupportedQuadrature', ...
            'Supported triangle quadrature degrees are 1 and 2.');
    end

    if polynomialDegree == 1
        trianglePoints = [1 / 3, 1 / 3];
        triangleWeights = 1 / 2;
    else
        trianglePoints = [1 / 6, 1 / 6; 2 / 3, 1 / 6; 1 / 6, 2 / 3];
        triangleWeights = [1 / 6; 1 / 6; 1 / 6];
    end

    % Two-point Gauss rule on t in [0,1] for each directed local edge.
    % It integrates polynomials through degree 3 along an edge.
    gaussParameters = (1 / 2) + (1 / (2 * sqrt(3))) * [-1; 1];
    edgePoints = [...
        1 - gaussParameters, gaussParameters * 0; ...
        gaussParameters, 1 - gaussParameters; ...
        gaussParameters * 0, 1 - gaussParameters];
    edgeWeights = 1 / 2 * ones(6, 1);
    edgeIds = [1; 1; 2; 2; 3; 3];

    varargout = {trianglePoints, triangleWeights, edgePoints, edgeWeights, edgeIds};
end

function validate_basis_inputs(element, rule)
    if ~isstruct(element) || ~isfield(element, 'coordinates')
        error('fem_basis:InvalidElement', ...
            'Element must be a struct containing coordinates.');
    end
    coordinates = element.coordinates;
    if ~isnumeric(coordinates) || ~isequal(size(coordinates), [3, 2]) || ...
            any(~isfinite(coordinates(:)))
        error('fem_basis:InvalidCoordinates', ...
            'Triangle coordinates must be a finite 3-by-2 array.');
    end
    if isfield(element, 'elementJacobianInverseTranspose')
        jacobianInverseTranspose = element.elementJacobianInverseTranspose;
        if ~isnumeric(jacobianInverseTranspose) || ...
                ~isequal(size(jacobianInverseTranspose), [2, 2]) || ...
                any(~isfinite(jacobianInverseTranspose(:)))
            error('fem_basis:InvalidJacobianInverseTranspose', ...
                'Element elementJacobianInverseTranspose must be a finite 2-by-2 matrix.');
        end
    end
    if isfield(element, 'elementJacobianDeterminant') && ...
            (~isnumeric(element.elementJacobianDeterminant) || ...
             ~isscalar(element.elementJacobianDeterminant) || ...
             ~isfinite(element.elementJacobianDeterminant) || ...
             element.elementJacobianDeterminant <= 0)
        error('fem_basis:InvalidJacobianDeterminant', ...
            'Element elementJacobianDeterminant must be finite and positive.');
    end
    if ~ischar(rule) && ~isstring(rule)
        error('fem_basis:InvalidRule', ...
            'Basis evaluation rule must be ''triangle'', ''edge'', or ''line''.');
    end
    rule = lower(char(rule));
    if ~any(strcmp(rule, {'triangle', 'edge', 'line'}))
        error('fem_basis:InvalidRule', ...
            'Basis evaluation rule must be ''triangle'', ''edge'', or ''line''.');
    end
    if strcmp(rule, 'line')
        edgeNormals = element_edge_normals(element);
        if ~isnumeric(edgeNormals) || ~isequal(size(edgeNormals), [3, 2]) || ...
                any(~isfinite(edgeNormals(:)))
            error('fem_basis:InvalidElementEdgeNormals', ...
                'Element elementEdgeNormals must be a finite 3-by-2 array.');
        end
    end
end

function edgeNormals = element_edge_normals(element)
    if isfield(element, 'elementEdgeNormals')
        edgeNormals = element.elementEdgeNormals;
    elseif isfield(element, 'elementNormals')
        edgeNormals = element.elementNormals;
    elseif isfield(element, 'normal')
        edgeNormals = element.normal;
    else
        error('fem_basis:MissingElementEdgeNormals', ...
            'Line-rule vector basis requires elementEdgeNormals or elementNormals.');
    end
end
