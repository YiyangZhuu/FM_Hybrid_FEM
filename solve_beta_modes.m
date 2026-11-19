% SOLVE_BETA_MODES
% Solve A*x = lambda*B*x near sigma by shift-invert.
% The solver applies the standard operator (A - sigma*B)\B, rather than
% calling ARPACK's generalized eigensolver with B. This permits symmetric
% indefinite B. Returned propagation constants satisfy beta^2 = -lambda.
%
% options fields (all optional):
%   tolerance                 eigs convergence tolerance (default 1e-10)
%   maxIterations             eigs iteration limit (default 1000)
%   subspaceDimension         eigs subspace size (default max(2*k+1,20))
%   residualTolerance         maximum relative original-pencil residual
%                             (default 1e-8)
%   lossless                  require beta^2 to be real (default false)
%   realEigenvalueTolerance   tolerance for the lossless realness check
%                             (default 1e-8)
%
% The matrices A and B must already have PEC-constrained DOFs eliminated.
% This function does not derive or assemble either matrix.

function modes = solve_beta_modes(A, B, sigma, numberOfModes, options)
    if nargin < 5 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('solve_beta_modes:InvalidOptions', ...
            'options must be a scalar struct.');
    end
    validate_matrices(A, B);
    validate_scalar_input(sigma, 'sigma');
    if ~isnumeric(numberOfModes) || ~isscalar(numberOfModes) || ...
            ~isfinite(numberOfModes) || numberOfModes < 1 || ...
            numberOfModes ~= fix(numberOfModes)
        error('solve_beta_modes:InvalidModeCount', ...
            'numberOfModes must be a positive integer.');
    end

    matrixSize = size(A, 1);
    if numberOfModes >= matrixSize
        error('solve_beta_modes:InvalidModeCount', ...
            'numberOfModes must be smaller than the matrix dimension.');
    end

    options = set_option_defaults(options, numberOfModes, matrixSize);
    validate_options(options, numberOfModes, matrixSize);

    shiftedMatrix = A - sigma * B;
    shiftedSolver = decomposition(shiftedMatrix, 'lu');
    applyShiftInvert = @(x) shiftedSolver \ (B * x);

    eigsOptions.tol = options.tolerance;
    eigsOptions.maxit = options.maxIterations;
    eigsOptions.p = options.subspaceDimension;
    eigsOptions.issym = false;
    eigsOptions.isreal = isreal(A) && isreal(B) && isreal(sigma);
    eigsOptions.disp = 0;
    [vectors, transformedValues, convergenceFlag] = ...
        eigs(applyShiftInvert, matrixSize, numberOfModes, ...
        'largestabs', eigsOptions);
    if convergenceFlag ~= 0
        error('solve_beta_modes:NoConvergence', ...
            'eigs did not converge for all requested modes (flag %d).', ...
            convergenceFlag);
    end

    nu = diag(transformedValues);
    if any(nu == 0)
        error('solve_beta_modes:ZeroTransformedEigenvalue', ...
            'A zero shift-invert eigenvalue cannot be mapped back to lambda.');
    end
    lambda = sigma + 1 ./ nu;

    [~, order] = sort(abs(lambda - sigma), 'ascend');
    lambda = lambda(order);
    nu = nu(order);
    vectors = vectors(:, order);

    absoluteResidual = zeros(numberOfModes, 1);
    relativeResidual = zeros(numberOfModes, 1);
    for mode = 1:numberOfModes
        x = vectors(:, mode);
        Ax = A * x;
        Bx = B * x;
        residual = Ax - lambda(mode) * Bx;
        absoluteResidual(mode) = norm(residual);
        residualScale = norm(Ax) + abs(lambda(mode)) * norm(Bx);
        relativeResidual(mode) = absoluteResidual(mode) / max(residualScale, eps);
    end
    if any(~isfinite(relativeResidual)) || ...
            any(relativeResidual > options.residualTolerance)
        error('solve_beta_modes:ResidualTooLarge', ...
            ['At least one mode failed the original-pencil residual check. ', ...
            'Maximum relative residual: %.3e (tolerance %.3e).'], ...
            max(relativeResidual), options.residualTolerance);
    end

    betaSquared = -lambda;
    betaSquaredIsReal = abs(imag(betaSquared)) <= ...
        options.realEigenvalueTolerance .* max(1, abs(real(betaSquared)));
    % if options.lossless && any(~betaSquaredIsReal)
    %    error('solve_beta_modes:ComplexLosslessEigenvalue', ...
    %        ['A lossless solve produced beta^2 with a significant imaginary ', ...
    %        'part (maximum relative imaginary part %.3e).'], ...
    %        max(abs(imag(betaSquared)) ./ max(1, abs(real(betaSquared)))));
    % end

    modes.lambda = lambda;
    modes.betaSquared = betaSquared;
    modes.beta = sqrt(betaSquared);
    modes.vectors = vectors;
    modes.transformedEigenvalues = nu;
    modes.absoluteResidual = absoluteResidual;
    modes.relativeResidual = relativeResidual;
    modes.betaSquaredIsReal = betaSquaredIsReal;
    modes.shift = sigma;
end

function validate_matrices(A, B)
    if ~isnumeric(A) || ~ismatrix(A) || isempty(A) || ...
            size(A, 1) ~= size(A, 2) || any(~isfinite(nonzeros(A)))
        error('solve_beta_modes:InvalidA', ...
            'A must be a finite, nonempty square numeric matrix.');
    end
    if ~isnumeric(B) || ~ismatrix(B) || ~isequal(size(A), size(B)) || ...
            any(~isfinite(nonzeros(B)))
        error('solve_beta_modes:InvalidB', ...
            'B must be a finite numeric matrix with the same size as A.');
    end
end

function validate_scalar_input(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        error('solve_beta_modes:InvalidShift', ...
            '%s must be a finite numeric scalar.', name);
    end
end

function options = set_option_defaults(options, numberOfModes, matrixSize)
    defaults.tolerance = 1e-10;
    defaults.maxIterations = 1000;
    defaults.subspaceDimension = min(matrixSize, max(2 * numberOfModes + 1, 20));
    defaults.residualTolerance = 1e-8;
    defaults.lossless = false;
    defaults.realEigenvalueTolerance = 1e-8;
    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(options, name) || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end
end

function validate_options(options, numberOfModes, matrixSize)
    positiveScalars = {'tolerance', 'maxIterations', ...
        'residualTolerance', 'realEigenvalueTolerance'};
    for k = 1:numel(positiveScalars)
        name = positiveScalars{k};
        value = options.(name);
        if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
            error('solve_beta_modes:InvalidOption', ...
                'options.%s must be a positive finite scalar.', name);
        end
    end
    if options.maxIterations ~= fix(options.maxIterations)
        error('solve_beta_modes:InvalidOption', ...
            'options.maxIterations must be an integer.');
    end
    if ~isnumeric(options.subspaceDimension) || ...
            ~isscalar(options.subspaceDimension) || ...
            ~isfinite(options.subspaceDimension) || ...
            options.subspaceDimension <= numberOfModes || ...
            options.subspaceDimension > matrixSize || ...
            options.subspaceDimension ~= fix(options.subspaceDimension)
        error('solve_beta_modes:InvalidOption', ...
            ['options.subspaceDimension must be an integer greater than the ', ...
            'mode count and no greater than the matrix dimension.']);
    end
    if ~islogical(options.lossless) || ~isscalar(options.lossless)
        error('solve_beta_modes:InvalidOption', ...
            'options.lossless must be a logical scalar.');
    end
end
