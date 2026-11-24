function plot_utils(action, varargin)
%PLOT_UTILS Dispatch mode-field plots and the combined plot_modes view.
%   plot_utils('plot_modes', mesh, dofs, pec, modes) creates all five plots.
%   Individual plot actions are 'plot_transverse_magnitude',
%   'plot_longitudinal_magnitude', 'plot_longitudinal_real',
%   'plot_transverse_vectors', and 'plot_total_magnitude'.

    if ~(ischar(action) && isrow(action)) && ...
            ~(isstring(action) && isscalar(action))
        error('plot_utils:InvalidAction', ...
            'action must be a character vector or scalar string.');
    end
    action = lower(char(action));
    validActions = {'plot_modes', 'plot_transverse_magnitude', ...
        'plot_longitudinal_magnitude', 'plot_longitudinal_real', ...
        'plot_transverse_vectors', ...
        'plot_total_magnitude'};
    if ~any(strcmp(action, validActions))
        error('plot_utils:UnknownAction', ...
            'Unknown plot action "%s".', action);
    end
    [mesh, dofs, pec, modes] = parse_plot_inputs(varargin);
    fields = reconstruct_mode_fields(mesh, dofs, pec, modes);
    interfaceEdges = find_material_interface_edges(mesh);

    switch action
        case 'plot_modes'
            plot_modes(mesh, fields, modes, interfaceEdges);
        case 'plot_transverse_magnitude'
            plot_transverse_magnitude(mesh, fields, modes, interfaceEdges);
        case 'plot_longitudinal_magnitude'
            plot_longitudinal_magnitude(mesh, fields, modes, interfaceEdges);
        case 'plot_longitudinal_real'
            plot_longitudinal_real(mesh, fields, modes, interfaceEdges);
        case 'plot_transverse_vectors'
            plot_transverse_vectors(mesh, fields, modes, interfaceEdges);
        case 'plot_total_magnitude'
            plot_total_magnitude(mesh, fields, modes, interfaceEdges);
    end
end

function [mesh, dofs, pec, modes] = parse_plot_inputs(inputs)
    if numel(inputs) ~= 4
        error('plot_utils:InvalidInputs', ...
            'Provide mesh, dofs, pec, and modes after the action.');
    end
    mesh = inputs{1};
    dofs = inputs{2};
    pec = inputs{3};
    modes = inputs{4};
    validate_plot_inputs(mesh, dofs, pec, modes);
end

function validate_plot_inputs(mesh, dofs, pec, modes)
    if ~isstruct(mesh) || ~all(isfield(mesh, ...
            {'nodes', 'elements', 'domains', 'edgeNodes', ...
             'elementEdges', 'bcEdgeIds', 'bcEdges', ...
             'elementJacobianInverseTranspose', ...
             'elementJacobianDeterminant'})) || ...
            ~isstruct(dofs) || ~all(isfield(dofs, ...
            {'totalCount', 'vectorRange', 'scalarRange', ...
             'vectorElementIds', 'vectorElementSigns'})) || ...
            ~isstruct(pec) || ~isfield(pec, 'freeGlobalIds') || ...
            ~isstruct(modes) || ~all(isfield(modes, {'vectors', 'n_eff'}))
        error('plot_utils:InvalidInputs', ...
            'Expected mesh, DOF, PEC, and mode data from the FEM workflow.');
    end
    if size(modes.vectors, 1) ~= numel(pec.freeGlobalIds) || ...
            size(modes.vectors, 2) < 1 || ...
            size(modes.vectors, 2) ~= numel(modes.n_eff) || ...
            numel(dofs.vectorRange) ~= size(mesh.edgeNodes, 1) || ...
            numel(dofs.scalarRange) ~= size(mesh.nodes, 1) || ...
            size(dofs.vectorElementIds, 1) ~= size(mesh.elements, 1)
        error('plot_utils:InconsistentInputs', ...
            'The mode vectors and mesh DOF mappings have inconsistent sizes.');
    end
end

function interfaceEdges = find_material_interface_edges(mesh)
    listedEdgeIds = unique(mesh.bcEdgeIds(:));
    interfaceEdgeIds = zeros(0, 1);
    for edgeIndex = 1:numel(listedEdgeIds)
        edgeId = listedEdgeIds(edgeIndex);
        adjacentElements = any(mesh.elementEdges == edgeId, 2);
        adjacentDomains = unique(mesh.domains(adjacentElements));
        if numel(adjacentDomains) > 1
            interfaceEdgeIds(end + 1, 1) = edgeId;
        end
    end
    interfaceEdges = mesh.edgeNodes(interfaceEdgeIds, :);
end

function fields = reconstruct_mode_fields(mesh, dofs, pec, modes)
    numberOfModes = size(modes.vectors, 2);
    numberOfElements = size(mesh.elements, 1);
    numberOfNodes = size(mesh.nodes, 1);
    fields.Ez = zeros(numberOfNodes, numberOfModes);
    fields.longitudinalMagnitude = zeros(numberOfNodes, numberOfModes);
    fields.Et = complex(zeros(numberOfElements, 2, numberOfModes));
    fields.transverseMagnitude = zeros(numberOfElements, numberOfModes);
    fields.totalMagnitude = zeros(numberOfElements, numberOfModes);

    for mode = 1:numberOfModes
        modeVector = zeros(dofs.totalCount, 1);
        modeVector(pec.freeGlobalIds) = modes.vectors(:, mode);
        % Eigenvectors have arbitrary global phase; align it for real-part plots.
        [~, referenceIndex] = max(abs(modeVector));
        if abs(modeVector(referenceIndex)) > 0
            modeVector = modeVector * exp(-1i * angle(modeVector(referenceIndex)));
        end

        ez = modeVector(dofs.scalarRange);
        fields.Ez(:, mode) = real(ez);
        fields.longitudinalMagnitude(:, mode) = abs(ez);
        transverse = complex(zeros(numberOfElements, 2));
        for elementIndex = 1:numberOfElements
            element.coordinates = mesh.nodes(mesh.elements(elementIndex, :), :);
            element.elementJacobianInverseTranspose = ...
                mesh.elementJacobianInverseTranspose(:, :, elementIndex);
            element.elementJacobianDeterminant = ...
                mesh.elementJacobianDeterminant(elementIndex);

            basis = fem_basis('vector', element);
            % Averaging the linear basis over symmetric quadrature gives its centroid value.
            meanBasis = squeeze(mean(basis, 1));
            localDofs = dofs.vectorElementIds(elementIndex, :);
            localSigns = dofs.vectorElementSigns(elementIndex, :).';
            localCoefficients = modeVector(dofs.vectorRange(localDofs));
            localCoefficients = localCoefficients(:) .* localSigns(:);
            transverse(elementIndex, :) = (meanBasis.' * localCoefficients).';
        end

        ezAtElementCenters = mean(ez(mesh.elements), 2);
        fields.Et(:, :, mode) = transverse;
        fields.transverseMagnitude(:, mode) = ...
            sqrt(sum(abs(transverse).^2, 2));
        fields.totalMagnitude(:, mode) = sqrt( ...
            fields.transverseMagnitude(:, mode).^2 + ...
            abs(ezAtElementCenters).^2);
    end
end

function plot_modes(mesh, fields, modes, interfaceEdges)
%PLOT_MODES Create all five mode-field figures.
    plot_transverse_magnitude(mesh, fields, modes, interfaceEdges);
    plot_longitudinal_magnitude(mesh, fields, modes, interfaceEdges);
    plot_longitudinal_real(mesh, fields, modes, interfaceEdges);
    plot_transverse_vectors(mesh, fields, modes, interfaceEdges);
    plot_total_magnitude(mesh, fields, modes, interfaceEdges);
end

function plot_transverse_magnitude(mesh, fields, modes, interfaceEdges)
%PLOT_TRANSVERSE_MAGNITUDE Plot |E_t| for all modes.
    plot_element_scalar_modes(mesh, fields.transverseMagnitude, modes, ...
        '|E_t|', 'Transverse field magnitude', interfaceEdges);
end

function plot_longitudinal_magnitude(mesh, fields, modes, interfaceEdges)
%PLOT_LONGITUDINAL_MAGNITUDE Plot |E_z| for all modes.
    plot_nodal_scalar_modes(mesh, fields.longitudinalMagnitude, modes, ...
        '|E_z|', 'Longitudinal field magnitude', interfaceEdges);
end

function plot_longitudinal_real(mesh, fields, modes, interfaceEdges)
%PLOT_LONGITUDINAL_REAL Plot Re(E_z) for all modes.
    plot_nodal_scalar_modes(mesh, fields.Ez, modes, ...
        'Re(E_z)', 'Real longitudinal field', interfaceEdges);
end

function plot_transverse_vectors(mesh, fields, modes, interfaceEdges)
%PLOT_TRANSVERSE_VECTORS Plot real transverse-field vectors for all modes.
    values = fields.Et;
    numberOfModes = size(values, 3);
    [numberOfRows, numberOfColumns] = tile_dimensions(numberOfModes);
    figure('Name', 'Transverse field vectors', 'NumberTitle', 'off');
    layout = tiledlayout(numberOfRows, numberOfColumns);
    layout.TileSpacing = 'compact';
    layout.Padding = 'compact';
    title(layout, 'Real transverse electric-field vectors');

    elementCenters = (mesh.nodes(mesh.elements(:, 1), :) + ...
        mesh.nodes(mesh.elements(:, 2), :) + ...
        mesh.nodes(mesh.elements(:, 3), :)) / 3;
    for mode = 1:numberOfModes
        nexttile(layout);
        quiver(elementCenters(:, 1), elementCenters(:, 2), ...
            real(values(:, 1, mode)), real(values(:, 2, mode)), ...
            'AutoScale', 'on');
        overlay_material_interfaces(mesh, interfaceEdges);
        axis equal tight;
        title(mode_title(mode, modes));
        xlabel('x');
        ylabel('y');
    end
end

function plot_total_magnitude(mesh, fields, modes, interfaceEdges)
%PLOT_TOTAL_MAGNITUDE Plot |E| for all modes.
    plot_element_scalar_modes(mesh, fields.totalMagnitude, modes, ...
        '|E|', 'Total electric-field magnitude', interfaceEdges);
end

function plot_element_scalar_modes(mesh, values, modes, fieldLabel, ...
        figureLabel, interfaceEdges)
    numberOfModes = size(values, 2);
    [numberOfRows, numberOfColumns] = tile_dimensions(numberOfModes);
    figure('Name', figureLabel, 'NumberTitle', 'off');
    layout = tiledlayout(numberOfRows, numberOfColumns);
    layout.TileSpacing = 'compact';
    layout.Padding = 'compact';
    title(layout, figureLabel);

    vertices = [mesh.nodes, zeros(size(mesh.nodes, 1), 1)];
    for mode = 1:numberOfModes
        nexttile(layout);
        patch('Faces', mesh.elements, 'Vertices', vertices, ...
            'FaceVertexCData', values(:, mode), 'FaceColor', 'flat', ...
            'EdgeColor', 'none');
        overlay_material_interfaces(mesh, interfaceEdges);
        view(2);
        axis equal tight;
        colorbarHandle = colorbar;
        title(mode_title(mode, modes));
        xlabel('x');
        ylabel('y');
        ylabel(colorbarHandle, fieldLabel);
    end
end

function plot_nodal_scalar_modes(mesh, values, modes, fieldLabel, ...
        figureLabel, interfaceEdges)
    numberOfModes = size(values, 2);
    [numberOfRows, numberOfColumns] = tile_dimensions(numberOfModes);
    figure('Name', figureLabel, 'NumberTitle', 'off');
    layout = tiledlayout(numberOfRows, numberOfColumns);
    layout.TileSpacing = 'compact';
    layout.Padding = 'compact';
    title(layout, figureLabel);

    for mode = 1:numberOfModes
        nexttile(layout);
        trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), ...
            zeros(size(mesh.nodes, 1), 1), values(:, mode), ...
            'EdgeColor', 'none');
        overlay_material_interfaces(mesh, interfaceEdges);
        view(2);
        axis equal tight;
        colorbarHandle = colorbar;
        title(mode_title(mode, modes));
        xlabel('x');
        ylabel('y');
        ylabel(colorbarHandle, fieldLabel);
    end
end

function overlay_material_interfaces(mesh, interfaceEdges)
    if isempty(interfaceEdges)
        return;
    end
    x = [mesh.nodes(interfaceEdges(:, 1), 1), ...
         mesh.nodes(interfaceEdges(:, 2), 1), ...
         nan(size(interfaceEdges, 1), 1)].';
    y = [mesh.nodes(interfaceEdges(:, 1), 2), ...
         mesh.nodes(interfaceEdges(:, 2), 2), ...
         nan(size(interfaceEdges, 1), 1)].';
    hold on;
    plot(x(:), y(:), 'k-', 'LineWidth', 1.2);
    hold off;
end

function [numberOfRows, numberOfColumns] = tile_dimensions(numberOfModes)
    numberOfColumns = ceil(sqrt(numberOfModes));
    numberOfRows = ceil(numberOfModes / numberOfColumns);
end

function label = mode_title(mode, modes)
    label = sprintf('%d. n_{eff} = %s', mode, ...
        num2str(modes.n_eff(mode), '%.5g'));
end
