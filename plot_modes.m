function plot_modes(mesh, dofs, pec, modes)
%PLOT_FIRST_MODE Plot the longitudinal field magnitude of all modes.

    numberOfModes = size(modes.vectors, 2);
    figure;
    title('|E_z|');
    tiledlayout(2, 3);

    for mode = 1:numberOfModes
        modeVector = zeros(dofs.totalCount, 1);
        modeVector(pec.freeGlobalIds) = modes.vectors(:, mode);
        ez = modeVector(dofs.scalarRange);

        nexttile;
        trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), ...
            zeros(size(mesh.nodes, 1), 1), abs(ez), ...
            'EdgeColor', 'none');
        view(2);
        axis equal tight;
        colorbar;
        title(sprintf('%d:, n_{eff} = %s', mode, ...
            num2str(modes.n_eff(mode), '%.4g')));
        xlabel('x');
        ylabel('y');
    end
end
