function results = compare_transonic_solvers(p, Minf, solvers, Grid_size, outDir, resDir)
% Compare algebraic solvers on a user-selected grid.
% Records nonlinear residual histories, computational time and convergence.
% Plots the figures residual vs iterations and residual vs computational
% time in linear and log scales.
% It prints the iterations, final residual, max(mach), computational time,
% and if it converged.

    % Study grid is supplied by the run script
    Nx = Grid_size(1,1);
    Ny = Grid_size(1,2);

    if ~(Minf >= 0.8 && Minf <= 1.2)
        error('Solver comparison requires a transonic Mach number.');
    end

    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(resDir, 'dir'), mkdir(resDir); end

    nSolver = numel(solvers);
    solutions = cell(nSolver, 1);

    solverName = strings(nSolver, 1);
    Time_s = zeros(nSolver, 1);
    Iterations = zeros(nSolver, 1);
    FinalResidual = zeros(nSolver, 1);
    Converged = false(nSolver, 1);
    MaxMach = zeros(nSolver, 1);

    % ==============================================================
    % Run each solver with identical numerical parameters
    % ==============================================================
    for j = 1:nSolver

        trial = p;
        trial.Nx = Nx;
        trial.Ny = Ny;
        trial.solver = solvers{j};

        trial.recordSolverHistory = true;
        trial.makePlots = false;
        trial.silentSolver = true;

        fprintf('Solver comparison: %s, grid %d x %d\n', ...
            trial.solver, Nx, Ny);

        s = cfm.solve_flow_fv(trial, Minf);

        % Store results independently
        solutions{j} = s;

        solverName(j) = string(trial.solver);
        Time_s(j) = s.computationalTime;
        Iterations(j) = s.iterations;
        FinalResidual(j) = s.residual;
        Converged(j) = s.converged;
        MaxMach(j) = max(s.Mlocal(:));

        if ~s.converged
            warning('%s did not converge (residual %.3e).', ...
                trial.solver, s.residual);
        end
    end

    % ==============================================================
    % Results table
    % ==============================================================
    results = table(solverName, Time_s, Iterations, ...
        FinalResidual, Converged, MaxMach, ...
        'VariableNames', {'Solver','Time_s', ...
        'NonlinearIterations','FinalResidual', ...
        'Converged','MaxMach'});

    writetable(results, ...
        fullfile(resDir, 's2_solver_comparison.csv'));

    disp(results);

    colors = lines(nSolver);
    lineStyles = {'-', '--', '-.', ':', '-'};


    % ==============================================================
    % FIGURE 1A: Residual vs Nonlinear Iteration - LINEAR SCALE
    % ==============================================================
    fig = figure('Color', 'w', 'Visible', 'on', ...
        'Position', [100 100 850 550]);

    hold on;

    for j = 1:nSolver
        s = solutions{j};

        plot(1:numel(s.residualHistory), ...
            s.residualHistory(:), ...
            'LineStyle', lineStyles{j}, ...
            'Color', colors(j,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', solvers{j});
    end

    hold off;
    set(gca, 'YScale', 'linear','FontSize', 15);
    grid on;

    xlabel('Iterations', 'FontSize', 19);
    ylabel('Residual', 'FontSize', 19);

    title(sprintf( ...
        'Residual vs Iteration (Linear Scale), M_\\infty = %.2f, %d x %d', ...
        Minf, Nx, Ny),'FontSize',16);

    legend('Location', 'best','FontSize',15);

    exportgraphics(fig, fullfile(outDir, ...
        's2_solver_residual_vs_iteration_linear.pdf'));


    % ==============================================================
    % FIGURE 1B: Residual vs Nonlinear Iteration - LOG SCALE
    % ==============================================================
    fig = figure('Color', 'w', 'Visible', 'on', ...
        'Position', [150 100 850 550]);

    hold on;

    for j = 1:nSolver
        s = solutions{j};

        semilogy(1:numel(s.residualHistory), ...
            max(s.residualHistory(:), realmin), ...
            'LineStyle', lineStyles{j}, ...
            'Color', colors(j,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', solvers{j});
    end

    hold off;
    set(gca, 'YScale', 'log','FontSize', 15);
    grid on;

    xlabel('Iterations', 'FontSize', 19);
    ylabel('Residual', 'FontSize', 19);

    title(sprintf( ...
        'Residual vs Iteration (Log Scale), M_\\infty = %.2f, %d x %d', ...
        Minf, Nx, Ny),'FontSize',16);

    legend('Location', 'best','FontSize',15);

    exportgraphics(fig, fullfile(outDir, ...
        's2_solver_residual_vs_iteration_log.pdf'));


    % ==============================================================
    % FIGURE 2A: Residual vs Computational Time - LINEAR SCALE
    % ==============================================================
    fig = figure('Color', 'w', 'Visible', 'on', ...
        'Position', [200 100 850 550]);

    hold on;

    for j = 1:nSolver
        s = solutions{j};

        plot(s.elapsedHistory(:), ...
            s.residualHistory(:), ...
            'LineStyle', lineStyles{j}, ...
            'Color', colors(j,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', solvers{j});
    end

    hold off;
    set(gca, 'YScale', 'linear','FontSize', 15);
    grid on;

    xlabel('Elapsed computational time (s)', 'FontSize', 19);
    ylabel('Residual', 'FontSize', 19);

    title(sprintf( ...
        'Residual vs Time (Linear Scale), M_\\infty = %.2f, %d x %d', ...
        Minf, Nx, Ny),'FontSize',16);

    legend('Location', 'best','FontSize',15);

    exportgraphics(fig, fullfile(outDir, ...
        's2_solver_residual_vs_time_linear.pdf'));


    % ==============================================================
    % FIGURE 2B: Residual vs Computational Time - LOG SCALE
    % ==============================================================
    fig = figure('Color', 'w', 'Visible', 'on', ...
        'Position', [250 100 850 550]);

    hold on;

    for j = 1:nSolver
        s = solutions{j};

        semilogy(s.elapsedHistory(:), ...
            max(s.residualHistory(:), realmin), ...
            'LineStyle', lineStyles{j}, ...
            'Color', colors(j,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', solvers{j});
    end

    hold off;
    set(gca, 'YScale', 'log','FontSize', 15);
    grid on;

    xlabel('Elapsed computational time (s)', 'FontSize', 19);
    ylabel('Residual', 'FontSize', 19);

    title(sprintf( ...
        'Residual vs Time (Log Scale), M_\\infty = %.2f, %d x %d', ...
        Minf, Nx, Ny),'FontSize',16);

    legend('Location', 'best','FontSize',15);

    exportgraphics(fig, fullfile(outDir, ...
        's2_solver_residual_vs_time_log.pdf'));

end
