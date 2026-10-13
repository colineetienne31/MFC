function out = grid_study(p, g, outDir)
    % grid_study L2 error (exact supersonic solution), number of cells and CPU time when refining in x, in y or in both.
    nf = numel(g.lists);
    out = cell(1, nf);
    for f = 1:nf
        L = g.lists{f};
        n = size(L, 1);
        res.grid = L;  res.N = L(:,1) .* L(:,2);
        res.err = zeros(n, 1);  res.time = zeros(n, 1);
        for k = 1:n
            q = p;  q.Nx = L(k,1);  q.Ny = L(k,2);
            t = inf;
            for rep = 1:g.reps                      % minimum over repetitions reduces timing noise
                tic; s = cfm.solve_flow_fv(q, g.Minf); t = min(t, toc);
            end
            e = cfm.error_l2(s.phiNorm - cfm.exact_supersonic(s.X, s.Y, q, g.Minf), 2);
            res.err(k) = e(1);  res.time(k) = t;
            fprintf('grid study (%s): %d x %d done\n', g.names{f}, q.Nx, q.Ny);
        end
        out{f} = res;
    end

    % LaTeX rows: family, grid, N, error, order (w.r.t. the refined mesh size), CPU time
    fprintf('\n%% --- Grid study rows ---\n');
    for f = 1:nf
        r = out{f};
        for k = 1:size(r.grid, 1)
            if k == 1, ps = '--'; else, ps = sprintf('%.2f', log2(r.err(k-1) / r.err(k))); end
            fprintf('%s & $%d\\times%d$ & %d & %.2e & %s & %.2f \\\\\n', g.names{f}, ...
                r.grid(k,1), r.grid(k,2), r.N(k), r.err(k), ps, r.time(k));
        end
    end

    % Figure: error against CPU time and against number of cells
    mk = {'o-', 's-', '^-'};
    figure('Color', 'w', 'Position', [100 100 900 360]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on;
    for f = 1:nf, plot(out{f}.time, out{f}.err, mk{f}); end
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
    xlabel('CPU time (s)'); ylabel('L_2 error of \phi/(U_\infty c)'); legend(g.names, 'Location', 'northeast');
    nexttile; hold on;
    for f = 1:nf, plot(out{f}.N, out{f}.err, mk{f}); end
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
    xlabel('number of cells N_x N_y'); ylabel('L_2 error of \phi/(U_\infty c)');
    exportgraphics(gcf, fullfile(outDir, 's2_grid_study_cost.pdf'));
end
