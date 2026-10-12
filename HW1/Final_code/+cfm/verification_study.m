function out = verification_study(p, v, outDir)
    % verification_study L2 errors and observed orders against h, Mach-wave angle, figure and LaTeX table rows.
    n = size(v.gridList, 1);
    h = p.xMax ./ v.gridList(:,1);
    eSup = zeros(n, 3);  phiSub = cell(n, 1);  betaNum = zeros(n, 1);
    for k = 1:n
        q = p;  q.Nx = v.gridList(k,1);  q.Ny = v.gridList(k,2);
        s = cfm.solve_flow_fv(q, v.Msup);
        eSup(k,:)  = cfm.error_l2(s.phiNorm - cfm.exact_supersonic(s.X, s.Y, q, v.Msup), v.nb);
        betaNum(k) = cfm.wave_angle(s.X, s.Y, s.phiNorm, v.yWave);
        s = cfm.solve_flow_fv(q, v.Msub);
        phiSub{k} = s.phiNorm;
        fprintf('verification: grid %d x %d done\n', q.Nx, q.Ny);
    end

    % Subsonic: difference between a grid and the next finer one (restricted 2x2 blocks)
    eSub = zeros(n-1, 3);
    for k = 1:n-1
        eSub(k,:) = cfm.error_l2(phiSub{k} - cfm.restrict2(phiSub{k+1}), v.nb);
    end

    % Observed orders (columns: all cells, interior, boundary band)
    pSup = log(eSup(2:end,:) ./ eSup(1:end-1,:)) ./ log(h(2:end) ./ h(1:end-1));
    pSub = log2(eSub(1:end-1,:) ./ eSub(2:end,:));       % three grids, mesh-size ratio 2
    beta = asind(1 / v.Msup);                            % theoretical Mach-wave angle
    out = struct('h', h, 'eSup', eSup, 'eSub', eSub, 'pSup', pSup, 'pSub', pSub, ...
                 'betaNum', betaNum, 'beta', beta);

    % Mach-wave angle
    fprintf('\nMach-wave angle, theory beta = %.2f deg\n', beta);
    for k = 1:n
        fprintf('  %d x %d : measured %.2f deg (difference %.2f deg)\n', ...
            v.gridList(k,1), v.gridList(k,2), betaNum(k), betaNum(k) - beta);
    end

    % Figure
    figure('Color', 'w', 'Position', [100 100 900 360]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    loglog(h, eSup(:,1), 'o-', h, eSup(:,2), 's-', h, eSup(:,3), '^-'); hold on;
    loglog(h, eSup(1,1)*(h/h(1)), 'k--', h, eSup(1,1)*(h/h(1)).^2, 'k:');  grid on;
    xlabel('h (m)'); ylabel('L_2 error of \phi/(U_\infty c)');
    title(sprintf('Supersonic M_\\infty=%.1f (exact solution)', v.Msup));
    legend('all cells', 'interior', 'boundary band', 'slope 1', 'slope 2', 'Location', 'southeast');
    nexttile;
    hc = h(1:n-1);
    loglog(hc, eSub(:,1), 'o-', hc, eSub(:,2), 's-', hc, eSub(:,3), '^-'); hold on;
    loglog(hc, eSub(1,1)*(hc/hc(1)), 'k--', hc, eSub(1,1)*(hc/hc(1)).^2, 'k:');  grid on;
    xlabel('h (m)'); ylabel('L_2 difference to next finer grid');
    title(sprintf('Subsonic M_\\infty=%.1f (successive grids)', v.Msub));
    exportgraphics(gcf, fullfile(outDir, 's2_verification_error.pdf'));

    % LaTeX table rows: grid, h, p (all), p (interior), p (boundary)
    fprintf('\n%% --- Supersonic rows ---\n');
    for k = 1:n
        if k == 1, s1 = '-- & -- & --'; else, s1 = sprintf('%.2f & %.2f & %.2f', pSup(k-1,:)); end
        fprintf('$%d\\times%d$ & %.4f & %s \\\\\n', v.gridList(k,1), v.gridList(k,2), h(k), s1);
    end
    fprintf('%% --- Subsonic rows ---\n');
    for k = 1:n
        if k < 3, s1 = '-- & -- & --'; else, s1 = sprintf('%.2f & %.2f & %.2f', pSub(k-2,:)); end
        fprintf('$%d\\times%d$ & %.4f & %s \\\\\n', v.gridList(k,1), v.gridList(k,2), h(k), s1);
    end
end
