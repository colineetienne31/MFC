%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Fluid Mechanics - Homework 1
% Section 1.2: Steady-State Convection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Output directories setup
out_dir = fullfile(pwd, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% =========================================================================
% 1. Parameter Block
% =========================================================================
p.L       = 40;                 % Domain length
p.a_val   = 1.0;                % Convection velocity
p.u_inlet = 1.0;                % Inlet Boundary Condition at x=0
p.N_ref   = 40;                 % Baseline grid size for single run
p.FULL    = true;               % Flag to toggle mesh convergence study

% Exact analytical solution and source term definitions
p.u_exact_fn = @(x) exp(sin(0.2 * x));
p.S_fn       = @(x) p.a_val * (0.2 * (cos(0.2 * x)) .* exp(sin(0.2 * x)));

% =========================================================================
% 2. Scheme Selection
% =========================================================================
scheme_test = 'backward'; % Options: 'backward', 'central', 'all'

if strcmpi(scheme_test, 'backward')
    scheme_list = {'backward'};
elseif strcmpi(scheme_test, 'central')
    scheme_list = {'central'};
elseif strcmpi(scheme_test, 'all')
    scheme_list = {'backward', 'central'};
end

% -------------------------------------------------------------------------
% Manual overrides for parametric testing
% (Uncomment any line below to override defaults)
% -------------------------------------------------------------------------
% p.N_ref = 80;         % Run baseline test with a finer grid
% p.FULL  = false;      % Skip mesh convergence study to save time
% p.a_val = -1.0;       % Reverse flow direction (requires outlet condition handling!)

% =========================================================================
% 3. Execution of Test Cases
% =========================================================================
for m = 1:length(scheme_list)
    current_scheme = scheme_list{m};
    fprintf('\n==========================================================\n');
    fprintf('Running Convection with Scheme: %s\n', upper(current_scheme));
    fprintf('==========================================================\n');
    
    % --- Baseline Grid Solution ---
    [x_ref, dx_ref, u_num_ref] = cfm.solve_advection_1d( ...
        p.L, p.N_ref, p.a_val, p.u_inlet, p.S_fn, current_scheme);
    
    u_exact_ref = p.u_exact_fn(x_ref);
    local_error_ref = abs(u_num_ref - u_exact_ref);

    % --- Figure 1: Numerical vs Exact Solution ---
    fig1 = figure('Color', 'w', 'Name', ['Solution: ' current_scheme]);
    plot(x_ref, u_exact_ref, 'g-', 'LineWidth', 1.8, 'DisplayName', 'Analytical'); hold on;
    plot(x_ref, u_num_ref, 'ro--', 'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', 'Numerical');
    grid on; xlabel('Axial Position x'); ylabel('u(x)'); 
    title(sprintf('Numerical vs Analytical (N = %d, %s)', p.N_ref, current_scheme));
    legend('Location', 'best');
    exportgraphics(fig1, fullfile(out_dir, sprintf('s1_conv_sol_%s.pdf', current_scheme)));

    % --- Figure 2: Local Absolute Error Distribution ---
    fig2 = figure('Color', 'w', 'Name', ['Error: ' current_scheme]);
    plot(x_ref, local_error_ref, 'b-d', 'LineWidth', 1.4, 'MarkerSize', 4);
    grid on; xlabel('Axial Position x'); ylabel('|e_i| = |u_{num} - u_{exact}|'); 
    title(sprintf('Pointwise Error Distribution (N = %d, %s)', p.N_ref, current_scheme));
    exportgraphics(fig2, fullfile(out_dir, sprintf('s1_conv_err_%s.pdf', current_scheme)));

    % --- Mesh Convergence Study ---
    if p.FULL
        N_values = [20, 40, 80, 160, 320];
        num_grids = length(N_values);
        dx_arr = zeros(num_grids, 1);
        err_L2 = zeros(num_grids, 1);
        err_Linf = zeros(num_grids, 1);

        for k = 1:num_grids
            [x_k, dx_k, u_num_k] = cfm.solve_advection_1d( ...
                p.L, N_values(k), p.a_val, p.u_inlet, p.S_fn, current_scheme);
            e_k = abs(u_num_k - p.u_exact_fn(x_k));
            dx_arr(k)   = dx_k;
            err_Linf(k) = max(e_k);
            err_L2(k)   = sqrt(mean(e_k.^2));
        end

        % Figure 3: Grid Convergence (Log-Log)
        fig3 = figure('Color', 'w', 'Name', ['Convergence: ' current_scheme]);
        loglog(dx_arr, err_L2, 's-', 'LineWidth', 1.5, 'DisplayName', 'L_2 Norm'); hold on;
        loglog(dx_arr, err_Linf, 'o-', 'LineWidth', 1.5, 'DisplayName', 'L_\infty Norm');
        
        % Plot theoretical slope 1 reference line
        ref_x = [dx_arr(end), dx_arr(1)];
        ref_y = err_L2(end) * (ref_x / ref_x(1)); 
        loglog(ref_x, ref_y, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Theoretical Order 1 Slope');
        
        grid on; xlabel('\Delta x (Grid Spacing)'); ylabel('Global Error Norm'); 
        title(sprintf('Mesh Refinement Study (%s)', current_scheme));
        legend('Location', 'northwest');
        exportgraphics(fig3, fullfile(out_dir, sprintf('s1_conv_mesh_%s.pdf', current_scheme)));

        % Print Convergence Summary Table
        p_L2   = log(err_L2(1:end-1) ./ err_L2(2:end)) ./ log(dx_arr(1:end-1) ./ dx_arr(2:end));
        p_Linf = log(err_Linf(1:end-1) ./ err_Linf(2:end)) ./ log(dx_arr(1:end-1) ./ dx_arr(2:end));

        fprintf('  N       dx         L2 Error      p(L2)     Linf Error    p(Linf)\n');
        fprintf('------------------------------------------------------------------\n');
        for k = 1:num_grids
            if k == 1
                fprintf('%4d   %8.4f   %12.4e     --      %12.4e      --\n', ...
                    N_values(k), dx_arr(k), err_L2(k), err_Linf(k));
            else
                fprintf('%4d   %8.4f   %12.4e   %6.2f    %12.4e   %6.2f\n', ...
                    N_values(k), dx_arr(k), err_L2(k), p_L2(k-1), err_Linf(k), p_Linf(k-1));
            end
        end
    end
end