%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Fluid Mechanics - Homework 1
% Section 1.2: Steady-State Convection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Create figures directory if it doesn't exist
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Parameter Block
L = 40;                 % Domain length
a_val = 1.0;            % Convection velocity
u_inlet = 1.0;          % Inlet Boundary Condition at x=0
N_ref = 40;             % Baseline grid size
scheme_x = 'backward';  % Scheme selection by name ('backward' or 'central')
FULL = true;            % Flag to toggle the time-consuming mesh convergence study

% Exact solution and corresponding source term
u_exact_fn = @(x) exp(sin(0.2 * x));
S_fn       = @(x) a_val * (0.2 * (cos (0.2 * x)) .* exp(sin(0.2 * x)));

%% Baseline Grid Solution (N = 40)
[x_ref, dx_ref, u_num_ref] = solve_advection_1d(L, N_ref, a_val, u_inlet, S_fn, scheme_x);

u_exact_ref = u_exact_fn(x_ref);
local_error_ref = abs(u_num_ref - u_exact_ref);

%% Deliverable Plots for Baseline
% Directory creation
out_dir = fullfile(pwd, 'figures');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Figure 1: Numerical vs Exact Solution
fig1 = figure('Color', 'w');
plot(x_ref, u_exact_ref, 'g-', 'LineWidth', 1.8, 'DisplayName', 'Analytical'); hold on;
plot(x_ref, u_num_ref, 'ro--', 'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', ['Numerical']);
grid on; xlabel('Axial Position x'); ylabel('u(x)'); title('Numerical Solution vs Analytical Reference (N = 40)');
legend('Location', 'best');

%Graphics Export Figure 1
exportgraphics(fig1, fullfile(out_dir, 's1_convection_solution.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(out_dir, 's1_convection_solution.jpg'), 'Resolution', 300);

% Figure 2: Local Absolute Error Distribution
fig2 = figure('Color', 'w');
plot(x_ref, local_error_ref, 'b-d', 'LineWidth', 1.4, 'MarkerSize', 4);
grid on; xlabel('Axial Position x'); ylabel('|e_i| = |u_{num} - u_{exact}|'); title('Local Pointwise Error Distribution (N = 40)');

%Graphics Export Figure 2
exportgraphics(fig2, fullfile(out_dir, 's1_convection_error.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(out_dir, 's1_convection_error.jpg'), 'Resolution', 300);

%% Mesh Convergence Study
if FULL
    N_values = [20, 40, 80, 160, 320];
    num_grids = length(N_values);
    dx_arr    = zeros(num_grids, 1);
    err_L2    = zeros(num_grids, 1);
    err_Linf  = zeros(num_grids, 1);

    for k = 1:num_grids
        [x_k, dx_k, u_num_k] = solve_advection_1d(L, N_values(k), a_val, u_inlet, S_fn, scheme_x);
        e_k = abs(u_num_k - u_exact_fn(x_k));

        dx_arr(k)   = dx_k;
        err_Linf(k) = max(e_k);
        err_L2(k)   = sqrt(mean(e_k.^2));
    end

    % Figure 3: Grid Convergence (Log-Log)
    fig3 = figure('Color', 'w');
    loglog(dx_arr, err_L2, 's-', 'LineWidth', 1.5, 'DisplayName', 'L_2 Norm'); hold on;
    loglog(dx_arr, err_Linf, 'o-', 'LineWidth', 1.5, 'DisplayName', 'L_\infty Norm');
    ref_x = [dx_arr(end), dx_arr(1)];
    ref_y = err_L2(end) * (ref_x / ref_x(1)); % Slope = 1
    loglog(ref_x, ref_y, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Theoretical Order 1 Slope');
    grid on; xlabel('\Delta x (Grid Spacing)'); ylabel('Global Error Norm'); title('Mesh Refinement Study (Log-Log)');
    legend('Location', 'northwest');
    
    %Graphics Export Figure 3
    exportgraphics(fig3, fullfile(out_dir, 's1_convection_convergence.pdf'), 'ContentType', 'vector');
    exportgraphics(fig3, fullfile(out_dir, 's1_convection_convergence.jpg'), 'Resolution', 300);

    % Compute empirical convergence order p between consecutive grids
    p_L2   = log(err_L2(1:end-1) ./ err_L2(2:end)) ./ log(dx_arr(1:end-1) ./ dx_arr(2:end));
    p_Linf = log(err_Linf(1:end-1) ./ err_Linf(2:end)) ./ log(dx_arr(1:end-1) ./ dx_arr(2:end));

    fprintf('================ MESH CONVERGENCE SUMMARY ================\n');
    fprintf('  N       dx         L2 Error      p(L2)     Linf Error    p(Linf)\n');
    fprintf('----------------------------------------------------------\n');
    for k = 1:num_grids
        if k == 1
            fprintf('%4d   %8.4f   %12.4e     --      %12.4e      --\n', ...
                N_values(k), dx_arr(k), err_L2(k), err_Linf(k));
        else
            fprintf('%4d   %8.4f   %12.4e   %6.2f    %12.4e   %6.2f\n', ...
                N_values(k), dx_arr(k), err_L2(k), p_L2(k-1), err_Linf(k), p_Linf(k-1));
        end
    end
    fprintf('==========================================================\n');

end