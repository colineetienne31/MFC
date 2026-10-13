%% Section 2 — finite-volume potential flow
clear; close all; clc;

% =========================================================================
% Parameter Block
% =========================================================================
p.Tinf = 288; 
p.pinf = 101325; 
p.gamma = 1.4; 
p.R = 287;
p.c = 1; 
p.tHalf = 0.03 * p.c;
p.xAirfoil = 1; 
p.xMax = 3; 
p.yMax = 2;
p.Nx = 120; 
p.Ny = 80;

% Solver Settings
p.solver = 'backslash';       % Options: 'backslash', 'gmres', 'bicgstab'
p.tolLinear = 1e-10; 
p.maxitLinear = 1000;
p.tolNonlinear = 1e-7; 
p.maxitNonlinear = 300;
p.relaxation = 0.35;
p.makePlots = true;

% Output directories
outDir = fullfile(pwd, 'figures'); 
if ~exist(outDir, 'dir'), mkdir(outDir); end

resDir = fullfile(pwd, 'results'); 
if ~exist(resDir, 'dir'), mkdir(resDir); end

% =========================================================================
% Regime Selection & Physical Configurations
% =========================================================================
regime = 'subsonic'; % Options: 'subsonic', 'supersonic', 'transonic', 'all'

% Invariant boundary conditions and transverse scheme
p.scheme_y   = 'central';     % Central differences along y
p.bc_inlet   = 'dirichlet';   % Boundary A (Inlet, x = 0): phi = 0
p.bc_top     = 'neumann';     % Boundary B (Top, y = yMax): dphi/dy = 0
p.bc_bottom  = 'tangent';     % Boundary D/E (Bottom / Airfoil): flow tangency

% Regime-dependent physical defaults
if strcmpi(regime, 'subsonic')
    Mlist       = [0.3];
    p.scheme_x  = 'central';
    p.bc_outlet = 'neumann';

elseif strcmpi(regime, 'supersonic')
    Mlist       = [2.2];
    p.scheme_x  = 'backward';
    p.bc_outlet = 'extrapolation';

elseif strcmpi(regime, 'transonic')
    Mlist       = [0.91];
    p.scheme_x  = 'transonic_switch';
    p.bc_outlet = 'neumann';

elseif strcmpi(regime, 'all')
    Mlist       = [0.3, 0.91, 2.2];
end

% -------------------------------------------------------------------------
% Manual overrides for parametric testing
% (Uncomment any line below to override the regime default settings)
% -------------------------------------------------------------------------
% p.scheme_x  = 'backward';       % Force backward/upwind scheme (e.g. test instability in subsonic)
 p.scheme_x  = 'central';        % Force central scheme
% p.scheme_x  = 'transonic_switch';% Force dynamic switch
% p.bc_outlet = 'extrapolation';  % Force free outflow (no boundary constraint imposed)
% p.bc_outlet = 'neumann';        % Force zero normal derivative (dphi/dx = 0)

sol = cell(size(Mlist));

% =========================================================================
% EXECUTION OF TEST CASES (Calls to +cfm module functions)
% =========================================================================
for m = 1:numel(Mlist)
    
    % Só substitui as tuas escolhas se estiveres a correr o regime 'all'
    if strcmpi(regime, 'all')
        if Mlist(m) < 0.8
            p.scheme_x = 'central';
            p.bc_outlet = 'neumann';
        elseif Mlist(m) > 1.2
            p.scheme_x = 'backward';
            p.bc_outlet = 'extrapolation';
        else
            p.scheme_x = 'transonic_switch'; 
            p.bc_outlet = 'neumann';
        end
    end
    
    sol{m} = cfm.solve_flow_fv(p, Mlist(m));
    
    fprintf(['M_inf=%.2f: solver=%s, converged=%d, ' ...
        'iterations=%d, residual=%.3e, computational time=%.3f s\n'], ...
        Mlist(m), ...
        sol{m}.solver, ...
        sol{m}.converged, ...
        sol{m}.iterations, ...  
        sol{m}.residual, ...
        sol{m}.computationalTime);
end


% =========================================================================
% POST-PROCESSING AND GRAPHICS
% =========================================================================
if p.makePlots
    
    for m = 1:numel(Mlist)
        s = sol{m};
        figure('Color', 'w', 'Position', [100, 200, 1200, 350]);
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % 1. Potential
        nexttile;
        contourf(s.X, s.Y, s.phiNorm, 30, 'LineColor', 'none'); colorbar;
        hold on; cfm.plot_airfoil(p); hold off;
        axis equal tight; xlim([0 p.xMax]); ylim([0 p.yMax]);
        xlabel('x (m)'); ylabel('y (m)');
        title('\phi/(U_\infty c)');
        
        % 2. Local Mach
        nexttile;
        contourf(s.X, s.Y, s.Mlocal, 30, 'LineColor', 'none'); colorbar;
        hold on; cfm.plot_airfoil(p); hold off;
        axis equal tight; xlim([0 p.xMax]); ylim([0 p.yMax]);
        xlabel('x (m)'); ylabel('y (m)'); title('Local Mach number');
        
        % 3. Relative density 
        nexttile;
        contourf(s.X, s.Y, s.rhoRatio, 30, 'LineColor', 'none'); colorbar;
        hold on; cfm.plot_airfoil(p); hold off;
        axis equal tight; xlim([0 p.xMax]); ylim([0 p.yMax]);
        xlabel('x (m)'); ylabel('y (m)'); title('\rho/\rho_\infty');
        
        sgtitle(sprintf('Finite-volume potential flow, M_\\infty=%.2f', Mlist(m)));
        
        % Export figure
        exportgraphics(gcf, fullfile(outDir, sprintf('s2_flow_M_%g.pdf', Mlist(m))));
    end
end


% =========================================================================
% SUMMARY TABLE
% =========================================================================
nCases = numel(Mlist);

Regime = strings(nCases, 1);
Fluxes = strings(nCases, 1);
Switch = strings(nCases, 1);
BoundaryA = strings(nCases, 1);
BoundaryB = repmat("dphi/dy=0", nCases, 1);
BoundaryC = strings(nCases, 1);
BoundaryD = repmat("dphi/dy=0", nCases, 1);
BoundaryE = repmat("dphi/dy=Uinf*dy_s/dx", nCases, 1);
Solver = strings(nCases, 1);
Convergence = strings(nCases, 1);
Iterations = zeros(nCases, 1);
Residual = zeros(nCases, 1);
Computational_time = zeros(nCases, 1);

for m = 1:nCases

    % Identify flow regime
    if Mlist(m) < 0.8
        Regime(m) = "Subsonic";
        Fluxes(m) = "Central x, central y";
        Switch(m) = "None";
        BoundaryA(m) = "phi=0";
        BoundaryC(m) = "dphi/dx=0";

    elseif Mlist(m) > 1.2
        Regime(m) = "Supersonic";
        Fluxes(m) = "Backward x, central y";
        Switch(m) = "None";
        BoundaryA(m) = "phi=0";
        BoundaryC(m) = "No imposed condition";

    else
        Regime(m) = "Transonic";
        Fluxes(m) = ...
            "Central x if A>=0, backward x if A<0; central y";
        Switch(m) = ...
            "A=1-Minf^2-(gamma+1)*Minf^2/Uinf*dphi/dx";
        BoundaryA(m) = "phi=0";
        BoundaryC(m) = "dphi/dx=0";
    end

    % Solver and performance results
    Solver(m) = string(sol{m}.solver);
    Convergence(m) = string(sol{m}.converged);
    Iterations(m) = sol{m}.iterations;
    Residual(m) = sol{m}.residual;
    Computational_time(m) = sol{m}.computationalTime;
end

Mach = Mlist(:);

summary = table(Regime, Mach, Fluxes, Switch, ...
    BoundaryA, BoundaryB, BoundaryC, BoundaryD, BoundaryE, ...
    Solver, Convergence, Iterations, Residual, Computational_time);

% Export results
writetable(summary, ...
    fullfile(resDir, 's2_method_boundary_solver.csv'));

disp(summary);
