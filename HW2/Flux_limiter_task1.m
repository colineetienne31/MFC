%% Flux limiter: own design + classical limiters, Sweby diagram and extras (standalone)
clear; close all; clc;
 
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);
 
%% Parameters (same as Task_1.m)
L = 200;
xmin = -100;
xmax = 100;
dx = 1;
a = 1;
nu = 0.8;
N = round(L/dx);
dt = nu*dx/a;
T1 = L/a;
T2 = 2*L/a;
N1_steps = round(T1/dt);
N2_steps = round(T2/dt);
b_blend = 0.6; % 0.6 at nu=0.8: leading phase-error cancellation, also used by the own limiter
 
%% Output folders (inside a Figures folder next to this .m file)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
output_dir = fullfile(script_dir, 'Figures');
folders.solutions = fullfile(output_dir, 'Solutions');
folders.sim_errors = fullfile(output_dir, 'Simulation_Errors');
folders.sweby = fullfile(output_dir, 'Sweby_Diagram');
folders.limiter = fullfile(output_dir, 'Flux_Limiter');
folder_names = fieldnames(folders);
for j = 1:numel(folder_names)
    folder = folders.(folder_names{j});
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end
 
%% Grid and initial conditions
x = (xmin + ((1:N)-0.5)*dx)'; % column of cell centers
u_acoustic = 0.5*exp(-log(2)*(x/3).^2);
u_square = 0.5*double(abs(x)<=10);
initial_conditions = {u_acoustic, u_square};
pulse_names = {'Acoustic pulse', 'Square pulse'};
 
%% Flux limiters: psi(r), r = (u_i - u_{i-1})/(u_{i+1} - u_i)
% Own limiter: TVD-limited version of the optimal blend psi = (1-b) + b*r
% (LW for b=0, BW for b=1, Fromm for b=0.5), clipped to the Sweby TVD region.
limiter_names = {'Minmod', 'Van Leer', 'Van Albada', 'MC', 'Superbee', ...
                 sprintf('Own limiter (b = %.2f)',b_blend)};
psi_fun = {@psi_minmod, @psi_vanleer, @psi_vanalbada, @psi_mc, ...
           @psi_superbee, @(r) psi_own(r,b_blend)};
n_lim = numel(psi_fun);
scheme_names = [{'Explicit Upwind', sprintf('LW-BW blend (b = %.2f)',b_blend)}, ...
                limiter_names];
n_schemes = numel(scheme_names); % 1 = upwind, 2 = blend, 3..8 = limiters
 
%% Sweby diagram
r = linspace(0, 4, 801);
fig_sweby = square_figure();
ax = axes('Parent',fig_sweby);
hold(ax,'on');
% TVD region: 0 <= psi <= min(2r,2)
fill(ax,[0 1 4 4],[0 2 2 0],[0.85 0.9 1],'EdgeColor','none','FaceAlpha',0.6);
% Second-order TVD region: between psi = r (BW) and psi = 1 (LW)
fill(ax,[0 0.5 1],[0 1 1],[1 0.8 0.6],'EdgeColor','none','FaceAlpha',0.8);
fill(ax,[1 2 4 4],[1 2 2 1],[1 0.8 0.6],'EdgeColor','none','FaceAlpha',0.8);
plot(ax,r,ones(size(r)),'k--','LineWidth',1);   % Lax-Wendroff
plot(ax,r,r,'k:','LineWidth',1);                 % Beam-Warming
sweby_colors = lines(n_lim);
h_sweby = gobjects(1,n_lim);
for j = 1:n_lim
    lw = 1.5 + 1.5*(j == n_lim);
    h_sweby(j) = plot(ax,r,psi_fun{j}(r),'-','Color',sweby_colors(j,:),'LineWidth',lw);
end
xlim(ax,[0 4]); ylim(ax,[0 3]);
xlabel(ax,'r'); ylabel(ax,'\psi(r)');
title(ax,'Sweby diagram');
legend(ax,h_sweby,limiter_names,'Location','northwest');
grid(ax,'on'); hold(ax,'off');
save_figure(fig_sweby,folders.sweby,'Sweby_Diagram');
close(fig_sweby);
 
%% Solutions, simulation-error evolution and total variation
time = (0:N2_steps)'*dt;
for p = 1:2
    u0 = initial_conditions{p}(:);
    U = repmat(u0,1,n_schemes);
    error_L1 = zeros(N2_steps+1,n_schemes);
    error_Linf = zeros(N2_steps+1,n_schemes);
    TV = zeros(N2_steps+1,n_schemes);
    TV0 = sum(abs(u0-circshift(u0,1)));
    if p == 1
        pulse_file = 'Acoustic';
    else
        pulse_file = 'Square';
    end
 
    for n = 0:N2_steps
        t = n*dt;
        u_exact = exact_solution(x,t,p,a,L,xmin);
 
        for s = 1:n_schemes
            difference = abs(U(:,s)-u_exact);
            error_L1(n+1,s) = mean(difference);
            error_Linf(n+1,s) = max(difference);
            TV(n+1,s) = sum(abs(U(:,s)-circshift(U(:,s),1)));
        end
 
        if n == N1_steps || n == N2_steps
            rev = 1 + (n == N2_steps);
 
            % Reference plot: upwind, LW-BW blend and own limiter
            fig_ref = square_figure();
            ax = axes('Parent',fig_ref);
            hold(ax,'on');
            plot(ax,x,u_exact,'k--','LineWidth',2);
            ref_columns = [1 2 n_schemes];
            ref_colors = lines(3);
            for j = 1:3
                plot(ax,x,U(:,ref_columns(j)),'-','Color',ref_colors(j,:), ...
                    'LineWidth',1.5);
            end
            xlim(ax,[-35,35]);
            xlabel(ax,'x'); ylabel(ax,'u');
            title(ax,sprintf('%s - own limiter - %d revolution(s), t = %.0f', ...
                pulse_names{p},rev,t));
            legend(ax,[{'Exact solution'},scheme_names(ref_columns)],'Location','best');
            grid(ax,'on'); hold(ax,'off');
            save_figure(fig_ref,folders.limiter, ...
                sprintf('Own_Limiter_%s_%d_Revolution',pulse_file,rev));
            close(fig_ref);
 
            % Classical limiters versus own limiter
            fig_lim = square_figure();
            ax = axes('Parent',fig_lim);
            hold(ax,'on');
            plot(ax,x,u_exact,'k--','LineWidth',2);
            lim_colors = lines(n_lim);
            for j = 1:n_lim
                plot(ax,x,U(:,j+2),'-','Color',lim_colors(j,:), ...
                    'LineWidth',1.3 + 0.7*(j == n_lim));
            end
            xlim(ax,[-35,35]);
            xlabel(ax,'x'); ylabel(ax,'u');
            title(ax,sprintf('%s - flux limiters - %d revolution(s), t = %.0f', ...
                pulse_names{p},rev,t));
            legend(ax,[{'Exact solution'},limiter_names],'Location','best');
            grid(ax,'on'); hold(ax,'off');
            save_figure(fig_lim,folders.solutions, ...
                sprintf('Solutions_Limiters_%s_%d_Revolution',pulse_file,rev));
            close(fig_lim);
        end
 
        if n == N2_steps
            break;
        end
        U(:,1) = explicit_upwind(U(:,1),nu);
        U(:,2) = lw_bw_blend(U(:,2),nu,b_blend);
        for j = 1:n_lim
            U(:,j+2) = flux_limiter_step(U(:,j+2),nu,psi_fun{j});
        end
    end
 
    % Error evolution and total variation (all schemes)
    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_L1,'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Mean absolute error');
    title(ax,[pulse_names{p},' - limiters L1 error evolution']);
    legend(ax,scheme_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_L1_Limiters_',pulse_file]);
    close(fig_errors);
 
    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_Linf,'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Maximum absolute error');
    title(ax,[pulse_names{p},' - limiters Linf error evolution']);
    legend(ax,scheme_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_Linf_Limiters_',pulse_file]);
    close(fig_errors);
 
    fig_tv = square_figure();
    ax = axes('Parent',fig_tv);
    hold(ax,'on');
    plot(ax,time,TV,'LineWidth',1.5);
    yline(ax,TV0,'k--','HandleVisibility','off');
    xlabel(ax,'Time'); ylabel(ax,'Total variation TV(t)');
    title(ax,[pulse_names{p},' - total variation']);
    legend(ax,scheme_names,'Location','best'); grid(ax,'on');
    hold(ax,'off');
    save_figure(fig_tv,folders.sim_errors,['Total_Variation_',pulse_file]);
    close(fig_tv);
 
    % Summary table at t = 2 revolutions (overshoot/undershoot = monotonicity check)
    overshoot = (max(U,[],1) - max(u0))';
    undershoot = (-min(U,[],1))';
    summary = table(string(scheme_names'), error_L1(end,:)', error_Linf(end,:)', ...
        overshoot, undershoot, TV(end,:)', ...
        'VariableNames',{'Scheme','L1','Linf','Overshoot','Undershoot','TV'});
    fprintf('\n%s - errors after 2 revolutions (TV exact = %.3f)\n',pulse_names{p},TV0);
    disp(summary);
    writetable(summary,fullfile(folders.limiter,['Summary_Limiters_',pulse_file,'.csv']));
end
 
%% Extra 1: grid refinement (observed order), nu = 0.8, 1 revolution
dx_values = [2, 1, 0.5, 0.25];
grid_names = {'Explicit Upwind', 'LW-BW blend', 'MC', 'Own limiter'};
grid_steps = {@(u,nu_) explicit_upwind(u,nu_), ...
              @(u,nu_) lw_bw_blend(u,nu_,b_blend), ...
              @(u,nu_) flux_limiter_step(u,nu_,psi_fun{4}), ...
              @(u,nu_) flux_limiter_step(u,nu_,psi_fun{6})};
err_grid = zeros(numel(dx_values),numel(grid_steps),2);
for d = 1:numel(dx_values)
    dx_g = dx_values(d);
    N_g = round(L/dx_g);
    dt_g = nu*dx_g/a;
    nsteps_g = round(T1/dt_g);
    x_g = (xmin + ((1:N_g)-0.5)*dx_g)';
    for p = 1:2
        u_g0 = exact_solution(x_g,0,p,a,L,xmin);
        for s = 1:numel(grid_steps)
            u_g = u_g0;
            for n = 1:nsteps_g
                u_g = grid_steps{s}(u_g,nu);
            end
            err_grid(d,s,p) = mean(abs(u_g-exact_solution(x_g,nsteps_g*dt_g,p,a,L,xmin)));
        end
    end
end
 
for p = 1:2
    if p == 1
        pulse_file = 'Acoustic';
    else
        pulse_file = 'Square';
    end
    fig_grid = square_figure();
    ax = axes('Parent',fig_grid);
    hold(ax,'on');
    grid_colors = lines(numel(grid_steps));
    grid_legend = cell(1,numel(grid_steps));
    for s = 1:numel(grid_steps)
        fit = polyfit(log(dx_values),log(err_grid(:,s,p)'),1);
        loglog(ax,dx_values,err_grid(:,s,p),'o-','Color',grid_colors(s,:), ...
            'LineWidth',1.5,'MarkerSize',8);
        grid_legend{s} = sprintf('%s (order %.2f)',grid_names{s},fit(1));
        fprintf('%s - %s: observed L1 order = %.2f\n',pulse_names{p},grid_names{s},fit(1));
    end
    set(ax,'XScale','log','YScale','log');
    xlabel(ax,'\Deltax'); ylabel(ax,'Mean absolute error (1 revolution)');
    title(ax,[pulse_names{p},' - grid refinement']);
    legend(ax,grid_legend,'Location','southeast'); grid(ax,'on');
    hold(ax,'off');
    save_figure(fig_grid,folders.limiter,['Grid_Refinement_',pulse_file]);
    close(fig_grid);
end
 
%% Extra 2: Courant number sweep for the own limiter (square pulse, 1 revolution)
nu_sweep = [0.2, 0.5, 0.8, 1];
fig_nu = square_figure();
ax = axes('Parent',fig_nu);
hold(ax,'on');
plot(ax,x,exact_solution(x,T1,2,a,L,xmin),'k--','LineWidth',2);
nu_colors = lines(numel(nu_sweep));
for k = 1:numel(nu_sweep)
    nu_k = nu_sweep(k);
    dt_k = nu_k*dx/a;
    nsteps_k = round(T1/dt_k);
    u_k = u_square;
    for n = 1:nsteps_k
        u_k = flux_limiter_step(u_k,nu_k,psi_fun{n_lim});
    end
    plot(ax,x,u_k,'-','Color',nu_colors(k,:),'LineWidth',1.5);
end
xlim(ax,[-35,35]);
xlabel(ax,'x'); ylabel(ax,'u');
title(ax,'Own limiter - Courant number sweep, 1 revolution');
legend(ax,[{'Exact solution'},arrayfun(@(v) sprintf('\\nu = %.1f',v),nu_sweep, ...
    'UniformOutput',false)],'Location','best');
grid(ax,'on'); hold(ax,'off');
save_figure(fig_nu,folders.limiter,'Own_Limiter_Courant_Sweep');
close(fig_nu);
 
fprintf('Flux limiter plots saved under: %s\n', output_dir);
 
%% Local functions
function u_new = explicit_upwind(u,nu)
    u_new = u-nu*(u-circshift(u,1));
end
 
% Conservative periodic Lax-Wendroff update
function u_new = lax_wendroff(u,nu)
    up = circshift(u,-1);
    um = circshift(u,1);
    u_new = u - (nu/2)*(up-um) ...
              + (nu^2/2)*(up-2*u+um);
end
 
% Conservative periodic Beam-Warming update, a > 0
function u_new = beam_warming(u,nu)
    um = circshift(u,1);
    um2 = circshift(u,2);
    u_new = u - (nu/2)*(3*u-4*um+um2) ...
              + (nu^2/2)*(u-2*um+um2);
end
 
% LW-BW constant-b blend
function u_new = lw_bw_blend(u,nu,b)
    u_new = (1-b)*lax_wendroff(u,nu) + b*beam_warming(u,nu);
end
 
% Flux-limited update, a > 0 (TVD for 0 < nu <= 1):
% F_{i+1/2}/a = u_i + (1-nu)/2 * psi(r_{i+1/2}) * (u_{i+1}-u_i)
% psi = 0 -> upwind, psi = 1 -> Lax-Wendroff, psi = r -> Beam-Warming
function u_new = flux_limiter_step(u,nu,psi)
    up = circshift(u,-1);
    um = circshift(u,1);
    dp = up-u;
    dm = u-um;
    F = u + 0.5*(1-nu)*psi(limiter_ratio(dm,dp)).*dp;
    u_new = u - nu*(F-circshift(F,1));
end
 
% Ratio of consecutive gradients; safe against division by zero
function r = limiter_ratio(dm,dp)
    den = dp;
    den(abs(den)<1e-12) = 1e-12;
    r = dm./den;
end
 
% Classical limiters
function psi = psi_minmod(r)
    psi = max(0,min(1,r));
end
 
function psi = psi_vanleer(r)
    psi = (r+abs(r))./(1+abs(r));
end
 
function psi = psi_vanalbada(r)
    psi = max(0,(r.^2+r)./(r.^2+1));
end
 
function psi = psi_mc(r)
    psi = max(0,min(min(2*r,(1+r)/2),2));
end
 
function psi = psi_superbee(r)
    psi = max(0,max(min(2*r,1),min(r,2)));
end
 
% Own limiter: blend (1-b)*LW + b*BW, clipped to the TVD region
function psi = psi_own(r,b)
    psi = max(0,min(min(2*r,(1-b)+b*r),2));
end
 
function u_exact = exact_solution(x,t,p,a,L,xmin)
    x_shifted = mod(x-a*t-xmin,L)+xmin;
    if p == 1
        u_exact = 0.5*exp(-log(2)*(x_shifted/3).^2);
    else
        u_exact = 0.5*double(abs(x_shifted)<=10);
    end
    u_exact = u_exact(:);
end
 
function fig = square_figure()
    fig = figure('Visible','off','Color','w', ...
        'Units','pixels','Position',[100 100 900 900], ...
        'PaperPositionMode','auto');
end
 
function save_figure(fig,folder,basename)
    exportgraphics(fig,fullfile(folder,[basename,'.png']),'Resolution',300);
    savefig(fig,fullfile(folder,[basename,'.fig']));
end