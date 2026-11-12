%% LW-BW blend: stability limit nu_max(b) and choice of b (standalone)
clear; close all; clc;
 
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 19);
set(groot, 'defaultLegendFontSize', 19);
 
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
b_blend = 0.6; % value used in Task_1.m, checked against the analysis below
tol_stab = 1e-10; % tolerance on max|G| <= 1
 
%% Output folders (inside a Figures folder next to this .m file)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
output_dir = fullfile(script_dir, 'Figures');
folders.stability = fullfile(output_dir, 'Stability_Diagrams');
folders.blend = fullfile(output_dir, 'Blend_Selection');
folders.solutions = fullfile(output_dir, 'Solutions');
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
 
%% Stability limit: scan in (nu, b)
% Stable <=> max_theta |G| <= 1. |G(-theta)| = |G(theta)|, so theta in [0, pi] suffices.
theta = linspace(0, pi, 1001);
nu_scan = 0:0.005:2.5;
b_scan = 0:0.01:1;
G_max = zeros(numel(nu_scan),numel(b_scan));
for i = 1:numel(nu_scan)
    [G_LW,G_BW] = amplification_LW_BW(nu_scan(i),theta);
    G_all = (1-b_scan')*G_LW + b_scan'*G_BW; % rows: b, columns: theta
    G_max(i,:) = max(abs(G_all),[],2)';
end
 
% Stability limit = largest nu such that the scheme is stable for all smaller nu
% (coarse scan, then bisection for the exact value)
nu_max = zeros(size(b_scan));
for j = 1:numel(b_scan)
    first_unstable = find(G_max(:,j) > 1+tol_stab, 1);
    if isempty(first_unstable)
        nu_max(j) = nu_scan(end);
    elseif first_unstable == 1
        nu_max(j) = 0;
    else
        nu_max(j) = stability_limit(b_scan(j),theta, ...
            nu_scan(first_unstable-1),nu_scan(first_unstable),tol_stab);
    end
end
 
[~,j_blend] = min(abs(b_scan-b_blend));
fprintf('Stability limit: LW (b=0): nu_max = %.3f | BW (b=1): nu_max = %.3f\n', ...
    nu_max(1),nu_max(end));
fprintf('Stability limit for b = %.2f: nu_max = %.3f (used nu = %.2f)\n', ...
    b_scan(j_blend),nu_max(j_blend),nu);
if nu > nu_max(j_blend)
    warning('nu = %.2f is outside the stability region for b = %.2f',nu,b_blend);
end
 
% Map of max|G| (stable region: max|G| = 1) with the stability limit on top
fig_map = square_figure();
ax = axes('Parent',fig_map);
imagesc(ax,b_scan,nu_scan,min(G_max,1.2));
set(ax,'YDir','normal');
hold(ax,'on');
plot(ax,b_scan,nu_max,'w-','LineWidth',2);
plot(ax,b_blend,nu,'wo','MarkerFaceColor','r','MarkerSize',9);
cb = colorbar(ax);
ylabel(cb,'max_\theta |G| (saturated at 1.2)');
xlabel(ax,'b'); ylabel(ax,'\nu');
title(ax,'LW-BW stability map');
hold(ax,'off');
save_figure(fig_map,folders.stability,'Stability_Map_nu_b');
close(fig_map);
 
% Stability limit curve nu_max(b)
fig_limit = square_figure();
ax = axes('Parent',fig_limit);
hold(ax,'on');
plot(ax,b_scan,nu_max,'LineWidth',1.5);
plot(ax,b_blend,nu_max(j_blend),'ro','MarkerFaceColor','r','MarkerSize',8);
yline(ax,nu,'k--','HandleVisibility','off');
xlabel(ax,'b'); ylabel(ax,'\nu_{max}');
title(ax,'LW-BW stability limit');
legend(ax,{'\nu_{max}(b)',sprintf('b = %.2f',b_blend)},'Location','best');
grid(ax,'on'); hold(ax,'off');
save_figure(fig_limit,folders.stability,'Stability_Limit_nu_max');
close(fig_limit);
 
%% Choice of b from the error analysis (nu fixed)
% eps_D = |G|, eps_phi = -arg(G)/(nu*theta); ideal value = 1 for both.
theta_fine = linspace(0, pi, 2001);
theta_cut = pi/2; % band of interest for the smooth (acoustic) pulse
in_band = theta_fine <= theta_cut;
b_fine = 0:0.005:1;
J_diff = zeros(size(b_fine));
J_disp = zeros(size(b_fine));
J_pred = zeros(numel(b_fine),2); % spectrally weighted RMS error after 1 revolution
 
% Spectral weights |u_hat(k)|^2 of the two pulses (k = theta/dx)
k_wave = theta_fine/dx;
w_spectrum = zeros(2,numel(theta_fine));
w_spectrum(1,:) = exp(-9*k_wave.^2/(2*log(2)));       % Gaussian exp(-ln2 (x/3)^2)
w_spectrum(2,:) = (sin(10*k_wave)./k_wave).^2;        % square |x| <= 10
w_spectrum(2,1) = 100;                                % limit k -> 0
 
G_exact = exp(-1i*nu*theta_fine);
[G_LW,G_BW] = amplification_LW_BW(nu,theta_fine);
for j = 1:numel(b_fine)
    G = (1-b_fine(j))*G_LW + b_fine(j)*G_BW;
    E_diff = abs(G);
    E_disp = -unwrap(angle(G))./(nu*theta_fine);
    E_disp(1) = 1; % continuous low-wavenumber limit
    J_diff(j) = trapz(theta_fine(in_band),abs(1-E_diff(in_band)));
    J_disp(j) = trapz(theta_fine(in_band),abs(1-E_disp(in_band)));
    for p = 1:2
        err_spec = abs(G.^N1_steps - G_exact.^N1_steps).^2;
        J_pred(j,p) = sqrt(trapz(theta_fine,w_spectrum(p,:).*err_spec) ...
                          /trapz(theta_fine,w_spectrum(p,:)));
    end
end
b_star = (1+nu)/3; % theoretical value (modified equation), to be verified
[~,i_diff] = min(J_diff);
[~,i_disp] = min(J_disp);
[~,i_pred1] = min(J_pred(:,1));
[~,i_pred2] = min(J_pred(:,2));
fprintf('\nChoice of b at nu = %.2f:\n',nu);
fprintf('  Theory (modified equation) b* = (1+nu)/3 : %.3f\n',b_star);
fprintf('  Min dispersive error (theta <= %.0f deg) : %.3f\n',rad2deg(theta_cut),b_fine(i_disp));
fprintf('  Min diffusive error  (theta <= %.0f deg) : %.3f\n',rad2deg(theta_cut),b_fine(i_diff));
fprintf('  Min predicted error, acoustic spectrum   : %.3f\n',b_fine(i_pred1));
fprintf('  Min predicted error, square spectrum     : %.3f\n',b_fine(i_pred2));
 
% Band-limited diffusive and dispersive error versus b
fig_b = square_figure();
ax = axes('Parent',fig_b);
hold(ax,'on');
plot(ax,b_fine,J_disp,'LineWidth',1.5);
plot(ax,b_fine,J_diff,'LineWidth',1.5);
xline(ax,b_star,'k--','HandleVisibility','off');
xlabel(ax,'b');
ylabel(ax,sprintf('\\int_0^{%.0f^\\circ} |1-\\epsilon|\\, d\\theta',rad2deg(theta_cut)));
title(ax,sprintf('Error versus b - nu = %.2f',nu));
legend(ax,{'Dispersive \epsilon_\phi','Diffusive \epsilon_D'},'Location','best');
grid(ax,'on'); hold(ax,'off');
save_figure(fig_b,folders.blend,'Error_Analysis_vs_b');
close(fig_b);
 
% Spectrally weighted prediction versus b
fig_pred = square_figure();
ax = axes('Parent',fig_pred);
hold(ax,'on');
plot(ax,b_fine,J_pred(:,1),'LineWidth',1.5);
plot(ax,b_fine,J_pred(:,2),'LineWidth',1.5);
xline(ax,b_star,'k--','HandleVisibility','off');
xlabel(ax,'b'); ylabel(ax,'Predicted RMS error after 1 revolution');
title(ax,sprintf('Spectral prediction versus b - nu = %.2f',nu));
legend(ax,pulse_names,'Location','best');
grid(ax,'on'); hold(ax,'off');
save_figure(fig_pred,folders.blend,'Spectral_Prediction_vs_b');
close(fig_pred);
 
%% Verification with simulations: error versus b after 1 and 2 revolutions
b_sim = 0:0.05:1;
err_sim = zeros(numel(b_sim),4,2); % columns: L1 (1 rev), L1 (2 rev), Linf (1 rev), Linf (2 rev)
for p = 1:2
    u0 = initial_conditions{p}(:);
    for j = 1:numel(b_sim)
        u = u0;
        for n = 1:N2_steps
            u = lw_bw_blend(u,nu,b_sim(j));
            if n == N1_steps || n == N2_steps
                diff_n = abs(u-exact_solution(x,n*dt,p,a,L,xmin));
                col = 1 + (n == N2_steps);
                err_sim(j,col,p) = mean(diff_n);
                err_sim(j,col+2,p) = max(diff_n);
            end
        end
    end
 
    if p == 1
        pulse_file = 'Acoustic';
    else
        pulse_file = 'Square';
    end
    [~,i_best] = min(err_sim(:,2,p));
    fprintf('%s: simulated best b (L1, 2 rev) = %.2f\n',pulse_names{p},b_sim(i_best));
 
    fig_sim = square_figure();
    ax = axes('Parent',fig_sim);
    semilogy(ax,b_sim,err_sim(:,:,p),'o-','LineWidth',1.5);
    hold(ax,'on');
    xline(ax,b_star,'k--','HandleVisibility','off');
    xlabel(ax,'b'); ylabel(ax,'Absolute error');
    title(ax,[pulse_names{p},' - LW-BW error versus b']);
    legend(ax,{'L1, 1 rev','L1, 2 rev','Linf, 1 rev','Linf, 2 rev'},'Location','best');
    grid(ax,'on'); hold(ax,'off');
    save_figure(fig_sim,folders.blend,['Simulation_Error_vs_b_',pulse_file]);
    close(fig_sim);
end
 
%% Solutions for several b after 2 revolutions
b_plot = [0, 0.4, b_star, 0.8, 1];
for p = 1:2
    u_exact = exact_solution(x,N2_steps*dt,p,a,L,xmin);
    fig_b_sol = square_figure();
    ax = axes('Parent',fig_b_sol);
    hold(ax,'on');
    plot(ax,x,u_exact,'k--','LineWidth',2);
    b_colors = lines(numel(b_plot));
    for j = 1:numel(b_plot)
        u = initial_conditions{p}(:);
        for n = 1:N2_steps
            u = lw_bw_blend(u,nu,b_plot(j));
        end
        plot(ax,x,u,'-','Color',b_colors(j,:),'LineWidth',1.5);
    end
    xlim(ax,[-35,35]);
    xlabel(ax,'x'); ylabel(ax,'u');
    title(ax,sprintf('%s - LW-BW, 2 revolutions, t = %.0f',pulse_names{p},N2_steps*dt));
    legend(ax,[{'Exact solution'},arrayfun(@(b) sprintf('b = %.2f',b),b_plot, ...
        'UniformOutput',false)],'Location','best');
    grid(ax,'on'); hold(ax,'off');
    if p == 1
        pulse_file = 'Acoustic';
    else
        pulse_file = 'Square';
    end
    save_figure(fig_b_sol,folders.solutions,['Solutions_LW_BW_b_Sweep_',pulse_file]);
    close(fig_b_sol);
end
 
fprintf('Blend selection plots saved under: %s\n', output_dir);
 
%% Local functions
% Amplification factors of Lax-Wendroff and Beam-Warming
function [G_LW,G_BW] = amplification_LW_BW(nu,theta)
    G_LW = 1 - 1i*nu*sin(theta) + nu^2*(cos(theta)-1);
    G_BW = 1 - (nu/2)*(3 - 4*exp(-1i*theta) + exp(-2i*theta)) ...
             + (nu^2/2)*(1 - 2*exp(-1i*theta) + exp(-2i*theta));
end
 
% Bisection between a stable (nu_lo) and an unstable (nu_hi) Courant number
function nu_lim = stability_limit(b,theta,nu_lo,nu_hi,tol)
    for it = 1:40
        nu_mid = 0.5*(nu_lo+nu_hi);
        [G_LW,G_BW] = amplification_LW_BW(nu_mid,theta);
        G = (1-b)*G_LW + b*G_BW;
        if max(abs(G)) > 1+tol
            nu_hi = nu_mid;
        else
            nu_lo = nu_mid;
        end
    end
    nu_lim = nu_lo;
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
