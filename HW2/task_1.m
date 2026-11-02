clear; close all; clc;

set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);

%% Each plot is its own invisible 900x900 figure, saved to folders, and closed.
%% Parameters
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
b_blend = (1+nu)/3; % 0.6 at nu=0.8: leading phase-error cancellation

%% Output folders (inside a Figures folder next to this .m file)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
output_dir = fullfile(script_dir, 'Figures');
folders.solutions = fullfile(output_dir, 'Solutions');
folders.sim_errors = fullfile(output_dir, 'Simulation_Errors');
folders.stability = fullfile(output_dir, 'Stability_Diagrams');
folders.diffusion = fullfile(output_dir, 'Diffusive_Errors');
folders.dispersion = fullfile(output_dir, 'Dispersive_Errors');
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
scheme_names = {'Explicit Upwind', 'RK3 Central', 'RK4 Central', ...
                'Implicit Upwind', 'Implicit Central'};
blend_names = {'Lax-Wendroff', 'Beam-Warming', ...
               sprintf('LW-BW blend (b = %.2f)',b_blend)};

%% Build and factorize implicit systems only once
A_upwind = spalloc(N,N,2*N);
A_central = spalloc(N,N,3*N);
for i = 1:N
    left = mod(i-2,N)+1;
    right = mod(i,N)+1;
    A_upwind(i,i) = 1+nu;
    A_upwind(i,left) = -nu;
    A_central(i,i) = 1;
    A_central(i,left) = -nu/2;
    A_central(i,right) = nu/2;
end
solver_upwind = decomposition(A_upwind,'lu');
solver_central = decomposition(A_central,'lu');

%% Solutions and simulation-error evolution
time = (0:N2_steps)'*dt;
for p = 1:2
    % Each output is its own invisible square figure.
    u0 = initial_conditions{p}(:);
    U = repmat(u0,1,8);
    error_L1 = zeros(N2_steps+1,8);
    error_Linf = zeros(N2_steps+1,8);

    for n = 0:N2_steps
        t = n*dt;
        x_shifted = mod(x-a*t-xmin,L)+xmin;
        if p == 1
            u_exact = 0.5*exp(-log(2)*(x_shifted/3).^2);
        else
            u_exact = 0.5*double(abs(x_shifted)<=10);
        end
        u_exact = u_exact(:);

        for s = 1:8
            difference = abs(U(:,s)-u_exact);
            error_L1(n+1,s) = mean(difference);
            error_Linf(n+1,s) = max(difference);
        end

        if n == N1_steps || n == N2_steps
            rev = 1 + (n == N2_steps);
            if p == 1
                pulse_file = 'Acoustic';
            else
                pulse_file = 'Square';
            end

            % Original plot: only the five time/space schemes
            fig_solutions = square_figure();
            ax = axes('Parent',fig_solutions);
            hold(ax,'on');
            plot(ax,x,u_exact,'k--','LineWidth',2);
            colors = lines(5);
            for s = 1:5
                plot(ax,x,U(:,s),'-','Color',colors(s,:), ...
                    'LineWidth',1.3);
            end
            xlim(ax,[-35,35]);
            xlabel(ax,'x'); ylabel(ax,'u');
            title(ax,sprintf('%s - %d revolution(s), t = %.0f',pulse_names{p},rev,t));
            legend(ax,[{'Exact solution'},scheme_names],'Location','best');
            grid(ax,'on'); hold(ax,'off');
            save_figure(fig_solutions,folders.solutions, ...
                sprintf('Solutions_%s_%d_Revolution',pulse_file,rev));
            close(fig_solutions);

            % Separate plot: only LW, BW and their constant-b blend
            fig_blend = square_figure();
            ax = axes('Parent',fig_blend);
            hold(ax,'on');
            plot(ax,x,u_exact,'k--','LineWidth',2);
            blend_colors = lines(3);
            for j = 1:3
                plot(ax,x,U(:,j+5),'-','Color',blend_colors(j,:), ...
                    'LineWidth',1.5);
            end
            xlim(ax,[-35,35]);
            xlabel(ax,'x'); ylabel(ax,'u');
            title(ax,sprintf('%s - LW/BW - %d revolution(s), t = %.0f', ...
                pulse_names{p},rev,t));
            legend(ax,[{'Exact solution'},blend_names],'Location','best');
            grid(ax,'on'); hold(ax,'off');
            save_figure(fig_blend,folders.solutions, ...
                sprintf('Solutions_LW_BW_%s_%d_Revolution',pulse_file,rev));
            close(fig_blend);
        end

        if n == N2_steps
            break;
        end
        U(:,1) = explicit_upwind(U(:,1),nu);
        U(:,2) = RK3_central(U(:,2),dt,a,dx);
        U(:,3) = RK4_central(U(:,3),dt,a,dx);
        U(:,4) = solver_upwind\U(:,4);
        U(:,5) = solver_central\U(:,5);
        U(:,6) = lax_wendroff(U(:,6),nu);
        U(:,7) = beam_warming(U(:,7),nu);
        U(:,8) = lw_bw_blend(U(:,8),nu,b_blend);
    end

    % Save original five-scheme L1 and Linf plots separately
    if p == 1
        pulse_file = 'Acoustic';
    else
        pulse_file = 'Square';
    end

    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_L1(:,1:5),'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Mean absolute error');
    title(ax,[pulse_names{p},' - L1 error evolution']);
    legend(ax,scheme_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_L1_',pulse_file]);
    close(fig_errors);

    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_Linf(:,1:5),'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Maximum absolute error');
    title(ax,[pulse_names{p},' - Linf error evolution']);
    legend(ax,scheme_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_Linf_',pulse_file]);
    close(fig_errors);

    % Additional, independent error figures for LW, BW and blend
    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_L1(:,6:8),'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Mean absolute error');
    title(ax,[pulse_names{p},' - LW/BW L1 error evolution']);
    legend(ax,blend_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_L1_LW_BW_',pulse_file]);
    close(fig_errors);

    fig_errors = square_figure();
    ax = axes('Parent',fig_errors);
    plot(ax,time,error_Linf(:,6:8),'LineWidth',1.5);
    xlabel(ax,'Time'); ylabel(ax,'Maximum absolute error');
    title(ax,[pulse_names{p},' - LW/BW Linf error evolution']);
    legend(ax,blend_names,'Location','best'); grid(ax,'on');
    save_figure(fig_errors,folders.sim_errors,['Error_Linf_LW_BW_',pulse_file]);
    close(fig_errors);

end

%% LW-BW stability, diffusion and dispersion
% Each nu gets one figure per category, using the same b values.
theta = linspace(-pi,pi,2001);
theta_pos = linspace(0,pi,1001);
nu_values = [0.8,1,1.5,2];
b_values = linspace(0,1,5);
legend_b = arrayfun(@(b)sprintf('b = %.2f',b),b_values,'UniformOutput',false);

for k = 1:length(nu_values)
    nu_k = nu_values(k);
    fig_stability = square_figure(); hold on;
    plot(cos(theta),sin(theta),'k--','LineWidth',1.5);
    fig_diffusion = square_figure(); hold on;
    fig_dispersion = square_figure(); hold on;

    for j = 1:length(b_values)
        b = b_values(j);
        G_LW = 1-1i*nu_k*sin(theta)+nu_k^2*(cos(theta)-1);
        G_BW = 1-nu_k/2*(3-4*exp(-1i*theta)+exp(-2i*theta)) ...
                 +nu_k^2/2*(1-2*exp(-1i*theta)+exp(-2i*theta));
        G = (1-b)*G_LW+b*G_BW;
        figure(fig_stability);
        plot(real(G),imag(G),'LineWidth',1.5);

        G_pos = G(theta>=0);
        E_diff = 1-abs(G_pos);
        phase_num = unwrap(angle(G_pos));
        E_disp = phase_num+nu_k*theta_pos;
        E_disp(abs(G_pos)<1e-10) = NaN;
        figure(fig_diffusion);
        plot(theta_pos,E_diff,'LineWidth',1.5);
        figure(fig_dispersion);
        plot(theta_pos,E_disp,'LineWidth',1.5);
    end

    figure(fig_stability);
    axis equal; grid on;
    xlabel('Re(G)'); ylabel('Im(G)');
    title(sprintf('LW-BW Stability - nu = %.2f',nu_k));
    legend([{'|G| = 1'},legend_b],'Location','best');
    hold off;
    save_figure(fig_stability,folders.stability,sprintf('Stability_nu_%0.2f',nu_k));

    figure(fig_diffusion);
    yline(0,'k--','HandleVisibility','off');
    xlabel('\theta = k\Deltax'); ylabel('1 - |G|');
    title(sprintf('LW-BW Diffusive Error - nu = %.2f',nu_k));
    legend(legend_b,'Location','best'); grid on; hold off;
    save_figure(fig_diffusion,folders.diffusion,sprintf('Diffusive_Error_nu_%0.2f',nu_k));

    figure(fig_dispersion);
    yline(0,'k--','HandleVisibility','off');
    xlabel('\theta = k\Deltax'); ylabel('Phase error (rad/step)');
    title(sprintf('LW-BW Dispersive Error - nu = %.2f',nu_k));
    legend(legend_b,'Location','best'); grid on; hold off;
    save_figure(fig_dispersion,folders.dispersion,sprintf('Dispersive_Error_nu_%0.2f',nu_k));
    close(fig_stability);
    close(fig_diffusion);
    close(fig_dispersion);
end
fprintf('All plots saved under: %s\n',output_dir);

%% Local functions (MATLAB requires these at the end of a script)
function u_new = explicit_upwind(u,nu)
    u_new = u-nu*(u-circshift(u,1));
end

function Lu = RK_operator(a,u,dx)
    Lu = -a*(circshift(u,-1)-circshift(u,1))/(2*dx);
end

function u_new = RK4_central(u,dt,a,dx)
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u+dt*k2/2,dx);
    k4 = RK_operator(a,u+dt*k3,dx);
    u_new = u+(dt/6)*(k1+2*k2+2*k3+k4);
end

function u_new = RK3_central(u,dt,a,dx)
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u-dt*k1+2*dt*k2,dx);
    u_new = u+(dt/6)*(k1+4*k2+k3);
end

function fig = square_figure()
    fig = figure('Visible','off','Color','w', ...
        'Units','pixels','Position',[100 100 900 900], ...
        'PaperPositionMode','auto');
end

<<<<<<< HEAD
%% Error evolution over time
fig_solutions = figure;
%fig_errors = figure;  %supressed bc its crushing the other graph (line 225)
=======
function save_figure(fig,folder,basename)
    exportgraphics(fig,fullfile(folder,[basename,'.png']),'Resolution',300);
    savefig(fig,fullfile(folder,[basename,'.fig']));
end
>>>>>>> 0280085 (updated the code for ex1 and added the figures)

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

<<<<<<< HEAD
% Repeat for both initial conditions
for p = 1:2

    u0 = initial_conditions{p};
    u0 = u0(:);

    % Initialize solutions
    U = repmat(u0, 1, 5);

    % Store errors (rows = time, columns = scheme)
    error_L1 = zeros(N2_steps+1, 5);
    error_Linf = zeros(N2_steps+1, 5);

    for n = 0:N2_steps

        % Current time
        t = n * dt;

        % Exact solution with periodic boundaries
        x_shifted = mod(x - a*t + 100, L) - 100;

        if p == 1
            u_exact = 0.5 * exp( ...
                -log(2)*(x_shifted/3).^2);
            u_exact = u_exact(:);
        else
            u_exact = 0.5 * double(abs(x_shifted) <= 10);
            u_exact = u_exact(:);
        end

        % Compute errors for each scheme
        for s = 1:5

            difference = abs(U(:,s) - u_exact);

            error_L1(n+1,s) = mean(difference);
            error_Linf(n+1,s) = max(difference);

        end
        
        if n==N1_steps
            % Plot solutions for 1 revolution
            figure(fig_solutions);
            subplot(2,2,(p-1)*2 + 1);
            hold on;
    
            plot(x, u_exact, 'k--', 'LineWidth', 2);
    
            plot(x, U(:,1), 'b-', 'LineWidth', 1.2);
            plot(x,U(:,2), 'r-', 'LineWidth', 1.2);
            plot(x,U(:,3), 'g-', 'LineWidth', 1.2);
            plot(x,U(:,4), 'm-', 'LineWidth', 1.2);
            plot(x,U(:,5),'c-', 'LineWidth', 1.2);
            xlim([-35 35]);
            xlabel('x');
            ylabel('u');
    
            title(sprintf('%s, t = %.0fs', ...
                pulse_names{p}, n*dt));
    
            legend('Exact solution', ...
                   'Explicit Upwind', ...
                   'RK3 Central', ...
                   'RK4 Central', ...
                   'Implicit Upwind', ...
                   'Implicit Central', ...
                   'Location', 'best');
            grid on;
            hold off;           
        end

        if n==N2_steps
            % Plot solutions for 2 revolutions
            subplot(2,2,(p-1)*2 + 2);
            hold on;
    
            plot(x, u_exact, 'k--', 'LineWidth', 2);
    
            plot(x, U(:,1), 'b-', 'LineWidth', 1.2);
            plot(x,U(:,2), 'r-', 'LineWidth', 1.2);
            plot(x,U(:,3), 'g-', 'LineWidth', 1.2);
            plot(x,U(:,4), 'm-', 'LineWidth', 1.2);
            plot(x,U(:,5),'c-', 'LineWidth', 1.2);
            xlim([-35 35]);
            xlabel('x');
            ylabel('u');
    
            title(sprintf('%s, t = %.0fs', ...
                pulse_names{p}, n*dt));
    
            legend('Exact solution', ...
                   'Explicit Upwind', ...
                   'RK3 Central', ...
                   'RK4 Central', ...
                   'Implicit Upwind', ...
                   'Implicit Central', ...
                   'Location', 'best');
            grid on;
            hold off;
            break;
            
        end

        % Advance all schemes one time step
        U(:,1) = explicit_upwind(U(:,1), nu);

        U(:,2) = RK3_central(U(:,2), dt, a, dx);

        U(:,3) = RK4_central(U(:,3), dt, a, dx);

        U(:,4) = implicit_upwind(U(:,4),nu,N);

        U(:,5) = implicit_central(U(:,5),nu,N);

    end
    
    %figure(fig_errors);  %comment line 99 
    figure;
    subplot(1,2,1);
    plot(time, error_L1, 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Mean absolute error');
    title([pulse_names{p}, ' - L1 error']);
    legend(scheme_names, 'Location', 'best');
    grid on;

    subplot(1,2,2);
    plot(time, error_Linf, 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Maximum absolute error');
    title([pulse_names{p}, ' - L infinity error']);
    legend(scheme_names, 'Location', 'best');
    grid on;

    sgtitle(['Error evolution - ', pulse_names{p}]);
=======
% LW-BW constant-b blend
function u_new = lw_bw_blend(u,nu,b)
    u_new = (1-b)*lax_wendroff(u,nu) + b*beam_warming(u,nu);
>>>>>>> 0280085 (updated the code for ex1 and added the figures)
end
