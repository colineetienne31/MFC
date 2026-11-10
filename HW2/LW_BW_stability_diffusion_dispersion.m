%% LW-BW stability, diffusion and dispersion (standalone)
clear; close all; clc;

set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);

%% Output folders (next to this .m file)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
output_dir = fullfile(script_dir, 'Figures');
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

%% LW-BW stability, diffusion and dispersion
% Same nu and b values in all three analyses.
theta = linspace(-pi, pi, 2001);
theta_pos = linspace(0, pi, 1001);
theta_deg = rad2deg(theta_pos);
nu_values = [0.6, 0.8, 1, 1.2, 1.5,2,2.1];
b_values = [0.4, 0.5, 0.6, 0.7, 0.8];
legend_b = arrayfun(@(b) sprintf('b = %.2f', b), ...
    b_values, 'UniformOutput', false);

for k = 1:numel(nu_values)
    nu_k = nu_values(k);

    % Separate, invisible, square figures
    fig_stability = square_figure();
    ax_stability = axes('Parent', fig_stability);
    hold(ax_stability, 'on');
    plot(ax_stability, cos(theta), sin(theta), 'k--', 'LineWidth', 1.5);

    fig_diffusion = square_figure();
    ax_diffusion = axes('Parent', fig_diffusion);
    hold(ax_diffusion, 'on');

    fig_dispersion = square_figure();
    ax_dispersion = axes('Parent', fig_dispersion);
    hold(ax_dispersion, 'on');

    for j = 1:numel(b_values)
        b = b_values(j);

        % Amplification factors: Lax-Wendroff and Beam-Warming
        G_LW = 1 - 1i*nu_k*sin(theta) + nu_k^2*(cos(theta)-1);
        G_BW = 1 - (nu_k/2)*(3 - 4*exp(-1i*theta) + exp(-2i*theta)) ...
                 + (nu_k^2/2)*(1 - 2*exp(-1i*theta) + exp(-2i*theta));
        G = (1-b)*G_LW + b*G_BW;

        % Polar stability plot
        plot(ax_stability, real(G), imag(G), 'LineWidth', 1.5);

        % Positive wavenumbers (including theta=0)
        G_pos = G(theta >= 0);

        % Relative diffusive error: ideal value = 1
        E_diff = abs(G_pos);
        plot(ax_diffusion, theta_deg, E_diff, 'LineWidth', 1.5);

        % Relative dispersive error: ideal value = 1
        % Phase is undefined at G=0; phase branch is unwrapped.
        phase_num = unwrap(angle(G_pos));
        E_disp = -phase_num ./ (nu_k*theta_pos);
        E_disp(1) = 1; % continuous low-wavenumber limit
        E_disp(abs(G_pos) < 1e-10) = NaN;
        plot(ax_dispersion, theta_deg, E_disp, 'LineWidth', 1.5);
    end

    %% Format and save stability diagram
    axis(ax_stability, 'equal');
    grid(ax_stability, 'on');
    xlabel(ax_stability, 'Re(G)');
    ylabel(ax_stability, 'Im(G)');
    title(ax_stability, sprintf('LW-BW Stability - nu = %.2f', nu_k));
    legend(ax_stability, [{'|G| = 1'}, legend_b], 'Location', 'best');
    hold(ax_stability, 'off');
    save_figure(fig_stability, folders.stability, ...
        sprintf('Stability_nu_%0.2f', nu_k));
    close(fig_stability);

    %% Format and save relative diffusive error
    yline(ax_diffusion, 1, 'k--', 'HandleVisibility', 'off');
    xlabel(ax_diffusion, 'Phase angle \theta (degrees)');
    xlim(ax_diffusion, [0 180]);
    xticks(ax_diffusion, 0:30:180);
    ylabel(ax_diffusion, '\epsilon_D = |G|');
    title(ax_diffusion, sprintf('Relative Diffusive Error - nu = %.2f', nu_k));
    legend(ax_diffusion, legend_b, 'Location', 'best');
    grid(ax_diffusion, 'on');
    hold(ax_diffusion, 'off');
    save_figure(fig_diffusion, folders.diffusion, ...
        sprintf('Diffusive_Error_nu_%0.2f', nu_k));
    close(fig_diffusion);

    %% Format and save relative dispersive error
    yline(ax_dispersion, 1, 'k--', 'HandleVisibility', 'off');
    xlabel(ax_dispersion, 'Phase angle \theta (degrees)');
    xlim(ax_dispersion, [0 180]);
    xticks(ax_dispersion, 0:30:180);
    ylabel(ax_dispersion, '\epsilon_\phi = arg(G)/(-\nu\theta)');
    title(ax_dispersion, sprintf('Relative Dispersive Error - nu = %.2f', nu_k));
    legend(ax_dispersion, legend_b, 'Location', 'best');
    grid(ax_dispersion, 'on');
    hold(ax_dispersion, 'off');
    save_figure(fig_dispersion, folders.dispersion, ...
        sprintf('Dispersive_Error_nu_%0.2f', nu_k));
    close(fig_dispersion);
end

fprintf('LW-BW plots saved under: %s\n', output_dir);

%% Local functions
function fig = square_figure()
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [100 100 900 900], ...
        'PaperPositionMode', 'auto');
end

function save_figure(fig, folder, basename)
    exportgraphics(fig, fullfile(folder, [basename, '.png']), 'Resolution', 300);
    savefig(fig, fullfile(folder, [basename, '.fig']));
end
