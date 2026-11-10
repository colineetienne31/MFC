%% Computational cost and accuracy: RK3 versus RK4

clear; clc;

L = 200; xmin = -100; a = 1;
nu = 0.8;
dx_values = [2, 1, 0.5]; % Set to 1 for just the homework grid
n_revolutions = 2;
num_repeats_note = 'timeit automatically chooses repeated measurements';

pulse_names = {'Acoustic', 'Square'};
results = table();

for d = 1:numel(dx_values)
    dx = dx_values(d);
    N = round(L/dx);
    assert(abs(N*dx-L) < 1e-10, 'L/dx must be an integer');
    dt = nu*dx/a;
    nsteps = round(n_revolutions*L/(a*dt));
    % Exact final time, rather than assuming rounded nsteps is precisely 2 revolutions
    t_final = nsteps*dt;
    x = xmin + ((0:N-1)' + 0.5)*dx;

    for p = 1:numel(pulse_names)
        if p == 1
            u0 = 0.5*exp(-log(2)*(x/3).^2);
        else
            u0 = 0.5*double(abs(x)<=10);
        end

        % Function handles for benchmarking
        run3 = @() advance_RK(u0, nsteps, dt, a, dx, 3);
        run4 = @() advance_RK(u0, nsteps, dt, a, dx, 4);

        % Warm-up calls (MATLAB JIT compilation/caches)
        run3(); run4();

        t3 = timeit(run3);
        t4 = timeit(run4);

        u3 = run3();
        u4 = run4();

        % Periodically shift the exact center-sampled initial profile
        xshift = mod(x-a*t_final-xmin,L)+xmin;
        if p == 1
            exact = 0.5*exp(-log(2)*(xshift/3).^2);
        else
            exact = 0.5*double(abs(xshift)<=10);
        end

        e3_L1 = mean(abs(u3-exact));
        e4_L1 = mean(abs(u4-exact));
        e3_Linf = max(abs(u3-exact));
        e4_Linf = max(abs(u4-exact));

        next = table(dx, N, nsteps, string(pulse_names{p}), t3, t4, ...
            t4/t3, 100*(t4/t3-1), e3_L1, e4_L1, e3_Linf, e4_Linf, ...
            'VariableNames', {'dx','Ncells','Nsteps','Pulse', ...
            'RK3_seconds','RK4_seconds','TimeRatio_RK4_RK3', ...
            'ExtraTime_percent','RK3_L1','RK4_L1','RK3_Linf','RK4_Linf'});
        results = [results; next]; %#ok<AGROW>
    end
end

disp(results);
fprintf('\nRK3 uses 3 spatial operator evaluations/step; RK4 uses 4.\n');
fprintf('The theoretical operator-work ratio is 4/3 = %.3f.\n',4/3);
fprintf('Measured timings depend on hardware, MATLAB, and problem size.\n');

% Save tabular output without opening a plot
output_folder = fullfile(fileparts(mfilename('fullpath')), 'Figures', 'Computational_Cost');
if ~exist(output_folder,'dir'); mkdir(output_folder); end
writetable(results, fullfile(output_folder,'RK3_RK4_comparison.csv'));
fprintf('Saved comparison table to: %s\n',output_folder);

%% Local helper for timing only
function u = advance_RK(u0,nsteps,dt,a,dx,order)
    u = u0;
    switch order
        case 3
            for n = 1:nsteps
                u = RK3_central(u,dt,a,dx);
            end
        case 4
            for n = 1:nsteps
                u = RK4_central(u,dt,a,dx);
            end
        otherwise
            error('order must be 3 or 4');
    end
end


%% FUNCTIONS 

function Lu = RK_operator(a,u,dx)
% Second-order centered spatial derivative with periodic boundaries
    Lu = -a*(circshift(u,-1)-circshift(u,1))/(2*dx);
end

function u_new = RK3_central(u,dt,a,dx)
% Classical three-stage third-order Runge-Kutta + centered spatial derivative
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u-dt*k1+2*dt*k2,dx);
    u_new = u+(dt/6)*(k1+4*k2+k3);
end

function u_new = RK4_central(u,dt,a,dx)
% Classical four-stage fourth-order Runge-Kutta + centered spatial derivative
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u+dt*k2/2,dx);
    k4 = RK_operator(a,u+dt*k3,dx);
    u_new = u+(dt/6)*(k1+2*k2+2*k3+k4);
end
