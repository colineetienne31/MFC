clear; close all; clc;

%% parameters
L = 200; % domain length
xmin = -100;
xmax = 100;
dx = 1;
a = 1;
nu = 0.8;
N = L/dx;          % Number of cells
dt = nu*dx/a;
T1 = L/a; % time for 1 revolution
T2 = (L/a)*2; % time for 2 revolutions
N1_steps = T1/dt; % 1st revolution
N2_steps = T2/dt; % 2nd revolution
func_type = 'acoustic'; % options: square or acoustic

%% grid and initial conditions
% Cell center coordinates
x = xmin + ((1:N) - 0.5)*dx;

% Acoustic pulse
u_acoustic = 0.5 * exp(-log(2) * (x/3).^2);

% Square pulse
u_square = 0.5 * double(abs(x) <= 10);

initial_conditions = {u_acoustic, u_square};
pulse_names = {'Acoustic pulse', 'Square pulse'};


%% scheme functions %%%%%%%%%%%%%%

function u_new = explicit_upwind(u,nu)
    u_left = circshift(u,1);
    u_new = u-nu*(u-u_left);
end

function u_new = implicit_upwind(u,nu,N)
    A = sparse(N,N);
    for i=1:N
        i_left=i-1;
        if i_left < 1
            i_left = 200;
        end
        A(i,i)=1+nu;
        A(i,i_left)=-nu;
    end
    A_factorized = decomposition(A, 'lu');
    u = u(:); % transpose u
    u_new = A_factorized\u;
end

function u_new = implicit_central(u,nu,N)
    A = sparse(N,N);
    for i=1:N
        i_right=i+1;
        i_left=i-1;
        if i_left < 1
            i_left = 200;
        end
        if i_right > 200
            i_right = 1;
        end
        A(i,i)=1;
        A(i,i_right)=nu/2;
        A(i,i_left)=-nu/2;
    end
    A_factorized = decomposition(A, 'lu');
    u = u(:); % transpose u
    u_new = A_factorized\u;
end

% spatial operator L
function L = RK_operator(a,u,dx)
    u_left = circshift(u,1);
    u_right = circshift(u,-1);
    L = -a*(u_right-u_left)/(2*dx);
end

function u_new = RK4_central(u,dt,a,dx)
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u+dt*k2/2,dx);
    k4 = RK_operator(a,u+dt*k3,dx);
    u_new = u + (dt/6)*(k1+2*k2+2*k3+k4);
end

function u_new = RK3_cental(u,dt,a,dx)
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u-dt*k1+2*dt*k2,dx);
    u_new = u + (dt/6)*(k1+4*k2+k3);
end


%% 4. Run simulations

figure;

for p = 1:2

    u0 = initial_conditions{p};

    for rev = 1:2

        if rev == 1
            nsteps = N1_steps;
        else
            nsteps = N2_steps;
        end

        % Initialize each numerical solution
        u_exp_upwind = u0;
        u_rk4_central = u0;
        u_imp_upwind = u0;
        u_imp_central = u0;
        u_rk3_central = u0;

        % Time integration
        for n = 1:nsteps

            % Explicit Euler + Upwind
            u_exp_upwind = explicit_upwind(u_exp_upwind, nu);

            % RK4 + Central differences
            u_rk4_central = RK4_central(u_rk4_central, dt, a, dx);

            % Implicit Euler + Upwind
            u_imp_upwind = implicit_upwind(u_imp_upwind,nu,N);

            % Implicit Euler + Central
            u_imp_central = implicit_central(u_imp_central,nu,N);

            % RK3 + central
            u_rk3_central = RK3_cental(u_rk3_central,dt,a,dx);

        end

        % Plot solutions
        subplot(2,2,(p-1)*2 + rev);
        hold on;

        plot(x, u0, 'k--', 'LineWidth', 2);

        plot(x, u_exp_upwind, 'b-', 'LineWidth', 1.2);
        plot(x, u_rk4_central, 'r-', 'LineWidth', 1.2);
        plot(x, u_imp_upwind, 'g-', 'LineWidth', 1.2);
        plot(x, u_imp_central, 'm-', 'LineWidth', 1.2);
        plot(x, u_rk3_central, 'c-', 'LineWidth', 1.2);

        xlim([-35 35]);
        xlabel('x');
        ylabel('u');

        title(sprintf('%s, t = %.0f', ...
            pulse_names{p}, nsteps*dt));

        legend('Exact solution', ...
               'Explicit Upwind', ...
               'RK4 Central', ...
               'Implicit Upwind', ...
               'Implicit Central', ...
               'RK3 Central', ...
               'Location', 'best');

        grid on;
        hold off;

    end
end