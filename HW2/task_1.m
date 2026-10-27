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

function u_new = RK3_central(u,dt,a,dx)
    k1 = RK_operator(a,u,dx);
    k2 = RK_operator(a,u+dt*k1/2,dx);
    k3 = RK_operator(a,u-dt*k1+2*dt*k2,dx);
    u_new = u + (dt/6)*(k1+4*k2+k3);
end


%% Error evolution over time
fig_solutions = figure;
%fig_errors = figure;  %supressed bc its crushing the other graph (line 225)

time = (0:N2_steps)' * dt;

scheme_names = {'Explicit Upwind', ...
                'RK3 Central', ...
                'RK4 Central', ...
                'Implicit Upwind', ...
                'Implicit Central'};

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
end
