function [x, dx, u] = solve_advection_1d(L, N, a, u_bc_inlet, S_func, scheme_x)
% Computes 1D steady advection using sparse matrix assembly for a specified finite-difference scheme.

dx = L / N;
x  = linspace(0, L, N+1)';
b  = S_func(x);

if strcmp(scheme_x, 'backward')

    c = a / dx;
    main_diag = c * ones(N+1, 1);
    sub_diag  = -c * ones(N+1, 1);

    % Boundary Condition (Dirichlet at x=0)
    main_diag(1) = 1.0; 
    b(1) = u_bc_inlet;

    % Build lower bidiagonal matrix
    A = spdiags([sub_diag, main_diag], [-1, 0], N+1, N+1);

elseif strcmp(scheme_x, 'central')
    c = a / (2 * dx);
    super_diag = c * ones(N+1, 1);
    sub_diag   = -c * ones(N+1, 1);
    main_diag  = zeros(N+1, 1);

    % Inflow Boundary Condition
    main_diag(1) = 1.0;
    super_diag(1) = 0.0;
    b(1) = u_bc_inlet;

    % Outflow Boundary closure (Upwind at last node to prevent instability)
    main_diag(end) = a / dx;
    sub_diag(end)  = -a / dx;
    super_diag(end) = 0.0;

    % Build tridiagonal matrix
    A = spdiags([sub_diag, main_diag, super_diag], [-1, 0, 1], N+1, N+1);
else
    error('Invalid scheme. Use ''backward'' or ''central''.');
end

% Solve linear system without loops
u = A \ b;
end