function sol = solve_flow_fv(p, Minf)
    % This file creates a mesh grid and solves the matrix system for a given mach number.
    % For the transonic case, an iteration scheme is used.
    % It returns a structure that contains the values needed to plot the
    % figures.
    
    % Grid and geometry
    dx = p.xMax / p.Nx;
    dy = p.yMax / p.Ny;
    x_faces = linspace(0, p.xMax, p.Nx + 1);
    y_faces = linspace(0, p.yMax, p.Ny + 1);
    x_centers = x_faces(1:end-1) + dx / 2;
    y_centers = y_faces(1:end-1) + dy / 2;
    [X, Y] = meshgrid(x_centers, y_centers);
    
    % slope of the airfoil profile
    slope = zeros(1, p.Nx);
    for i = 1:p.Nx
        x_pos = x_centers(i);
        if (x_pos >= p.xAirfoil) && (x_pos <= p.xAirfoil + p.c)
            if (x_pos <= p.c/2 + p.xAirfoil) 
                slope(i) = p.tHalf / (p.c/2);
            else
                slope(i) = -p.tHalf / (p.c/2);
            end 
        end
    end
    
    grid.Nx = p.Nx; grid.Ny = p.Ny; grid.dx = dx; grid.dy = dy; grid.slope = slope;
  
    U_inf = Minf * sqrt(p.gamma * p.R * p.Tinf);
    phi_current = zeros(p.Nx * p.Ny, 1);
    
    converged = false;
    iterations = 0;
    residual = 1;
    tic;
    % Iteration for the non-linear problem (transonic)
    while ~converged && iterations < p.maxitNonlinear
        iterations = iterations + 1;
        
        % create matrix A and vector b
        [A, b] = cfm.assemble_flow_matrix(grid, p, Minf, U_inf, phi_current);
        
        % Solve the linear system
        [phi_new, ~, ~, ~] = cfm.solve_sparse_system(A, b, p);
        
        if Minf >= 0.8 && Minf <= 1.2 
            phi_new = p.relaxation*phi_new + (1-p.relaxation)*phi_current;
            disp('done');
        end

        phi_current = phi_new;
        % 6. update matrix A coefficients
        [A_check, b_check] = cfm.assemble_flow_matrix( ...
            grid, p, Minf, U_inf, phi_current);
        % residual
        residual = norm(A_check*phi_current - b_check, inf);
        
        if residual < p.tolNonlinear || Minf < 0.8 || Minf > 1.2
            converged = true;
        end
    end
    
    % Calculate the temperature, pressure, density, and mach for each cell
    phi_matrix = reshape(phi_current, [p.Nx, p.Ny])';
    [v, u] = gradient(phi_matrix, dx, dy);
    q_sq = (U_inf + u).^2 + v.^2;
    
    T_ratio = 1 - ((p.gamma - 1) / 2) * Minf^2 * ((q_sq - U_inf^2) / U_inf^2);
    T_local = T_ratio * p.Tinf;
    
    % save the results in the output structure
    sol.X = X;
    sol.Y = Y;
    sol.phiNorm = phi_matrix / (U_inf * p.c);
    sol.Mlocal = sqrt(q_sq) ./ sqrt(p.gamma * p.R * T_local);
    sol.rhoRatio = T_ratio .^ (1 / (p.gamma - 1));
    sol.solver = p.solver;
    sol.converged = converged;
    sol.iterations = iterations;
    sol.residual = residual;
    sol.computationalTime = toc;
end
