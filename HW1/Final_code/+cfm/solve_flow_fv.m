function sol = solve_flow_fv(p, Minf)
    % solve_flow_fv Resolves the finite-volume potential flow for a given Mach number.
    
    % 1. Geração da Malha e Geometria
    dx = p.xMax / p.Nx;
    dy = p.yMax / p.Ny;
    x_faces = linspace(0, p.xMax, p.Nx + 1);
    y_faces = linspace(0, p.yMax, p.Ny + 1);
    x_centers = x_faces(1:end-1) + dx / 2;
    y_centers = y_faces(1:end-1) + dy / 2;
    [X, Y] = meshgrid(x_centers, y_centers);
    
    % Declive do perfil (condição de fronteira no fundo)
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
    
    % Guardar dados da grelha numa estrutura para passar às outras funções
    grid.Nx = p.Nx; grid.Ny = p.Ny; grid.dx = dx; grid.dy = dy; grid.slope = slope;
    
    % 2. Parâmetros Base e Inicialização
    U_inf = Minf * sqrt(p.gamma * p.R * p.Tinf);
    phi_current = zeros(p.Nx * p.Ny, 1);
    
    converged = false;
    iterations = 0;
    residual = 1;
    tic;
    % Ciclo iterativo para a não-linearidade (Transónico)
    while ~converged && iterations < p.maxitNonlinear
        iterations = iterations + 1;
        
        % Montar a matriz A e vetor b com base no phi atual (para calcular o dphi/dx no switch)
        [A, b] = cfm.assemble_flow_matrix(grid, p, Minf, U_inf, phi_current);
        
        % Resolver o sistema linear
        [phi_new, ~, ~, ~] = cfm.solve_sparse_system(A, b, p);
        
        % Aplicar relaxação
        phi_new = p.relaxation * phi_new + (1 - p.relaxation) * phi_current;
        
        % Calcular resíduo (variação máxima)
        residual = max(abs(phi_new - phi_current));
        phi_current = phi_new;
        
        if residual < p.tolNonlinear || Minf < 0.8 || Minf > 1.2
            % Se for estritamente subsónico ou supersónico, é linear, logo converge numa iteração
            converged = true;
        end
    end
    
    % 3. Pós-Processamento
    phi_matrix = reshape(phi_current, [p.Nx, p.Ny])';
    [v, u] = gradient(phi_matrix, dx, dy);
    q_sq = (U_inf + u).^2 + v.^2;
    
    T_ratio = 1 - ((p.gamma - 1) / 2) * Minf^2 * ((q_sq - U_inf^2) / U_inf^2);
    T_local = T_ratio * p.Tinf;
    
    % Guardar os resultados na estrutura de saída esperada pelo script principal
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