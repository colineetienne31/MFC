function [phi, flag, relres, iter] = solve_sparse_system(A, b, p)
    % solve_sparse_system Solves the linear system based on the chosen solver algorithm.
    
    if strcmpi(p.solver, 'backslash')
        phi = A \ b;
        flag = 0; relres = 0; iter = 1;
    elseif strcmpi(p.solver, 'gmres')
        [phi, flag, relres, iter] = gmres(A, b, [], p.tolLinear, p.maxitLinear);
    elseif strcmpi(p.solver, 'bicgstab')
        [phi, flag, relres, iter] = bicgstab(A, b, p.tolLinear, p.maxitLinear);
    else
        % Fallback to the direct solver
        phi = A \ b;
        flag = 0; relres = 0; iter = 1;
    end
end