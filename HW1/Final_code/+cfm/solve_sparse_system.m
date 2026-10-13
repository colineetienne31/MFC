function [phi, flag, relres, iter] = solve_sparse_system(A, b, p)
    % Solves the linear system based on the chosen solver algorithm.

    if strcmpi(p.solver, 'backslash')
        phi = A \ b;
        flag = 0; relres = 0; iter = 1;

    elseif strcmpi(p.solver, 'gmres')
        [phi, flag, relres, iter] = gmres(A, b, [], p.tolLinear, p.maxitLinear);

    elseif strcmpi(p.solver, 'bicgstab')
        [phi, flag, relres, iter] = bicgstab(A, b, p.tolLinear, p.maxitLinear);

    elseif strcmpi(p.solver, 'gauss-seidel')
        % Gauss-Seidel using sparse triangular splitting.
        % Tolerances and maximum iterations are defined in run_s2_flow.m.

        phi = zeros(size(b));
        d = diag(A);

        if any(d == 0)
            error('Gauss-Seidel requires a nonzero diagonal.');
        end

        if isfield(p, 'maxitStationary')
            maxit = p.maxitStationary;
        else
            maxit = p.maxitLinear;
        end

        if isfield(p, 'tolStationary')
            tol = p.tolStationary;
        else
            tol = p.tolLinear;
        end

        if maxit < 1 || tol <= 0
            error('Gauss-Seidel maximum iterations and tolerance must be positive.');
        end

        % A = D + L + U
        % Gauss-Seidel uses the lower triangular matrix D + L.
        Lgs = tril(A, -1) + spdiags(d, 0, numel(d), numel(d));

        normb = norm(b);
        if normb == 0
            normb = 1;
        end

        flag = 1;
        relres = norm(b - A*phi) / normb;
        iter = 0;

        while iter < maxit && relres > tol

            % Gauss-Seidel update
            phi = phi + Lgs \ (b - A*phi);

            iter = iter + 1;
            relres = norm(b - A*phi) / normb;

            if ~isfinite(relres)
                flag = 2;
                break;
            end
        end

        if relres <= tol
            flag = 0;
        end

    else
        error('Unknown solver: %s', p.solver);
    end
end
