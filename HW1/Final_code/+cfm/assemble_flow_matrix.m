function [A, b] = assemble_flow_matrix(grid, p, Minf, U_inf, phi_prev)
    % assemble_flow_matrix Assembles the sparse system for potential flow.
    
    N_total = grid.Nx * grid.Ny;
    A = spalloc(N_total, N_total, 6 * N_total);
    b = zeros(N_total, 1);
    
    Cy = 1 / grid.dy^2;
    
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            k = i + (j - 1) * grid.Nx;
            b_val = 0;
            
            % -------------------------------------------------------------
            % 1. FRONTEIRA DE ENTRADA (Boundary A: Dirichlet phi = 0)
            % Prioridade máxima na coluna i = 1
            % -------------------------------------------------------------
            if i == 1
                A(k, k) = 1;
                b(k) = 0;
                continue; % Salta diretamente para a próxima célula
            end
            
            % -------------------------------------------------------------
            % 2. ESQUEMA EM X E COMPRESSIBILIDADE
            % -------------------------------------------------------------
            if i > 1 && i < grid.Nx
                dphi_dx = (phi_prev(k + 1) - phi_prev(k - 1)) / (2 * grid.dx);
            else
                dphi_dx = (phi_prev(k) - phi_prev(k - 1)) / grid.dx;
            end
            
            A_local = 1 - Minf^2 - ((p.gamma + 1) * Minf^2 / U_inf) * dphi_dx;
            
            % Seleção do esquema longitudinal
            if strcmpi(p.scheme_x, 'central')
                Cx = A_local / grid.dx^2;
                a_W = Cx; a_E = Cx; a_N = Cy; a_S = Cy;
                a_P = -2 * (Cx + Cy); a_WW = 0;
                
            elseif strcmpi(p.scheme_x, 'backward')
                Cx = -A_local / grid.dx^2; 
                a_W = 2 * Cx; a_E = 0; a_N = Cy; a_S = Cy;
                a_P = -(Cx + 2 * Cy); a_WW = -Cx;
                
            elseif strcmpi(p.scheme_x, 'transonic_switch')
                if A_local >= 0
                    Cx = A_local / grid.dx^2;
                    a_W = Cx; a_E = Cx; a_N = Cy; a_S = Cy;
                    a_P = -2 * (Cx + Cy); a_WW = 0;
                else
                    Cx = -A_local / grid.dx^2; 
                    a_W = 2 * Cx; a_E = 0; a_N = Cy; a_S = Cy;
                    a_P = -(Cx + 2 * Cy); a_WW = -Cx;
                end
            else
                error('Esquema em X não reconhecido!');
            end
            
            % -------------------------------------------------------------
            % 3. CONDIÇÕES DE FRONTEIRA
            % -------------------------------------------------------------
            
            % FUNDO E PERFIL (Boundary D e E: j = 1) -> dphi/dy = U_inf * slope
            if j == 1
                if strcmpi(p.bc_bottom, 'tangent')
                    fluxo = U_inf * grid.slope(i);
                    % Elimina o vizinho fantasma inferior (S):
                    % dphi/dy = (phi_P - phi_S)/dy = fluxo => phi_S = phi_P - dy*fluxo
                    a_P = a_P + a_S;
                    b_val = b_val + a_S * grid.dy * fluxo;
                    a_S = 0;
                end
            end
            
            % TOPO (Boundary B: j = Ny) -> dphi/dy = 0
            if j == grid.Ny
                if strcmpi(p.bc_top, 'neumann')
                    % phi_N = phi_P
                    a_P = a_P + a_N;
                    a_N = 0;
                end
            end
            
            % SAÍDA (Boundary C: i = Nx)
            if i == grid.Nx
                if strcmpi(p.bc_outlet, 'neumann')
                    if A_local >= 0
                        % phi_E = phi_P
                        a_P = a_P + a_E;
                        a_E = 0;
                    end
                elseif strcmpi(p.bc_outlet, 'extrapolation')
                    % Zero tração/condição natural
                    a_E = 0;
                end
            end
            
            % -------------------------------------------------------------
            % 4. MONTAGEM NA MATRIZ ESPARSA
            % -------------------------------------------------------------
            A(k, k) = a_P;
            if a_W ~= 0,         A(k, k - 1)       = a_W;  end
            if a_E ~= 0,         A(k, k + 1)       = a_E;  end
            if a_S ~= 0 && j > 1, A(k, k - grid.Nx) = a_S;  end
            if a_N ~= 0 && j < grid.Ny, A(k, k + grid.Nx) = a_N;  end
            if a_WW ~= 0 && i > 2
                A(k, k - 2) = a_WW;
            end
            
            b(k) = b_val;
            
        end
    end
end