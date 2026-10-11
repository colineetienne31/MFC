clear; clc; close all;

%% ==========================================================
% REGIME
%% ==========================================================

regime = 'transonic';      % 'subsonic' | 'transonic' | 'supersonic'

switch regime
    case 'subsonic'
        Minf = 0.3;

    case 'transonic'
        Minf = 0.91;

    case 'supersonic'
        Minf = 2.2;
end

%% ==========================================================
% DATA
%% ==========================================================

T_inf = 288;
p_inf = 101325;

c = 1;
t_2 = 0.03*c;

R = 287;
gamma = 1.4;

U_inf = Minf*sqrt(gamma*R*T_inf);
rho_inf = p_inf/(R*T_inf);

x_max = 3;
y_max = 2;

x_airfoil = 1;

Nx = 60;
Ny = 40;

%% ==========================================================
% MESH
%% ==========================================================

x_faces = linspace(0,x_max,Nx+1);
y_faces = linspace(0,y_max,Ny+1);

dx = x_max/Nx;
dy = y_max/Ny;

x_centers = x_faces(1:end-1) + dx/2;
y_centers = y_faces(1:end-1) + dy/2;

[X,Y] = meshgrid(x_centers,y_centers);

%% ==========================================================
% DIAMOND AIRFOIL SLOPE
%% ==========================================================

slope = zeros(1,Nx);

for i = 1:Nx

    x_pos = x_centers(i);

    if (x_pos >= x_airfoil) && ...
       (x_pos <= x_airfoil + c)

        if x_pos <= x_airfoil + c/2
            slope(i) = t_2/(c/2);
        else
            slope(i) = -t_2/(c/2);
        end

    end
end

%% ==========================================================
% LINEAR CASES
%% ==========================================================

if ~strcmp(regime,'transonic')

    N_total = Nx*Ny;

    A = spalloc(N_total,N_total,7*N_total);
    b = zeros(N_total,1);

    for j = 1:Ny

        for i = 1:Nx

            k = i + (j-1)*Nx;

            if strcmp(regime,'subsonic')

                Cx = (1-Minf^2)/dx^2;
                Cy = 1/dy^2;

                aW = Cx;
                aE = Cx;
                aS = Cy;
                aN = Cy;
                aP = -2*(Cx+Cy);

                aWW = 0;

            else

                % Supersonic
                Cx = (Minf^2-1)/dx^2;
                Cy = 1/dy^2;

                aWW = Cx;
                aW  = -2*Cx;
                aE  = 0;

                aS = Cy;
                aN = Cy;

                aP = Cx - 2*Cy;

            end

            b_val = 0;

            %% LEFT BOUNDARY

            if i == 1

                aP = 1;
                aW = 0;
                aE = 0;
                aS = 0;
                aN = 0;
                aWW = 0;

                b_val = 0;

            end

            %% SECOND COLUMN FOR SUPERSONIC

            if strcmp(regime,'supersonic')
                if i == 2

                    aWW = 0;

                end
            end

            %% RIGHT BOUNDARY

            if i == Nx

                if ~strcmp(regime,'supersonic')

                    aP = aP + aE;
                    aE = 0;

                end
            end

            %% TOP BOUNDARY

            if j == Ny

                aP = aP + aN;
                aN = 0;

            end

            %% BOTTOM BOUNDARY

            if j == 1

                airfoil_flow = U_inf*slope(i);

                aP = aP + aS;

                b_val = b_val + aS*dy*airfoil_flow;

                aS = 0;

            end

            %% ASSEMBLY

            A(k,k) = aP;

            if aW~=0 && i>1
                A(k,k-1) = aW;
            end

            if aE~=0 && i<Nx
                A(k,k+1) = aE;
            end

            if aS~=0 && j>1
                A(k,k-Nx) = aS;
            end

            if aN~=0 && j<Ny
                A(k,k+Nx) = aN;
            end

            if aWW~=0 && i>2
                A(k,k-2) = aWW;
            end

            b(k) = b_val;

        end
    end

    phi_vector = A\b;
    phi = reshape(phi_vector,[Nx Ny])';

%% ==========================================================
% TRANSONIC
%% ==========================================================

else

    phi = zeros(Ny,Nx);

    omega = 0.4;

    tol = 1e-7;

    maxIter = 500;

    for iter = 1:maxIter

        phi_old = phi;

        N_total = Nx*Ny;

        A = spalloc(N_total,N_total,7*N_total);
        b = zeros(N_total,1);

        [dphidx, ~] = gradient(phi, dx, dy);

        Acoef = 1 ...
            - Minf^2 ...
            - (gamma+1)*Minf^2.*dphidx/U_inf;

        for j = 1:Ny

            for i = 1:Nx

                k = i + (j-1)*Nx;

                Cy = 1/dy^2;

                aWW = 0;

                if Acoef(j,i) > 0

                    Cx = Acoef(j,i)/dx^2;

                    aW = Cx;
                    aE = Cx;
                    aP = -2*Cx;

                else

                    Cx = abs(Acoef(j,i))/dx^2;

                    aWW = Cx;
                    aW  = -2*Cx;
                    aE  = 0;
                    aP  = Cx;

                end

                aS = Cy;
                aN = Cy;

                aP = aP - 2*Cy;

                b_val = 0;

                %% LEFT

                if i == 1

                    aP = 1;
                    aW = 0;
                    aE = 0;
                    aN = 0;
                    aS = 0;
                    aWW = 0;

                end

                %% TOP

                if j == Ny

                    aP = aP + aN;
                    aN = 0;

                end

                %% RIGHT

                if i == Nx

                    aP = aP + aE;
                    aE = 0;

                end

                %% BOTTOM

                if j == 1

                    airfoil_flow = U_inf*slope(i);

                    aP = aP + aS;

                    b_val = b_val + aS*dy*airfoil_flow;

                    aS = 0;

                end

                %% ASSEMBLY

                A(k,k) = aP;

                if aW~=0 && i>1
                    A(k,k-1)=aW;
                end

                if aE~=0 && i<Nx
                    A(k,k+1)=aE;
                end

                if aS~=0 && j>1
                    A(k,k-Nx)=aS;
                end

                if aN~=0 && j<Ny
                    A(k,k+Nx)=aN;
                end

                if aWW~=0 && i>2
                    A(k,k-2)=aWW;
                end

                b(k)=b_val;

            end
        end

        phi_new = reshape(A\b,[Nx Ny])';

        phi = omega*phi_new + (1-omega)*phi_old;

        err = max(abs(phi(:)-phi_old(:)));

        fprintf('Iter %4d   Error = %.3e\n',iter,err);

        if err < tol
            fprintf('Converged\n');
            break
        end

    end

end

%% ==========================================================
% POST-PROCESSING
%% ==========================================================

[u,v] = gradient(phi,dy,dx);

q_sq = (U_inf + u).^2 + v.^2;

T_ratio = 1 ...
    - ((gamma-1)/2) ...
    * Minf^2 ...
    * ((q_sq-U_inf^2)/U_inf^2);

T_local = T_ratio*T_inf;

rho_ratio = T_ratio.^(1/(gamma-1));

M_local = sqrt(q_sq)./sqrt(gamma*R.*T_local);

phi_norm = phi/(U_inf*c);

%% ==========================================================
% FIGURES
%% ==========================================================

figure('Position',[100 100 900 900])

subplot(3,1,1)

contourf(X,Y,phi_norm,30,'LineColor','none')
colorbar

hold on
plot([x_airfoil x_airfoil+c],[0 0], ...
    'k','LineWidth',3)

title(['\phi/(U_\infty c)   -   ' regime])

axis equal
axis([0 x_max 0 y_max])

subplot(3,1,2)

contourf(X,Y,M_local,30,'LineColor','none')
colorbar

hold on
plot([x_airfoil x_airfoil+c],[0 0], ...
    'k','LineWidth',3)

title('Mach Number')

axis equal
axis([0 x_max 0 y_max])

subplot(3,1,3)

contourf(X,Y,rho_ratio,30,'LineColor','none')
colorbar

hold on
plot([x_airfoil x_airfoil+c],[0 0], ...
    'k','LineWidth',3)

title('\rho / \rho_\infty')

axis equal
axis([0 x_max 0 y_max])

colormap(jet)

%% ==========================================================
% SUPSERSONIC MACH ANGLE
%% ==========================================================

if strcmp(regime,'supersonic')

    beta = asind(1/Minf);

    fprintf('\n');
    fprintf('Mach angle = %.2f degrees\n',beta);

end