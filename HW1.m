% data

T_inf = 288;
p_inf = 101325;
M_inf = 0.3;
c = 1;
t_2 = 0.03 * c;

R = 287;
gama = 1.4;

U_inf = M_inf* sqrt(gama * R * T_inf);
R_inf = p_inf / (R * T_inf);

x_max = 3;
y_max = 2;
x_airfoil = 1;

N_x = 60;
N_y = 40;

%mesh

x_faces = linspace(0,x_max,N_x);
y_faces = linspace(0,y_max,N_y);

dx = x_max/N_x;
dy = y_max/N_y;

x_centers = x_faces(1:end-1) + dx/2;
y_centers = y_faces(1:end-1) + dy/2;

[X,Y] = meshgrid(x_faces,y_faces);

slope = zeros(1,N_x);

for i = 1;N_x
    x_pos = x_centers(i);

    if (x_pos >= x_airfoil) && (x_pos <= x_airfoil + c)
    
        if (x_pos <= c/2 + x_airfoil) 
            slope(i) = t_2/(c/2)
        end

        if (x_pos >= c/2 + x_airfoil)
            slope(i)= -t_2/(c/2)
        end 

    end

end


