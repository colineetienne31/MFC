function phiN = exact_supersonic(X, Y, p, Minf)
    % exact_supersonic Exact linear supersonic potential phi/(Uinf*c) for the diamond airfoil.
    % phi = F(x - kappa*y), kappa = sqrt(Minf^2-1), F'(s) = -Uinf*ys'(s)/kappa.
    % Valid while the Mach waves stay inside the domain (xAirfoil + kappa*yMax > xMax).
    kappa = sqrt(Minf^2 - 1);
    s = X - kappa * Y;
    ys = zeros(size(s));
    front = (s >= p.xAirfoil) & (s <= p.xAirfoil + p.c/2);
    rear  = (s >  p.xAirfoil + p.c/2) & (s <= p.xAirfoil + p.c);
    ys(front) = p.tHalf / (p.c/2) * (s(front) - p.xAirfoil);
    ys(rear)  = p.tHalf - p.tHalf / (p.c/2) * (s(rear) - p.xAirfoil - p.c/2);
    phiN = -ys / (kappa * p.c);
end
