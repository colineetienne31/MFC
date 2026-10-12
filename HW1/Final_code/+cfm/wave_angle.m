function beta = wave_angle(X, Y, phiN, yRange)
    % wave_angle Mach-wave angle (degrees) measured on the computed field.
    % For each row, the wave front is the first cell where |phi| exceeds 10% of the field maximum;
    % a straight line x = a*y + b is fitted to the fronts and beta = atan(1/a).
    thr = 0.1 * max(abs(phiN(:)));
    ys = [];  xf = [];
    for j = find(Y(:,1) >= yRange(1) & Y(:,1) <= yRange(2))'
        i = find(abs(phiN(j,:)) > thr, 1, 'first');
        if ~isempty(i), ys(end+1) = Y(j,1); xf(end+1) = X(1,i); end %#ok<AGROW>
    end
    c = polyfit(ys, xf, 1);
    beta = atand(1 / c(1));
end
