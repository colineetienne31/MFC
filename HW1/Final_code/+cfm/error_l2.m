function e = error_l2(d, nb)
    % error_l2 Root-mean-square (L2) norm of the error field d (Ny x Nx).
    % Returns [all cells, interior cells, boundary band]; the band is the nb cells next to each boundary.
    band = false(size(d));
    band([1:nb, end-nb+1:end], :) = true;
    band(:, [1:nb, end-nb+1:end]) = true;
    e = [sqrt(mean(d(:).^2)), sqrt(mean(d(~band).^2)), sqrt(mean(d(band).^2))];
end
