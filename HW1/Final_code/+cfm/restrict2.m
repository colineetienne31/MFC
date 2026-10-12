function F = restrict2(f)
    % restrict2 Averages 2x2 blocks of fine cells onto the grid with half the cells per direction.
    F = 0.25 * (f(1:2:end, 1:2:end) + f(2:2:end, 1:2:end) ...
              + f(1:2:end, 2:2:end) + f(2:2:end, 2:2:end));
end
