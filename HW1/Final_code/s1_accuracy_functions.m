
function results = s1_accuracy_functions(k_values, x0, h_min, n_h, y_limits)
% Computes and plots finite-difference errors for the functions exp(sin(k*x)).

    results = struct([]);

    for m = 1:numel(k_values)
        k = k_values(m);
        h = logspace(log10(2*pi/k), log10(h_min), n_h);

        % Exact derivative and finite-difference approximations
        exact = f_der(x0, k);
        df_backward = backward_1(k, h, x0);
        df_central = central_2(k, h, x0);

        % Absolute errors
        err_backward = abs(df_backward - exact);
        err_central = abs(df_central - exact);

        % Leading-order truncation errors
        TE_backward = h/2 .* abs(f_der2(x0, k));
        TE_central = h.^2/6 .* abs(f_der3(x0, k));

        % Store results
        results(m).k = k;
        results(m).h = h;
        results(m).err_backward = err_backward;
        results(m).err_central = err_central;
        results(m).TE_backward = TE_backward;
        results(m).TE_central = TE_central;

        % Generate figure
        plot_accuracy(h, err_backward, err_central, ...
                      TE_backward, TE_central, m, y_limits);
    end
end


function y = f(x, k)
% Evaluates f(x) = exp(sin(k*x)).
    y = exp(sin(k*x));
end


function y = f_der(x, k)
% Evaluates the exact first derivative of f.
    y = k*cos(k*x)*exp(sin(k*x));
end


function y = f_der2(x, k)
% Evaluates the exact second derivative of f.
    y = k^2 * (cos(k*x)^2 - sin(k*x)) * exp(sin(k*x));
end


function y = f_der3(x, k)
% Evaluates the exact third derivative of f.
    y = k^3 * (cos(k*x)^3 ...
        - 3*sin(k*x)*cos(k*x) ...
        - cos(k*x)) * exp(sin(k*x));
end


function y = central_2(k, h, x0)
% Computes the second-order central finite-difference approximation.
    y = (f(x0 + h, k) - f(x0 - h, k)) ./ (2*h);
end


function y = backward_1(k, h, x0)
% Computes the first-order backward finite-difference approximation.
    y = (f(x0, k) - f(x0 - h, k)) ./ h;
end


function plot_accuracy(h, err_backward, err_central, ...
                       TE_backward, TE_central, m, y_limits)
% Plots finite-difference errors and leading-order truncation errors.

    figure;

    loglog(h, err_backward, '-', 'LineWidth', 1.5);
    hold on;
    loglog(h, err_central, '-', 'LineWidth', 1.5);
    loglog(h, TE_backward, '--', 'LineWidth', 1.5);
    loglog(h, TE_central, '--', 'LineWidth', 1.5);
    hold off;

    set(gca, 'XDir', 'reverse', 'FontSize', 14.5);
    xlim([min(h), max(h)]);
    ylim(y_limits);
    pbaspect([1 1 1]);

    xlabel('Grid size, h', 'FontSize', 20);
    ylabel('Absolute error, \epsilon', 'FontSize', 20);

    legend('Backward-1 error', ...
           'Central-2 error', ...
           'Backward-1 leading TE', ...
           'Central-2 leading TE', ...
           'Location', 'best', 'FontSize', 15);

    grid on;
    title(sprintf('Function f_%d: FD errors at x = %g', m, 40), ...
          'FontSize', 16);
end
