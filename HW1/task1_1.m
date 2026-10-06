%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Fluid Mechanics - Homework 1
% Problem 1.1: Plot of Finite-Differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; 
close all; 
clc;
format long


function y = f(x,k)
    y = exp(sin(k*x));
end

function y = f_der(x,k)
    y = k*cos(k*x)*exp(sin(k*x));
end

function y = f_der2(x,k)
    y = k^2 * (cos(k*x)^2 - sin(k*x)) * exp(sin(k*x));
end

function y = f_der3(x,k)
    y = k^3 * (cos(k*x)^3 ...
        - 3*sin(k*x)*cos(k*x) ...
        - cos(k*x)) * exp(sin(k*x));
end

k1 = 0.02;
k2 = 0.2;


%fprintf("Ref value of f'(1.5) = %.15f\n", f_ana(1.5));

steps1 = (2*pi/k1) : -1 : -30;
steps2 = (2*pi/k2) : -1 : -30;
h1 = logspace(log10(2*pi/k1), -20, 100);
h2 = logspace(log10(2*pi/k2), -20, 100);

function y=central_2(k,h,x0)
% 2nd-order central
    y = zeros(size(h));
    for i = 1:length(h)
        h_i = h(i);
        y(i) = (f(x0 + h_i,k) - f(x0 - h_i,k)) / (2*h_i);
    end
end

function y=backward_1(k,h,x0)
% 1st-order backward
    y = zeros(size(h));
    for i = 1:length(h)
        h_i = h(i); 
        y(i) = (f(x0,k)-f(x0-h_i,k)) / h_i;
    end
end

%% FIGURE FOR FUNCTION F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_backward_f1 = abs(backward_1(k1,h1,40) - f_der(40,k1)) / abs(f_der(40,k1));
err_central_f1 = abs(central_2(k1,h1,40) - f_der(40,k1)) / abs(f_der(40,k1));
TE_backward_f1 = h1/2 * abs(f_der2(40,k1));
TE_central_f1  = h1.^2/6 * abs(f_der3(40,k1));

figure;
loglog(h1, err_backward_f1, '-', 'LineWidth', 1.5); hold on;
loglog(h1, err_central_f1, '-', 'LineWidth', 1.5);
loglog(h1, TE_backward_f1, '--', 'LineWidth', 1.5);
loglog(h1, TE_central_f1, '--', 'LineWidth', 1.5);
hold off;

set(gca, 'XDir', 'reverse');

%axis equal;
pbaspect([1 1 1]);

ylim([1e-16 1e0])
xlim([h1(end) h1(1)]);
xlabel('h',FontSize=20);
ylabel('Relative error, ε',FontSize=20);
legend('1st order', '2nd order', 'Location', 'best',fontsize=15);
grid on;
title('Figure 1 - Relative error of FD approximations at x = 40',FontSize=15);


%% FIGURE FOR FUNCTION F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_backward_f2 = abs(backward_1(k2,h2,40) - f_der(40,k2)) / abs(f_der(40,k2));
err_central_f2 = abs(central_2(k2,h2,40) - f_der(40,k2)) / abs(f_der(40,k2));
TE_backward_f2 = h2/2 * abs(f_der2(40,k2));
TE_central_f2  = h2.^2/6 * abs(f_der3(40,k2));

figure;
loglog(h2, err_backward_f2, '-', 'LineWidth', 1.5); hold on;
loglog(h2, err_central_f2, '-', 'LineWidth', 1.5);
loglog(h2, TE_backward_f2, '--', 'LineWidth', 1.5);
loglog(h2, TE_central_f2, '--', 'LineWidth', 1.5);
hold off;

set(gca, 'XDir', 'reverse');

%axis equal;
pbaspect([1 1 1]);

ylim([1e-16 1e0])
xlim([h2(end) h2(1)])
xlabel('h',FontSize=20);
ylabel('Relative error, ε',FontSize=20);
legend('Backward-1 error', ...
       'Central-2 error', ...
       'Backward-1 leading TE', ...
       'Central-2 leading TE', ...
       'Location','best', fontsize=15);
grid on;
title('Function f2 - Relative error of FD approximations at x = 40',FontSize=15);