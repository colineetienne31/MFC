%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Fluid Mechanics - Homework 1
% Section 1.1: Accuracy of Finite-Difference Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;
format long;

% =========================================================================
% 1. Parameter Block
% =========================================================================
p.x0       = 40;            % Point of evaluation
p.h_min    = 1e-17;         % Minimum step size (h)
p.n_h      = 70;            % Number of step sizes to test
p.y_limits = [1e-14, 1e1];  % Y-axis limits for plots (fixed for comparison)

% =========================================================================
% 2. Test Case Selection (Wave number k)
% =========================================================================
test_case = 'low_k'; % Options: 'low_k', 'high_k', 'all'

if strcmpi(test_case, 'low_k')
    p.k_values = [0.02];
elseif strcmpi(test_case, 'high_k')
    p.k_values = [0.2];
elseif strcmpi(test_case, 'all')
    p.k_values = [0.02, 0.2];
end

% -------------------------------------------------------------------------
% Manual overrides for parametric testing
% (Uncomment any line below to override defaults)
% -------------------------------------------------------------------------
% p.k_values = [0.5, 1.0];  % Test custom high-frequency wave numbers
% p.x0       = 10;          % Evaluate error at a different spatial point
% p.n_h      = 150;         % Increase plot resolution (more h points)

% =========================================================================
% 3. Execution (Calls to +cfm module)
% =========================================================================
results = cfm.s1_accuracy_functions( ...
    p.k_values, p.x0, p.h_min, p.n_h, p.y_limits);