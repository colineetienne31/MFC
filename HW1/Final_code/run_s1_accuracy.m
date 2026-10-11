%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Fluid Mechanics - Homework 1
% Section 1.1: Accuracy of Finite-Difference Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;
format long;

%% Parameters

k_values = [0.02, 0.2];
x0 = 40;

h_min = 1e-17;
n_h = 70;

% Use the same vertical scale for both figures
y_limits = [1e-14, 1e1];

%% Compute errors and generate figures

results = s1_accuracy_functions( ...
    k_values, x0, h_min, n_h, y_limits);
