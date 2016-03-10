% Unit-test for additive dea model
clear
clc

% Tolerance
tol = 1e-5;

% Load data
load 'deadataFLS'

% ADDITIVE CRS
rhoX = ones(size(X));
rhoY = ones(size(Y));
add = deaaddit(X, Y, 'rhoX', rhoX, 'rhoY', rhoY);

% ADDITIVE VRS
rhoX = ones(size(X));
rhoY = ones(size(Y));
add_vrs = deaaddit(X, Y, 'rts', 'vrs', 'rhoX', rhoX, 'rhoY', rhoY);

%% Test additive-crs-eff
exp = [0.00000; 10.76923; 10.83784;  0.00000; 22.15385; 17.92308; 0.00000; 12.07692; 11.61379; 34.38462; 4.00000];
act = add.eff;
assert(all(abs(act - exp) < tol))

%% Test additive-vrs-eff
exp = [0.000000; 7.333333; 0.000000; 0.000000; 18.000000; 0.000000; 0.000000; 0.000000; 0.000000; 33.500000; 4.000000];
act = add_vrs.eff;
assert(all(abs(act - exp) < tol))

