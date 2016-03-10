% Unit-test for super-efficiency dea model
clear
clc

% Tolerance
tol = 1e-5;

% Load data. Du Liang and Zhu (2010) using Tone (2002) data
names = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};
X = [4 3; 7 3; 8 1; 4 2; 2 4; 10 1; 12 1];
Y = ones(7,1);

% SUPER EFFICIENCY CRS
rhoX = ones(size(X));
rhoY = ones(size(Y));
super_addit = deaadditsuper(X, Y, 'rhoX', rhoX, 'rhoY', rhoY);

%% Test super-efficiency-addit-crs-eff
exp = [0.125; 0.2; 0.5];
act = super_addit.eff(3:5);
assert(all(abs(act - exp) < tol))
