% Unit-test for super-efficiency dea model
clear
clc

% Tolerance
tol = 1e-5;

% Load data
load 'deadataFLS'

% SUPER EFFICIENCY CRS
super_io = deasuper(X, Y, 'orient', 'io');

%% Test super-efficiency-io-crs-eff
exp = [1.1302817; 0.6222897; 0.8198562; 1.1168385; 0.3103709; 0.5555556; 1.2449455; 0.7576691; 0.8201058; 0.5000000; 1.0000000];
act = super_io.eff;
assert(all(abs(act - exp) < tol))
