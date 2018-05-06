% Unit-test for allocative efficiency dea model
clear
clc

% Tolerance
tol = 1e-5;

% Data
X = [3 2; 1 3; 4 6];
Y = [3; 5; 6];
W = [4 2];
P = 6;

% Cost model CRS
cost = deaalloc(X, Y, 'Xprice', W);

% Cost model VRS
cost_vrs = deaalloc(X, Y, 'Xprice', W, 'rts', 'vrs');

% Revenue model CRS
revenue = deaalloc(X, Y, 'Yprice', P);

% Revenue model VRS
revenue_vrs = deaalloc(X, Y, 'Yprice', P, 'rts', 'vrs');

%% Test cost-eff-cost
exp = [0.3750000; 1.0000000; 0.4285714];
act = cost.eff.C;
assert(all(abs(act - exp) < tol))

%% Test cost-eff-alloc
exp = [0.4166667; 1.0000000; 0.7142857];
act = cost.eff.A;
assert(all(abs(act - exp) < tol))

%% Test cost-eff-cost-VRS
exp = [0.625; 1; 1];
act = cost_vrs.eff.C;
assert(all(abs(act - exp) < tol))

%% Test cost-eff-alloc-VRS
exp = [0.625; 1; 1];
act = cost_vrs.eff.A;
assert(all(abs(act - exp) < tol))

%% Test revenue-eff-revenue
exp = [0.9; 1; 0.6];
act = revenue.eff.R;
assert(all(abs(act - exp) < tol))

%% Test revenue-eff-alloc
exp = [1; 1; 1];
act = revenue.eff.A;
assert(all(abs(act - exp) < tol))

