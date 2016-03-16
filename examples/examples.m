%% Data Envelopment Analysis Toolbox examples
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 16, March, 2016
%

% Clear and clc
clear
clc

%% Section 3. Basic DEA models

load 'deadataFLS'

% 3.1 Radial input oriented model
io = dea(X, Y, 'orient', 'io');
deadisp(io);

io_vrs = dea(X, Y, 'orient', 'io', 'rts', 'vrs');
deadisp(io_vrs);

io_scale = deascale(X, Y, 'orient', 'io');
deadisp(io_scale);

% 3.2. Radial output oriented model
oo = dea(X, Y, 'orient', 'oo');
deadisp(oo);

oo_vrs = dea(X, Y, 'orient', 'oo', 'rts', 'vrs');
deadisp(oo_vrs);

oo_scale = deascale(X, Y, 'orient', 'oo');
deadisp(oo_scale);

% 3.3 The directional model
ddf = dea(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y);
deadisp(ddf);

ddf_vrs = dea(X, Y, 'orient', 'ddf', 'rts', 'vrs', 'Gx', X, 'Gy', Y);
deadisp(ddf_vrs);

ddf_scale = deascale(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y);
deadisp(ddf_scale);

% 3.4 The additive model
add_vrs = deaaddit(X, Y, 'rts', 'vrs');
deadisp(add_vrs);

% Range Adjusted Measure (RAM)
n = size(X, 1);
m = size(X, 2);
s = size(Y, 2);
rhoX = repelem(1 ./ ((m + s) * range(X, 1)), n, 1);
rhoY = repelem(1 ./ ((m + s) * range(Y, 1)), n, 1);
add_ram = deaaddit(X, Y, 'rts', 'vrs', 'rhoX', rhoX, 'rhoY', rhoY);
deadisp(add_ram);

% 3.5 Super-efficiency models
super = deasuper(X, Y, 'orient', 'io');
deadisp(super);

superddf = deasuper(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y);
deadisp(superddf);

additsuper = deaadditsuper(X, Y);
deadisp(additsuper);

%% Section 4: Productivity change

% Data for this section
X = [2; 3; 5; 4; 4];
X(:, :, 2) = [1; 2; 4; 3; 4];
X(:, :, 3) = [0.5; 1.5; 3; 2; 4];

Y = [1; 4; 6; 3; 5];
Y(:, :, 2) = [1; 4; 6; 3; 3];
Y(:, :, 3) = [  2;   4; 6; 3; 1];

malmquist = deamalm(X, Y, 'orient', 'io');
deadisp(malmquist)

malmquist1 = deamalm(X, Y, 'orient', 'oo', 'fixbaset', 1);
deadisp(malmquist1)

%% Section 5: Allocative models: Overall economic efficiency

% Data for this section
X = [3 2; 1 3; 4 6];
Y = [3; 5; 6];
W = [4 2];
P = 6;

% 5.1 Overall cost efficiency
cost = deaalloc(X, Y, 'Xprice', W);
deadisp(cost);

% 5.2 Overall revenue efficiency
revenue = deaalloc(X, Y, 'Yprice', P);
deadisp(revenue);

% 5.3 Overall profit efficiency
profit = deaalloc(X, Y, 'Xprice', W, 'Yprice', P);
deadisp(profit, 'names/X/Y/eff.T/eff.A/eff.P');

% 5.4 Profit efficiency: the weighted additive approach
addprofit = deaadditprofit(X, Y, 'rts', 'vrs', 'Xprice', W, 'Yprice', P);
deadisp(addprofit, 'names/X/Y/eff.T/eff.A/eff.P');

%% Section 6: Undesirable outputs

% Data for this section
X0 = ones(5,1);
Y0 = [7; 5; 1; 3; 4];
Yu0 = [2; 5; 3; 3; 2];

undesirable = deaund(X0, Y0, Yu0);
deadisp(undesirable);

%% Section 7: Productivity change: The Malmquist-Luenberger
X1 = ones(5,1);
Y1 = [8; 5.5; 2; 2; 4];
Yu1 = [1; 3; 2; 4; 1];

X = X0;
X(:, :, 2) = X1;
Y = Y0;
Y(:, :, 2) = Y1;
Yu = Yu0;
Yu(:, :, 2) = Yu1;

ml = deamalmluen(X, Y, Yu);
deadisp(ml);

%% Section 8: Bootstrapping DEA estimators

load 'deadataFLS'
rng(1234567); % Set seed for reproducibility

io_b = deaboot(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05);
deadisp(io_b);

% Test of RTS
rng(1234567); % Set seed for reproducibility
deatestrts(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05, 'disp', 1);

% Malmquist bootstrap
X = [2; 3; 5; 4; 4];
X(:, :, 2) = [1; 2; 4; 3; 4];

Y = [1; 4; 6; 3; 5];
Y(:, :, 2) = [1; 4; 6; 3; 3];

rng(1234567); % Set seed for reproducibility
malmquist = deamalmboot(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05);
deadisp(malmquist)

%% Section 9: Advanced options, displaying and exporting results

% 7.1 Specifying DMU names
X0 = ones(5,1);
Y0 = [7; 5; 1; 3; 4];
Yu0 = [2; 5; 3; 3; 2];
names = {'A','B','C','D','E'}';

undesirable = deaund(X0, Y0, Yu0, 'names', names);
deadisp(undesirable);

% 7.2 Changing the reference set
load 'deadataFLS'
Xref = X(2:end, :);
Yref = Y(2:end, :);
Xeval = X(1, :);
Yeval = Y(1, :);
io1 = dea(Xref, Yref, 'orient', 'io', 'Xeval', Xeval, 'Yeval', Yeval);
disp(io1.eff);

% 7.3 Custom display
load 'deadataFLS'
io = dea(X, Y, 'orient', 'io');
deadisp(io, 'names/eff/Xeff/Yeff');

% 7.4 Exporting results
load 'deadataFLS'
io = dea(X, Y, 'orient', 'io');
T = dea2table(io);
writetable(T, 'ioresults.csv');

