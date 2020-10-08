%% Data Envelopment Analysis Toolbox examples
%
%   Copyright 2019 Inmaculada C. Alvarez, Javier Barbero, Jose L. Zofio
%   http://www.deatoolbox.com
%
%   Version: 1.0
%
% Note: To run this examples the DEA Toolbox should be added to the
% MATLAB search path.
% https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html

% Clear and clc
clear
clc

%% Section 3. Basic DEA models

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04, FARMOUT04];

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% 3.1 Radial input oriented model
io = dea(X, Y, 'orient', 'io', 'names', statenames);
deadisp(io, 'names/eff/slack.X/slack.Y');

io_vrs = dea(X, Y, 'orient', 'io', 'rts', 'vrs', 'names', statenames);
deadisp(io_vrs, 'names/eff/slack.X/slack.Y');

io_scale = deascale(X, Y, 'orient', 'io', 'names', statenames);
deadisp(io_scale);

% 3.2. Radial output oriented model
oo = dea(X, Y, 'orient', 'oo', 'names', statenames);
deadisp(oo, 'names/eff/slack.X/slack.Y');

oo_vrs = dea(X, Y, 'orient', 'oo', 'rts', 'vrs', 'names', statenames);
deadisp(oo_vrs, 'names/eff/slack.X/slack.Y');

oo_scale = deascale(X, Y, 'orient', 'oo', 'names', statenames);
deadisp(oo_scale);

% 3.3 The directional model
ddf = dea(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y, 'names', statenames);
deadisp(ddf, 'names/eff/slack.X/slack.Y');

ddf_vrs = dea(X, Y, 'orient', 'ddf', 'rts', 'vrs', 'Gx', X, 'Gy', Y, 'names', statenames);
deadisp(ddf_vrs, 'names/eff/slack.X/slack.Y');

ddf_scale = deascale(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y, 'names', statenames);
deadisp(ddf_scale);

% 3.4 The additive model
add_vrs = deaaddit(X, Y, 'rts', 'vrs', 'names', statenames);
deadisp(add_vrs, 'names/eff/slack.X/slack.Y');

% Range Adjusted Measure (RAM)
n = size(X, 1);
m = size(X, 2);
s = size(Y, 2);
rhoX = repelem(1 ./ ((m + s) * range(X, 1)), n, 1);
rhoY = repelem(1 ./ ((m + s) * range(Y, 1)), n, 1);
add_ram = deaaddit(X, Y, 'rts', 'vrs', 'rhoX', rhoX, 'rhoY', rhoY, 'names', statenames);
deadisp(add_ram, 'names/eff/slack.X/slack.Y');

% 3.5 Super-efficiency models
super = deasuper(X, Y, 'orient', 'io', 'rts', 'vrs', 'names', statenames);
deadisp(super, 'names/eff/slack.X/slack.Y');

superddf = deasuper(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y, 'rts', 'vrs', 'names', statenames);
deadisp(superddf, 'names/eff/slack.X/slack.Y');

additsuper = deaadditsuper(X, Y, 'rts', 'vrs', 'names', statenames);
deadisp(additsuper, 'names/eff/slack.X/slack.Y');

%% Section 4: Cross efficiency

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04, FARMOUT04];

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% Cross efficiency
cross_ben = deacross(X, Y, 'objective', 'benevolent', 'mean', 'exclusive', 'names', statenames);
deadisp(cross_ben);

cross_agg = deacross(X, Y, 'objective', 'aggressive', 'mean', 'inclusive', 'names', statenames);
deadisp(cross_agg);

%% Section 5: Productivity change

% Data for this section
load 'DataAgriculture'

X00 = [CAPITAL00, LAND00, LABOR00, INTINP00];
Y00 = [LS00, CROPS00, FARMOUT00];

X04 = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y04 = [LS04, CROPS04, FARMOUT04];

X = X00;
X(:,:,2) = X04;

Y = Y00;
Y(:,:,2) = Y04;

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% Malmquist
malmquist = deamalm(X, Y, 'orient', 'io', 'names', statenames);
deadisp(malmquist)

%% Section 6: Allocative models: Overall economic efficiency

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04, FARMOUT04];

X = X ./ 1000;
Y = Y ./ 1000;

W = [WCAPITAL04, WLAND04, WLABOR04, WINTINP04];
P = [PLS04, PCROPS04, PFARMOUT04];

statenames = STATE_NAME04;

% 6.1 Overall cost efficiency
cost = deaalloc(X, Y, 'Xprice', W, 'names', statenames);
deadisp(cost, 'names/eff.T/eff.A/eff.C');

% 6.2 Overall revenue efficiency
revenue = deaalloc(X, Y, 'Yprice', P, 'names', statenames);
deadisp(revenue, 'names/eff.T/eff.A/eff.R');

% 6.3 Overall profit efficiency
profit = deaalloc(X, Y, 'Xprice', W, 'Yprice', P, 'rts', 'vrs', 'names', statenames);
deadisp(profit, 'names/eff.T/eff.A/eff.P');

% 6.4 Profit efficiency: the weighted additive approach
addprofit = deaadditprofit(X, Y, 'rts', 'vrs', 'Xprice', W, 'Yprice', P, 'names', statenames);
deadisp(addprofit, 'names/eff.T/eff.A/eff.P');

%% Section 7: Undesirable outputs

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04,FARMOUT04, ];
Yu = PESTEXPNCASES04;

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% DEA with undesirable outputs
undesirable = deaund(X, Y, Yu, 'names', statenames);
deadisp(undesirable, 'names/eff/slack.Yu');

%% Section 8: Productivity change: The Malmquist-Luenberger

% Data for this section
load 'DataAgriculture'

X00 = [CAPITAL00, LAND00, LABOR00, INTINP00];
Y00 = [LS00, CROPS00, FARMOUT00];
Yu00 = [PESTEXPNCASES00];

X04 = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y04 = [LS04, CROPS04, FARMOUT04];
Yu04 = [PESTEXPNCASES04];

X = X00;
X(:,:,2) = X04;

Y = Y00;
Y(:,:,2) = Y04;

Yu = Yu00;
Yu(:,:,2) = Yu04;

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% Malmquist-Luenberger
ml = deamalmluen(X, Y, Yu, 'names', statenames);
deadisp(ml);

%% Section 9: Bootstrapping DEA estimators

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04, FARMOUT04];

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

rng(1234567, 'twister'); % Set seed for reproducibility

% Start Parallel pool if Parallel Computing Toolbox is installed
try
    parpool;
catch
    disp("Paralell Computing Toolbox not available");    
end

% DEA Bootstrapping
tic
io_b = deaboot(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05, 'names', statenames);
toc
deadisp(io_b);

% Test of RTS
rng(1234567, 'twister'); % Set seed for reproducibility
tic
deatestrts(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05, 'disp', 1);
toc

% Malmquist bootstrap
load 'DataAgriculture'

X00 = [CAPITAL00, LAND00, LABOR00, INTINP00];
Y00 = [LS00, CROPS00, FARMOUT00];

X04 = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y04 = [LS04, CROPS04, FARMOUT04];

X = X00;
X(:,:,2) = X04;

Y = Y00;
Y(:,:,2) = Y04;

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

rng(1234567, 'twister'); % Set seed for reproducibility
tic
malmquist = deamalmboot(X, Y, 'orient', 'io', 'nreps', 200, 'alpha', 0.05, 'names', statenames);
toc
deadisp(malmquist)

% Shut Down Parallel Pool
try
    p = gcp;
    delete(p)
catch
    disp("Paralell Computing Toolbox not available");
end
    
%% Section 10: Advanced options, displaying and exporting results

% Data for this section
load 'DataAgriculture'

X = [CAPITAL04, LAND04, LABOR04, INTINP04];
Y = [LS04, CROPS04, FARMOUT04];

X = X ./ 1000;
Y = Y ./ 1000;

statenames = STATE_NAME04;

% 10.1 Advanced optimization options
opts = optimoptions('linprog','display','off', 'Algorithm','interior-point');
io_ip = dea(X, Y, 'orient', 'io', 'names', statenames, 'secondstep', 0, 'optimopts', opts);

% 10.2 Changing the reference set
Xref = X(2:end, :);
Yref = Y(2:end, :);
Xeval = X(1, :);
Yeval = Y(1, :);
io1 = dea(Xref, Yref, 'orient', 'io', 'Xeval', Xeval, 'Yeval', Yeval);
disp(io1.eff);

% 10.4 Exporting results
io = dea(X, Y, 'orient', 'io', 'names', statenames);
T = dea2table(io);
writetable(T, 'ioresults.csv');

