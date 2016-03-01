% Unit-test for basic dea models
clear
clc

% Tolerance
tol = 1e-5;

% Load data
load 'deadataFLS'

% INPUT-ORIENTED CRS
io = dea(X, Y, 'orient', 'io');

% INPUT-ORIENTED VRS
io_vrs = dea(X, Y, 'orient', 'io', 'rts', 'vrs');

% OUTPUT-ORIENTED CRS
oo = dea(X, Y, 'orient', 'oo');

% OUTPUT-ORIENTED VRS
oo_vrs = dea(X, Y, 'orient', 'oo', 'rts', 'vrs');


%% Test io-crs-eff
exp = [1.0000000; 0.6222897; 0.8198562; 1.0000000; 0.3103709; 0.5555556; 1.0000000; 0.7576691; 0.8201058; 0.5000000; 1.0000000];
act = io.eff;
assert(all(abs(act - exp) < tol))

%% Test io-crs-slackX
exp = zeros(11,2);
exp(6, 1) = 4.444444;
exp(9, 1) = 1.640212;
exp(11, 2) = 4;
act = io.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test io-crs-slackY
exp = zeros(11,1);
act = io.slack.Y;
assert(all(all(abs(act - exp) < tol)))


%% Test io-vrs-eff
exp = [1.0000000; 0.8699862; 1.0000000; 1.0000000; 0.7116402; 1.0000000; 1.0000000; 1.0000000; 1.0000000; 0.5000000; 1.0000000];
act = io_vrs.eff;
assert(all(abs(act - exp) < tol))

%% Test io-vrs-slackX
exp = zeros(11,2);
exp(11, 2) = 4;
act = io_vrs.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test io-vrs-slackY
exp = zeros(11,1);
exp(5, 1) = 2.698413;
act = io_vrs.slack.Y;

assert(all(all(abs(act - exp) < tol)))

%% Test oo-crs-eff
exp = [1.000000; 1.606969; 1.219726; 1.000000; 3.221951; 1.800000; 1.000000; 1.319837; 1.219355; 2.000000; 1.000000];
act = oo.eff;
assert(all(abs(act - exp) < tol))

%% Test oo-crs-slackX
exp = zeros(11,2);
exp(6, 1) = 8;
exp(9, 1) = 2;
exp(11, 2) = 4;
act = oo.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test oo-crs-slackY
exp = zeros(11,1);
act = oo.slack.Y;
assert(all(all(abs(act - exp) < tol)))


%% Test oo-vrs-eff
exp = [1.000000; 1.507519; 1.000000; 1.000000; 3.203947; 1.000000; 1.000000; 1.000000; 1.000000; 1.169811; 1.000000];
act = oo_vrs.eff;
assert(all(abs(act - exp) < tol))

%% Test oo-vrs-slackX
exp = zeros(11,2);
exp(10, 1) = 5;
exp(10, 2) = 11;
exp(11, 2) = 4;
act = oo_vrs.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test oo-vrs-slackY
exp = zeros(11,1);
act = oo_vrs.slack.Y;

assert(all(all(abs(act - exp) < tol)))

