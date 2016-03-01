% Unit-test for ddf dea models
clear
clc

% Tolerance
tol = 1e-5;

% Load data
load 'deadataFLS'

% DDF Inputs CRS
ddf_i = dea(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', 0);

% DDF Inputs VRS
ddf_i_vrs = dea(X, Y, 'orient', 'ddf', 'rts', 'vrs', 'Gx', X, 'Gy', 0);

% DDF Outputs CRS
ddf_o = dea(X, Y, 'orient', 'ddf', 'Gx', 0, 'Gy', Y);

% DDF Outputs VRS
ddf_o_vrs = dea(X, Y, 'orient', 'ddf', 'rts', 'vrs', 'Gx', 0, 'Gy', Y);

% DDF InputOutput CRS
ddf_io = dea(X, Y, 'orient', 'ddf', 'Gx', X, 'Gy', Y);

% DDF InputOutput VRS
ddf_io_vrs = dea(X, Y, 'orient', 'ddf', 'rts', 'vrs', 'Gx', X, 'Gy', Y);

%% Test ddf-inputs-crs-eff
exp = 1 - [1.0000000; 0.6222897; 0.8198562; 1.0000000; 0.3103709; 0.5555556; 1.0000000; 0.7576691; 0.8201058; 0.5000000; 1.0000000];
act = ddf_i.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-inputs-crs-slackX
exp = zeros(11,2);
exp(6, 1) = 4.444444;
exp(9, 1) = 1.640212;
exp(11, 2) = 4;
act = ddf_i.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-inputs-crs-slackY
exp = zeros(11,1);
act = ddf_i.slack.Y;
assert(all(all(abs(act - exp) < tol)))



%% Test ddf-inputs-vrs-eff
exp = 1 - [1.0000000; 0.8699862; 1.0000000; 1.0000000; 0.7116402; 1.0000000; 1.0000000; 1.0000000; 1.0000000; 0.5000000; 1.0000000];
act = ddf_i_vrs.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-inputs-vrs-slackX
exp = zeros(11,2);
exp(11, 2) = 4;
act = ddf_i_vrs.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-inputs-vrs-slackY
exp = zeros(11,1);
exp(5, 1) = 2.698413;
act = ddf_i_vrs.slack.Y;

assert(all(all(abs(act - exp) < tol)))



%% Test ddf-outputs-crs-eff
exp = [1.000000; 1.606969; 1.219726; 1.000000; 3.221951; 1.800000; 1.000000; 1.319837; 1.219355; 2.000000; 1.000000] - 1;
act = ddf_o.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-outputs-crs-slackX
exp = zeros(11,2);
exp(6, 1) = 8;
exp(9, 1) = 2;
exp(11, 2) = 4;
act = ddf_o.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-outputs-crs-slackY
exp = zeros(11,1);
act = ddf_o.slack.Y;
assert(all(all(abs(act - exp) < tol)))


%% Test ddf-outputs-vrs-eff
exp = [1.000000; 1.507519; 1.000000; 1.000000; 3.203947; 1.000000; 1.000000; 1.000000; 1.000000; 1.169811; 1.000000] - 1;
act = ddf_o_vrs.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-outputs-vrs-slackX
exp = zeros(11,2);
exp(10, 1) = 5;
exp(10, 2) = 11;
exp(11, 2) = 4;
act = ddf_o_vrs.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-outputs-vrs-slackY
exp = zeros(11,1);
act = ddf_o_vrs.slack.Y;

assert(all(all(abs(act - exp) < tol)))



%% Test ddf-InputOutput-crs-eff
exp = [1.000000; 1.232825; 1.098988; 1.000000; 1.526285; 1.285714; 1.000000; 1.137871; 1.098837; 1.333333; 1.000000] - 1;
act = ddf_io.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-InputOutput-crs-slackX
exp = zeros(11,2);
exp(6, 1) = 5.714286;
exp(9, 1) = 1.802326;
exp(11, 2) = 4;
act = ddf_io.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-InputOutput-crs-slackY
exp = zeros(11,1);
act = ddf_io.slack.Y;
assert(all(all(abs(act - exp) < tol)))


%% Test ddf-IputOutput-vrs-eff
exp = [1.000000; 1.107613; 1.000000; 1.000000; 1.288360; 1.000000; 1.000000; 1.000000; 1.000000; 1.162866; 1.000000] - 1;
act = ddf_io_vrs.eff;
assert(all(abs(act - exp) < tol))

%% Test ddf-InputOutput-vrs-slackX
exp = zeros(11,2);
exp(10, 2) = 5.456026;
exp(11, 2) = 4;
act = ddf_io_vrs.slack.X;
assert(all(all(abs(act - exp) < tol)))

%% Test ddf-InputOutput-vrs-slackY
exp = zeros(11,1);
exp(5, 1) = 0.3915344;
act = ddf_io_vrs.slack.Y;

assert(all(all(abs(act - exp) < tol)))
