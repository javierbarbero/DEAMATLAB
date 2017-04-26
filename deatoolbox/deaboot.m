function [ out ] = deaboot(  X, Y, varargin )
%DEABOOT Data envelopment analysis bootstrap.
%   Computes data envelopment analysis bootstrap following Simar and Wilson
%   (1998) and Bogetoft and Otto (2001)
%
%   out = DEABOOT(X, Y, Name, Value) computes data envelopment analysis
%   bootstrap model with inputs X and outputs Y. Model properties are 
%specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
%   directional distane function 'ddf'.
%   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
%   returns to sacle 'vrs'.
%   - 'names': DMU names.
%   - 'nreps': number of bootstrap replications. Default is 200.
%   - 'alpha': alpha value for confidence intervals. Default is 0.05.
%
%   Example
%     
%       io_b = deaboot(X, Y, 'orient', 'io', 'nreps', 200);
%       deadisp(io_b);
%
%   See also DEAOUT, DEA
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 26, April, 2017
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end    
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
    [n, m] = size(X);
    s = size(Y,2);
    
    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    
    % Orientation
    switch(options.orient)
        case {'io','input'}
            orient = 'io';
        case {'oo','output'}
            orient = 'oo';
        case {'ddf'}
            % orient = 'ddf';
            error('DEA Bootstrap not availalbe for ''ddf''');
        case {'none'}
            error('Radial Model is Oriented');
    end   
    
    % RETURNS TO SCALE
    rts = options.rts;
    
    % Xeval, X and Yeval, Y must be equal in this function
    if ~isempty(options.Xeval) && size(options.Xeval) ~= size(X)
        error('Xeval and X must be equal')
    end
    
    if ~isempty(options.Yeval) && size(options.Yeval) ~= size(Y)
        error('Yeval and Y must be equal')
    end   
    
    % Original efficiency estimates
    if isempty(options.effRef)
        efforig = dea(X, Y, varargin{:});
        efforig = efforig.eff;
    else
        % Use reference efficiencies (when computing the RTS test)
        efforig = options.effRef;        
    end

    % Invert efficiencies if 'io'
    if strcmp(orient, 'io')
        efforig_this = 1 ./ efforig;
    else
        efforig_this = efforig;
    end

    % Get only inefficient units
    e = 1e-5;
    effn = efforig(efforig_this > 1 + e);

    % Build the 2m set
    eff2m = [2 - effn; effn];

    % Compute hm (Optimal bandwith. After equation 3.24)

    %hm = 1.06 * min([eff2m_s, eff2m_iqr ./ 1.34]) * (length(eff2m)) ^(-1/5);
    hm = 0.9 * min([std(eff2m),iqr(eff2m) ./ 1.34]) * (length(eff2m)) ^(-1/5);

    % Adjust bandwith (Equation 3.26)
    % h = hm .*(length(eff2m) ./ length(effn)) * (std(efforig_this) ./ eff2m_s )
    h = hm .* ((length(eff2m) ./ length(efforig)))^(1/5) * (std(efforig_this) ./ std(eff2m) ); 

    % Get number of Bootstrap replications and significance
    nreps = options.nreps;
    alph = options.alpha;
    
    Eflag = nan(n, nreps);
    effBoot = nan(n, nreps);
    
    % For each replication
    parfor i=1:nreps
        % Dario and Simar (2007)

        % Invert efficiencies if 'io'
        efforigb = efforig_this;

        % [1]: Build the 2m set and sample n observations
        eff2m = [2 - efforigb; efforigb];
        effn = datasample(eff2m, n);

        % [2]: Perturbate
        effn_star = effn + h .* randn(n, 1);

        % [3]: Refine and correct for the mean and variance of smoothed values
        % eff_star = mean(effn) + (effn_star - mean(effn)) ./ (sqrt(1 + h^2 ./ var(effn)));
        eff_star = mean(effn) + (effn_star - mean(effn)) ./ (sqrt(1 + h^2 ./ var(eff2m)));

        % [4]: Reflect values
        eff_starr = eff_star;
        eff_starr(eff_star < 1) = 2 - eff_star(eff_star < 1);

        % Invert efficiencies if 'io'
        if strcmp(orient, 'io')
            eff_starr = 1 ./ eff_starr;
        end

        % [5] Generate inefficient inputs or output and perform DEA
        Xref = [];
        Yref = [];
        switch(orient)
            case 'io'
                Xref = repmat(efforig ./ eff_starr, 1, m) .* X;
                Yref = Y;
            case 'oo'
                Xref = X;
                Yref = repmat(efforig ./ eff_starr, 1, s) .* Y;
        end
        
        % Perform DEA with the subsample
        b = dea(Xref, Yref, 'orient', orient, varargin{:}, 'Xeval', X, 'Yeval', Y, 'secondstep', 0);
                
        % Get efficiency for the bootstrap sample
        effBoot(:, i) = b.eff;
        
        Eflag(:, i ) = b.exitflag(:, 1);
        
    end

    % Invert efficiencies if 'io'
    if strcmp(orient, 'io')
        effBoot = 1 ./ effBoot;
    end

    % Bootstrap efficiency
    eff.bias = mean(effBoot, 2) - efforig_this;
    eff.b = efforig_this - eff.bias;

    % Confidence interval
    eff.c = repelem(efforig_this, 1, 2) + quantile(repmat(efforig_this, 1, nreps) - effBoot, [0.5*alph, 1 - 0.5*alph], 2);

    % Store other results
    eff.o = efforig;
    eff.Boot = effBoot;
    
    % Invert efficiencies if 'io'
    if strcmp(orient, 'io')
        eff.Boot = 1 ./ effBoot;
        eff.b = 1 ./ eff.b;
        eff.bias = efforig - eff.b;
        eff.c = 1./ eff.c;
        eff.c = eff.c(:,[2,1]); % Change column order if 'io'
    end

    % Coompute variance
    eff.var = var(effBoot, 0, 2);
    
    % Extract some results
    neval = NaN;
    lambda = NaN;
    slack.X = NaN;
    slack.Y = NaN;
    Xeff = NaN;
    Yeff = NaN;
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'radial-bootstrap', 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/eff.o/eff.b/eff.c',...
        'nreps', nreps, 'alpha', alph);
    
    
end

