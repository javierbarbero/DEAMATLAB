function [ out ] = deamalm( X, Y, varargin )
%DEAMALM Data envelopment analysis Malmquist indices
%   Computes data envelopment analysis Malmquist indices
%
%   out = DEAMALM(X, Y, Name, Value) computes data envelopment analysis 
%   Malmquist indices with inputs X and outputs Y. Model properties are 
%   specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
%   - 'names': DMU names.
%   - 'fixbaset': previous year 0 (default), first year 1.
%   - 'geomean': compute geometric mean for technological change. Default
%   is 1.
%
%   Example
%     
%      iomalm = deamalm(X, Y, 'orient', 'io');
%      oomalmfixex = deamalm(X, Y, 'orient', 'oo', 'fixbaset', 1);
%
%   See also DEAOUT, DEA, DEAMALMLUEN
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 25, March, 2016
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    if size(X,3) ~= size(Y,3)
        error('Number of time periods in X and Y must be equal')
    end
    
    % Get number of DMUs (m), inputs (m) and outputs (s)
    [n, m, T] = size(X);
    s = size(Y,2);
    
    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    
    % Xeval, X and Yeval, Y must be equal in this function
    if ~isempty(options.Xeval) && size(options.Xeval) ~= size(X)
        error('Xeval and X must be equal')
    end
    
    if ~isempty(options.Yeval) && size(options.Yeval) ~= size(Y)
        error('Yeval and Y must be equal')
    end
    
    % Check RTS
    if ~strcmp(options.rts, 'crs')
        error('Malmquist index only available for ''crs'' returns to scale')
    end
    
    % Check orientation
    if strcmp(options.orient, 'ddf')
        error('Malmquist index for ''ddf'' not yet implemented')
    end
    
    % Create matrices to store results
    M = nan(n, T - 1);
    MTEC = nan(n, T - 1);
    MTC = nan(n, T - 1);
    if options.geomean
        Eflag = nan(n, (T - 1) * 4);
    else
        Eflag = nan(n, (T - 1) * 3);
    end
    
    % For each time period
    for t=1:T-1
        % Get base period
        if isempty(options.fixbaset) || options.fixbaset == 0
            tb = t;
        elseif options.fixbaset == 1
            tb = 1;
        end
        
        % Compute efficiency at base period
        temp_dea = dea(X(:,:,tb), Y(:,:,tb), varargin{:}, 'secondstep', 0);
        t_eff = temp_dea.eff;
        Eflag(:, (2*t - 1)) = temp_dea.exitflag(:, 1);
        
        % Compute efficiency at time t + 1
        temp_dea = dea(X(:,:,t + 1), Y(:,:,t + 1), varargin{:}, 'secondstep', 0);
        t1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 1) = temp_dea.exitflag(:, 1);
               
        % Evaluate each DMU at t + 1, with the others at base period                         
        temp_dea = dea(X(:,:,tb), Y(:,:,tb), varargin{:},...
                    'Xeval', X(:,:, t + 1),...
                    'Yeval', Y(:,:, t + 1), 'secondstep', 0);

        tevalt1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 2) = temp_dea.exitflag(:, 1);
        
        % If geomean
        if options.geomean
            % Evaluate each DMU at t + 1, with the others at base period                         
            temp_dea = dea(X(:,:,t + 1), Y(:,:,t + 1), varargin{:},...
                    'Xeval', X(:,:, tb),...
                    'Yeval', Y(:,:, tb), 'secondstep', 0);

            t1evalt_eff = temp_dea.eff;   
            Eflag(:, (2*t - 1) + 3) = temp_dea.exitflag(:, 1);
        else 
            t1evalt_eff = NaN;
        end
        
        % Inverse efficiencies if 'oo'
        if strcmp(options.orient, 'oo')
            t_eff = 1 ./ t_eff;
            t1_eff = 1 ./ t1_eff;
            tevalt1_eff = 1 ./ tevalt1_eff;
            t1evalt_eff = 1 ./ t1evalt_eff;
        end
        
        % Technical Efficiency
        MTEC(:, t) = t1_eff ./ t_eff;
        
        % Technological Change
        if options.geomean
            MTC(:, t) = ((tevalt1_eff ./ t1_eff) .* (t_eff ./ t1evalt_eff)).^(1/2);
        else
            MTC(:, t) = tevalt1_eff ./ t1_eff;
        end
        
        % Malmquist index
        M(:, t) = MTEC(:, t) .* MTC(:, t);               
        
    end
    
    % Store Malmquist results in the efficiency structure
    eff.M = M;
    eff.MTEC = MTEC;
    eff.MTC = MTC;
    eff.T = T;
    eff.fixbaset = options.fixbaset;
    
    
    % Extract some results
    neval = NaN;
    lambda = NaN;
    slack.X = NaN;
    slack.Y = NaN;
    Xeff = NaN;
    Yeff = NaN;
        
    % Save results
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'radial-malmquist', 'orient', options.orient, 'rts', options.rts,...
        'lambda', lambda, 'slack', slack, ...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/eff.M/eff.MTEC/eff.MTC' );
    

end

