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
%   - 'fixbaset': base year is previous year 0 (default), or always the 
%     first year 1.
%   - 'period': compute geometric mean of base and comparison periods for 
%     technological change ('geomean'), use base period as reference ('base'),
%     or use comparison period as reference ('comparison').
%     
%
%   Deprecated parameters:
%   - 'geomean': compute geometric mean for technological change. Default
%     is 1. 'geomean' parameter has been deprecated and will dissapear in a
%     future realse. Set the new 'period' parapeter to 'geomean' for the 
%     previous behavior of 'geomean' = 1. Set 'period' to 'base' for the 
%     preivous behaviour of 'geomean' = 0.
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
%   LAST UPDATE: 6, May, 2017
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    if size(X,3) ~= size(Y,3)
        error('Number of time periods in X and Y must be equal')
    end
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
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
    
    % Check orientation
    if strcmp(options.orient, 'ddf')
        error('Malmquist index for ''ddf'' not yet implemented')
    end
    
    % Create matrices to store results
    M = nan(n, T - 1);
    MTEC = nan(n, T - 1);
    MTC = nan(n, T - 1);
    Eflag = nan(n, (T - 1) * 4);
    
    % Check if 'geomean' and the old parameter 'period' are correct
    if ~isempty(options.geomean)
        warning('''geomean'' parameter has been deprecated and will dissapear in a future realse.\n Set the new ''period'' parameter to ''geomean'' for the previous behavior of ''geomean'' = 1.\n Set ''period'' to ''base'' for the preivous behaviour of ''geomean'' = 0. See help for more information.', 'DEATOOLBOX:deprecated');        
        if options.geomean
            if ~strcmp(options.period, 'geomean' )
                error('If ''geomean'' is set to 1, ''period'' must be set to ''geomean''')
            end
        else
            if ~strcmp(options.period, 'base' )                
                error('If ''geomean'' is set to 0, ''period'' must be set to ''base''')
            end
        end
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
        tb_eff = temp_dea.eff;
        Eflag(:, (2*t - 1)) = temp_dea.exitflag(:, 1);
        
        % Compute efficiency at time t + 1
        temp_dea = dea(X(:,:,t + 1), Y(:,:,t + 1), varargin{:}, 'secondstep', 0);
        t1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 1) = temp_dea.exitflag(:, 1);
               
        % Evaluate each DMU at t + 1, with the others at base period                         
        temp_dea = dea(X(:,:,tb), Y(:,:,tb), varargin{:},...
                    'Xeval', X(:,:, t + 1),...
                    'Yeval', Y(:,:, t + 1), 'secondstep', 0);

        tbevalt1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 2) = temp_dea.exitflag(:, 1);
        
        % Additional calculatiosn for 'geomean' or 'comparison' period
        switch(options.period)
            case {'geomean','comparison'}
                % Evaluate each DMU at base period, with the others at t + 1                        
                temp_dea = dea(X(:,:,t + 1), Y(:,:,t + 1), varargin{:},...
                        'Xeval', X(:,:, tb),...
                        'Yeval', Y(:,:, tb), 'secondstep', 0);

                t1evaltb_eff = temp_dea.eff;   
                Eflag(:, (2*t - 1) + 3) = temp_dea.exitflag(:, 1);
            case 'base' 
                t1evaltb_eff = NaN;
        end
        
        % Inverse efficiencies if 'oo'
        if strcmp(options.orient, 'oo')
            tb_eff = 1 ./ tb_eff;
            t1_eff = 1 ./ t1_eff;
            tbevalt1_eff = 1 ./ tbevalt1_eff;
            t1evaltb_eff = 1 ./ t1evaltb_eff;
        end
        
        % Technical Efficiency
        MTEC(:, t) = t1_eff ./ tb_eff;
        
        % Technological Change
        switch(options.period)
            case 'geomean'
                MTC(:, t) = ((tbevalt1_eff ./ t1_eff) .* (tb_eff ./ t1evaltb_eff) ).^(1/2);
            case 'base'
                MTC(:, t) = tbevalt1_eff ./ t1_eff ;
            case 'comparison'
                MTC(:, t) = tb_eff ./ t1evaltb_eff ;
        end
        
        % Malmquist index
        M(:, t) = MTEC(:, t) .* MTC(:, t);               
        
    end
    
    % Store Malmquist results in the efficiency structure
    eff.M = M;
    eff.MTEC = MTEC;
    eff.MTC = MTC;
    eff.T = T;
    
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
    
    out.period = options.period;
    out.fixbaset = options.fixbaset;    
    
    % Custom display texts
    out.disptext_text2 = 'Malmquist:';
    out.disptext_text4 = 'M = Malmquist. MTEC = Technical Efficiency Change. MTC = Technical Change.';
    

end

