function [ out ] = deamalmluen( X, Y, Yu, varargin )
%DEAMALMLUEN Data envelopment analysis Malmquist-Luenberger indices
%   Computes data envelopment analysis Malmquist-Luenberger indices
%
%   out = DEAMALMLUEN(X, Y, Yu, Name, Value) computes data envelopment 
%   analysis Malmquist-Luenberger indices with inputs X and outputs Y. 
%   Model properties are specified using one or more Name ,Value pair 
%   arguments.
%
%   Additional properties:
%   - 'names': DMU names.
%   - 'orient': orientation. Directional distane function with undesirable
%   outputs 'ddf' (Aparicio, Pastor and Zofio, 2013), default. Directional 
%   distance function with undesirable outputs 'ddf_ccf' (Chung, Fare and 
%   Grosskopf).
%   - 'fixbaset': previous year 0 (default), first year 1.
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
%      malmluen = deamalmluen(X, Y, Yu);
%      malmluenfixed = deamalmluen(X, Y, Yu, 'fixbaset', 1);
%
%   See also DEAOUT, DEA, DEAUND, DEAMALMLUEN
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
    
    if size(Yu,3) ~= size(Y,3)
        error('Number of time periods in Y and Yu must be equal')
    end
    
    % Get number of DMUs (n), inputs (m), outputs (s), and undesirable
    % outputs (r)
    [n, m, T] = size(X);
    s = size(Y,2);
    r = size(Yu,2);
    
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
        error('Malmquist-Luenberger index only available for ''crs'' returns to scale')
    end
    
    % Replace default 'none' orientation to 'ddf'
    orient = options.orient;
    if strcmp(orient, 'none')
        orient = 'ddf';
    end        
    
    % Check orientation
    if ~strcmp(orient, 'ddf') && ~strcmp(orient, 'ddf_cfg')
        error('Malmquist-Luenberger index is for ''ddf'' or ''ddf_ccf'' with undesirable outputs')
    end
        
    % Create matrices to store results
    ML = nan(n, T - 1);
    MLTEC = nan(n, T - 1);
    MLTC = nan(n, T - 1);
    if options.geomean
        Eflag = nan(n, (T - 1) * 4);
    else
        Eflag = nan(n, (T - 1) * 3);
    end
    
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
        temp_dea = deaund(X(:,:,tb), Y(:,:,tb), Yu(:,:,tb), varargin{:}, 'secondstep', 0);
        tb_eff = temp_dea.eff;
        Eflag(:, (2*t - 1)) = temp_dea.exitflag(:, 1);
        
        % Compute efficiency at time t + 1
        temp_dea = deaund(X(:,:,t + 1), Y(:,:,t + 1), Yu(:,:,t + 1), varargin{:}, 'secondstep', 0 );
        t1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 1) = temp_dea.exitflag(:, 1);
               
        % Evaluate each DMU at t + 1, with the others at base period
        temp_dea = deaund(X(:,:,tb), Y(:,:,tb), Yu(:,:,tb), varargin{:},...
                    'Xeval', X(:,:, t + 1),...
                    'Yeval', Y(:,:, t + 1),...
                    'Yueval', Yu(:,:, t + 1), 'secondstep', 0);
        tbevalt1_eff = temp_dea.eff;
        Eflag(:, (2*t - 1) + 2) = temp_dea.exitflag(:, 1);
        
        % Additional calculatiosn for 'geomean' or 'comparison' period
        switch(options.period)
            case {'geomean','comparison'}
                % Evaluate each DMU at t + 1, with the others at base period                         
                temp_dea = deaund(X(:,:,t + 1), Y(:,:,t + 1), Yu(:,:,t + 1), varargin{:},...
                        'Xeval', X(:,:, tb),...
                        'Yeval', Y(:,:, tb),...
                        'Yueval', Yu(:, :, tb), 'secondstep', 0);

                t1evaltb_eff = temp_dea.eff;   
                Eflag(:, (2*t - 1) + 3) = temp_dea.exitflag(:, 1);
            case 'base' 
                t1evaltb_eff = NaN;
        end
        
        % Technical Efficiency
        MLTEC(:, t) = (1 + tb_eff) ./ (1 + t1_eff);

        % Technological Change
        switch(options.period)
            case 'geomean'
                MLTC(:, t) = (( (1 + t1_eff) ./ (1 + tbevalt1_eff) ) .* ((1 + t1evaltb_eff) ./ (1 + tb_eff) )).^(1/2);
            case 'base'
                MLTC(:, t) = (1 + t1_eff) ./ (1 + tbevalt1_eff);
            case 'comparison'
                MLTC(:, t) = (1 + t1evaltb_eff) ./ (1 + tb_eff);    
        end
        
        % Malmquist-Luenberger index
        ML(:, t) = MLTEC(:, t) .* MLTC(:, t);
                
    end
    
    % Store Malmquist results in the efficiency structure
    eff.ML = ML;
    eff.MLTEC = MLTEC;
    eff.MLTC = MLTC;
    eff.T = T;
   
    
    % Extract some results
    neval = NaN;
    lambda = NaN;
    slack.X = NaN;
    slack.Y = NaN;
    slack.Yu = NaN;
    Xeff = NaN;
    Yeff = NaN;
        
    % Save results
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'directional-malmquist-luenberger', 'orient', orient, 'rts', options.rts,...
        'lambda', lambda, 'slack', slack, ...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/eff.ML/eff.MLTEC/eff.MLTC',...
        'r', r, 'Yu', Yu);
        
    out.period = options.period;
    out.fixbaset = options.fixbaset;  
    
    % Custom display texts
    out.disptext_text2 = 'Malmquist-Luenberger:';
    out.disptext_text4 = 'ML: Malmquist-Luenberger. MLTEC: Technical Efficiency Change. MLTC: Technical Change.';
    

end