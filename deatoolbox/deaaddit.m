function [ out ] = deaaddit( X, Y, varargin )
%DEAADDIT Data envelopment analysis additive model
%   Computes data envelopment analysis additive model
%
%   out = DEAADDIT(X, Y, Name, Value) computes data envelopment analysis
%   additive model with inputs X and outputs Y. Model properties are specified using 
%   one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
%   returns to sacle 'vrs'.
%   - 'names': DMU names.
%
%   Advanced parameters:
%   - 'Xeval: inputs to evaluate if different from X.
%   - 'Yeval': outputs to evaluate if different from Y.
%
%   Example
%     
%      add = deaaddit(X, Y, 'rts', 'vrs');
%
%   See also DEAOUT, DEA, DEAADDITSUPER
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 1, March, 2016
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    % Get number of DMUs (m), inputs (m) and outputs (s)
    [n, m] = size(X);
    s = size(Y,2);
    
    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    
    % Orientation
    orient = options.orient;
    if ~strcmp(options.orient, 'none')
        error('Additive model is non-oriented');
    end
        
    % If evaluate DMU at different X or Y
    if ~isempty(options.Xeval)
        Xeval = options.Xeval;
    else
        Xeval = X;
    end
    
    if ~isempty(options.Yeval)
        Yeval = options.Yeval;
    else
        Yeval = Y;
    end
    
    if size(Xeval,1) ~= size(Yeval,1)
        % Check size: rows
        error('Number of rows in Xref and Yref must be equal')
    end
    
    if size(Xeval,2) ~= size(X,2)
        % Check columns Xref
        error('Number of columns in Xref and X must be equal')
    end
    
    if size(Yeval,2) ~= size(Y,2)
        % Check columns Yref
        error('Number of columns in Yref and Y must be equal')
    end
    
    neval = size(Xeval,1);
    
    % RETURNS TO SCALE
    rts = options.rts;
    switch(rts)
        case 'crs'
            %AeqRTS1 = [];
            %beqRTS1 = [];
            
            AeqRTS2 = [];
            beqRTS2 = [];
        case 'vrs'
            %AeqRTS1 = [ones(1,n), 0];
            %beqRTS1 = 1;
            
            AeqRTS2 = [ones(1,n), zeros(1,m), zeros(1,s)];
            beqRTS2 = 1;
    end
        
    % OBJECTIVE FUNCTION
    % n zeros for \lambda
    % -ones for input slacks: m
    % -ones for output slacks: s    
    f = [zeros(1,n), -ones(1,m), -ones(1,s)];
    
    % LOWER BOUNDS
    % Zero for n, m, s
    lb = zeros(n + m + s, 1);
    
    % Create matrix to store results
    % Rows: n DMUs
    % Cols: n DMUs + m input slacks + s output slacks
    Z = zeros(neval, n + m + s);
    Eflag = NaN(neval, 1);
    
    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;    
    
    % OPTIMIZE: SOLVE LINEAR PROGRAMMING MODEL
    % Solve linear problem for each DMU
    for j=1:neval
        Aeq = [ X',   eye(m,m),  zeros(m,s);
                Y', zeros(s,m), -eye(s,s);               
               AeqRTS2];
        beq = [Xeval(j,:)'; Yeval(j,:)'; beqRTS2];
        [z, ~, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], [], optimopts);
        if exitflag ~= 1
            warning('Optimization exit flag: %i', exitflag)
        end
        if isempty(z)
            warning('Optimization doesn''t return a result. Results set to NaN.')
            z = nan(n + m + s, 1);                  
        end
        Z(j,:) = z;
        Eflag(j) = exitflag;
    end
    
    % Get results
    lambda = Z(:,1:n);
    slackX = Z(:, n + 1 : n + m);
    slackY = Z(:, n + m + 1 : n + m + s);      
    eff = sum(slackX, 2) + sum(slackY, 2);
    
    % Compute efficient inputs and outputs
    Xeff = Xeval - slackX;
    Yeff = Yeval + slackY;
    
    % Slacks structure
    slack.X = slackX;
    slack.Y = slackY;    
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'additive', 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag, ...
        'dispstr', 'names/X/Y/slack.X/slack.Y/eff');

end

