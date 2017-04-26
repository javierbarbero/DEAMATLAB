function [ out ] = deaadditprofit( X, Y, varargin )
%DEAALLOC Data envelopment analysis additive profit inefficiency
%   Computes data envelopment analysis additive profit inefficiency
%
%   out = DEAADDITPROFIT(X, Y, Name, Value) computes data envelopment analysis 
%   additive profit inefficiency model with inputs X and outputs Y.
%   Model properties are specified using one or more Name ,Value pair 
%   arguments.
%
%   Additional properties:
%   - 'Xprice': input prices.
%   - 'Yprice': output prices.
%   - 'rhoX': input slacks weights. Default is MIP: 1 ./ X
%   - 'rhoY': output slacks weights. Default is MIP: 1 ./ Y.
%   - 'names': DMU names.
%
%   Example
%     
%      addprofit = deaadditprofit(X, Y, 'Xprice', W, 'Yprice', C);
%
%   See also DEAOUT, DEA, DEAADDIT, DEAALLOC
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
    
    % RETURNS TO SCALE
    rts = options.rts;
    switch(rts)
        case 'crs'
            error('Only ''vrs'' are allowed')
        otherwise
            
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
        
    % Get costs and revenues
    W = options.Xprice;
    P = options.Yprice;
    
    % Get whick kind of model we are estimating
    if ~isempty(W) && isempty(P)
        error('Both ''Xprice'' and ''Yprice'' must be specified')
    end
    
    if isempty(W) && ~isempty(P)
        error('Both ''Xprice'' and ''Yprice'' must be specified')
    end
    
    if ~isempty(W) && ~isempty(P)
        dispstr = 'names/X/Y/eff.T/eff.A/eff.P';
    end
        
    % Expand W and P if needed (if all firms have same prices and costs)
    if ~isempty(W) && size(W,1) == 1
        W = repelem(W, neval, 1);
    end
    
    if ~isempty(P) && size(P,1) == 1
        P = repelem(P, neval, 1);
    end
    
    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;
    
    % SLACKS WEIGHTS
    rhoX = options.rhoX;
    rhoY = options.rhoY;
    
    if isempty(rhoX)
        % rhoX = ones(size(X));
        % MIP
        rhoX = 1 ./ X;
    end
    if isempty(rhoY)
        % rhoY = ones(size(Y));
        % MIP
        rhoY = 1 ./ Y;
    end
       

    % TECHNICAL efficiency
    tech = deaaddit(X, Y, varargin{:}, 'rts', 'vrs', 'rhoX', rhoX, 'rhoY', rhoY);
    eff.T = tech.eff;
    
    
    % PROFIT efficiency
    
    % Create variables to store results
    lambda = nan(neval, n);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    
    % For each DMU
    for j=1:neval

        % Objective function
        f = -[zeros(1,n), -W(j,:), P(j,:)];

        % Constraints
        A = [ X',  -eye(m,m), zeros(m,s);
             -Y', zeros(s,m),  eye(s,s)];
        b = [ zeros(m,1);
             -zeros(s,1)];
        Aeq = [ones(1,n), zeros(1,m), zeros(1, s)];
        beq = 1;
        lb = zeros(1, n + m + s);

        % Optimize
        z = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);

        % Get efficient inputs
        lambda(j,:) = z(1:n);
        Xeff(j,:) = z(n + 1 : n + m);
        Yeff(j,:) = z(n + m + 1 : end);

    end
    
    numerator = ((sum(Yeff .* P, 2) - sum(Xeff .* W, 2)) - (sum(Y .* P, 2) - sum(X .* W, 2)));
    denominator = min(min([W ./ rhoX, P ./ rhoY]));
    eff.P = numerator ./ denominator;
    
    % ALLOCATIVE
    eff.A = eff.P - eff.T;
    
    % Slacks structure
    slack.X = NaN;
    slack.Y = NaN;   
    
    Eflag = NaN;
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'additive-profit', 'orient', 'none', 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', dispstr,...
        'Xprice', W, 'Yprice', P);


end


