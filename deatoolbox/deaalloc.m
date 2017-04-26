function [ out ] = deaalloc( X, Y, varargin )
%DEAALLOC Data envelopment analysis allocative model
%   Computes data envelopment analysis allocative model: cost, revenue and
%   profit.
%
%   out = DEAALLOC(X, Y, Name, Value) computes data envelopment analysis 
%   allocative model (cost, revenue and profit) with inputs X and outputs Y.
%   Model properties are specified using one or more Name ,Value pair 
%   arguments.
%
%   Additional properties:
%   - 'Xprice': input prices.
%   - 'Yprice': output prices.
%   - 'Gx': input directions for profit model. Default is X.
%   - 'Gy': output directions for profit model. Default is Y.
%   - 'names': DMU names.
%
%   Example
%     
%      cost = deaalloc(X, Y, 'Xprice', W);
%      revenuew = deaalloc(X, Y, 'Yprice', P);
%      profit = deaalloc(X, Y, 'Xprice', W, 'Yprice', C);
%
%   See also DEAOUT, DEA, DEASCALE, DEAMALM, DEAADDIT, DEASUPER
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
            AeqRTS1 = [];
            beqRTS1 = [];
        otherwise
            error('''rts'' automatically set depending on the allocative model. ''crs'' for cost and revenue. ''vrs'' for profit.')
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
        model = 'allocative-cost';
        orient = 'io';
        dispstr = 'names/X/Xprice/Y/eff.T/eff.A/eff.C';
    end
    
    if isempty(W) && ~isempty(P)
        model = 'allocative-revenue';
        orient = 'oo';
        dispstr = 'names/X/Y/Yprice/eff.T/eff.A/eff.R';
    end
    
    if ~isempty(W) && ~isempty(P)
        model = 'allocative-profit';
        orient = 'ddf';
        dispstr = 'names/X/Xprice/Y/Yprice/eff.T/eff.A/eff.P';
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
    
    % Create variables to store results
    lambda = nan(neval, n);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    
    % Solve model depending on the program
    switch model
        
        % SOLO CRS
        case 'allocative-cost'
                
            % For each DMU
            for j=1:neval

                % Objective function
                f = [zeros(1,n), W(j,:)];

                % Constraints
                A = [ X', -eye(m,m);
                     -Y', zeros(s,m)];
                b = [zeros(m,1);
                    -Yeval(j,:)'];
                Aeq = AeqRTS1;
                beq = beqRTS1;
                lb = zeros(1, n + m);

                % Optimize
                z = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);

                % Get efficient inputs
                lambda(j,:) = z(1:n);
                Xeff(j,:) = z(n + 1 : end);

            end
            
            % Cost efficiency
            eff.C = sum(Xeff .* W, 2) ./ sum(X .* W, 2);
            
            % Technical efficiency. Input-oriented DEA model.
            tempdea = dea(X, Y, varargin{:}, 'orient', orient);
            eff.T = tempdea.eff;
            
            % Allocative Efficiency
            eff.A = eff.C ./ eff.T;
            
        % SOLO CRS
        case 'allocative-revenue'
            
            % For each DMU
            for j=1:neval

                % Objective function
                f = -[zeros(1,n), P(j,:)];

                % Constraints
                A = [ X', -zeros(m,s)
                     -Y',  eye(s,s)];
                b = [ Xeval(j,:)';
                     -zeros(s,1)];
                Aeq = AeqRTS1;
                beq = beqRTS1;
                lb = zeros(1, n + s);

                % Optimize
                z = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);

                % Get efficient inputs
                lambda(j,:) = z(1:n);
                Yeff(j,:) = z(n + 1 : end);

            end
            
            % Revenue efficiency
            eff.R = sum(Y .* P, 2) ./ sum(Yeff .* P, 2);
            
            % Technical efficiency. Output-oriented DEA model.
            tempdea = dea(X, Y, varargin{:}, 'orient', orient);
            eff.T = 1 ./ tempdea.eff;
            
            % Allocative efficiency.
            eff.A = eff.R ./ eff.T;
            
        % SOLO VRS
        case 'allocative-profit'
            
            % For each DMU
            for j=1:neval

                % Objective function
                f = -[zeros(1,n), -W(j,:), P(j,:)];

                % Constraints
                A = [ X',  -eye(m,m), zeros(m,s);
                     -Y', zeros(s,m),  eye(s,s)];
                %b = [ Xeval(j,:)';
                %     -Yeval(j,:)'];
                b = [ zeros(m,1);
                     -zeros(s,1)];
                %Aeq = AeqRTS1;
                %beq = beqRTS1;
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
            
            % Get directions
            Gx = options.Gx;
            Gy = options.Gy;
                        
            if length(Gx) == 1
                Gx = repmat(Gx, size(X,1), size(X,2));
            end
            
            if length(Gy) == 1
                Gy = repmat(Gy, size(Y,1), size(Y,2));
            end
            
            if isempty(Gx)
                Gx = X;
            end
            
            if isempty(Gy)
                Gy = Y;
            end
            
            % Profit efficiency
            % Cambiado signo + por -
            eff.P = ((sum(Yeff .* P, 2) - sum(Xeff .* W, 2)) - (sum(Y .* P, 2) - sum(X .* W, 2))) ./ (sum(P .* Gy, 2) + sum(W .* Gx, 2));
            %eff.P = (sum(Y .* R, 2) - sum(X .* C, 2)) ./ (sum(Yeff .* R, 2) - sum(Xeff .* C, 2));
            
            % Technical efficiency. DDF DEA model under VRS.
            tempdea = dea(X, Y, varargin{:}, 'orient', orient, 'rts', 'vrs', 'Gx', Gx, 'Gy', Gy);
            eff.T = tempdea.eff;
            
            % Allocative efficiency.
            eff.A = eff.P - eff.T; % Allocative
            
            % TODO: lambdas TEech or Profit???
            
    end

    
    % Slacks structure
    slack.X = NaN;
    slack.Y = NaN;   
    
    Eflag = NaN;
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', model, 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', dispstr,...
        'Xprice', W, 'Yprice', P);


end

