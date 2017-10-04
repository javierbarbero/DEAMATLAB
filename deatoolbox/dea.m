function [ out ] = dea( X, Y, varargin)
%DEA Data envelopment analysis radial and directional model
%   Computes data envelopment analysis radial and directional model
%
%   out = DEA(X, Y, Name, Value) computes data envelopment analysis model
%   with inputs X and outputs Y. Model properties are specified using 
%   one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
%   directional distane function 'ddf'.
%   - 'rts': returns to scale. Constant returns to scale 'crs', variable
%   returns to scale 'vrs'.
%   - 'Gx': input directions for 'ddf' orientation. Default is Xeval.
%   - 'Gy': output directions for 'ddf' orientation. Default is Yeval.
%   - 'names': DMU names.
%   - 'secondstep': 1 to compute input and output slacks. Default is 1.
%
%   Advanced parameters:
%   - 'Xeval: inputs to evaluate if different from X.
%   - 'Yeval': outputs to evaluate if different from Y.
%
%   Example
%     
%      io = dea(X, Y, 'orient', 'io');
%      oo_vrs = dea(X, Y, 'orient', 'oo', 'rts', 'vrs');
%      ddf = dea(X, Y, 'ddf', 'Gx', X, 'Gy', Y);
%
%   See also DEAOUT, DEASCALE, DEAMALM, DEAADDIT, DEASUPER
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 1, August, 2017
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
            orient = 'ddf';
        case {'none'}
            error('Radial Model is Oriented');
    end       
    
    % RETURNS TO SCALE
    rts = options.rts;
    switch(rts)
        case 'crs'
            AeqRTS1 = [];
            beqRTS1 = [];
            
            AeqRTS2 = [];
            beqRTS2 = [];
        case 'vrs'
            AeqRTS1 = [ones(1,n), 0];
            beqRTS1 = 1;
            
            AeqRTS2 = [ones(1,n), zeros(1,m), zeros(1,s)];
            beqRTS2 = 1;
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
        error('Number of rows in Xeval and Yeval must be equal')
    end
    
    if size(Xeval,2) ~= size(X,2)
        % Check columns Xref
        error('Number of columns in Xeval and X must be equal')
    end
    
    if size(Yeval,2) ~= size(Y,2)
        % Check columns Yref
        error('Number of columns in Yeval and Y must be equal')
    end
    
    neval = size(Xeval,1);

    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;
    
    % Create variables to store results
    lambda = nan(neval, n);
    slackX = nan(neval, m);
    slackY = nan(neval, s);
    eff = nan(neval, 1);
    Eflag = nan(neval, 2);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    dualineqlin = nan(neval, m + s);
    dualeqlin = nan(neval, 1);
    
    % OPTIMIZE: SOLVE LINEAR PROGRAMMING MODEL
    switch(orient)
        case 'io'
            
            % For each DMU
            for j=1:neval
                
                % FIRST STEP:
                % Objective function
                f = [zeros(1,n), 1];
                
                % Constraints
                A = [ X', -Xeval(j,:)';
                     -Y', zeros(s,1)];
                b = [zeros(m,1);
                    -Yeval(j,:)'];
                Aeq = AeqRTS1;
                beq = beqRTS1;
                lb = zeros(1, n + 1);
                                
                % Optimize
                [z, ~, exitflag, ~, dual] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                if isempty(z)
                    if options.warning
                        warning('DMU %i. First Step. Optimization doesn''t return a result. Efficiency set to NaN.', j)
                    end
                    z = nan(n + 1, 1);  
                    dual.ineqlin = nan(1, m + s);
                    dual.eqlin = nan(1, 1);
                end
                
                % Get efficiency
                theta = z(end);
                Eflag(j, 1) = exitflag;
                eff(j) = theta;
                
                % Get dual results
                dualineqlin(j, :) = dual.ineqlin;
                if ~isempty(beqRTS1)
                    dualeqlin(j, 1) = dual.eqlin;
                else
                    dualeqlin(j, 1) = NaN;
                end                
                
                % SECOND STEP
                
                if(options.secondstep) && ~isnan(theta)
                
                    % Objective function
                    f = [zeros(1, n), -ones(1, m + s)];

                    % Constraints
                    Aeq = [ X', eye(m,m)  , zeros(m,s);
                            Y', zeros(s,m), -eye(s,s);
                           AeqRTS2];
                    beq = [theta .* Xeval(j,:)';
                            Yeval(j,:)';
                            beqRTS2];
                    lb = zeros(n + s + m, 1);

                    % Optimize
                    [z, ~, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], [], optimopts);
                    if exitflag ~= 1
                        if options.warning
                            warning('DMU %i. Second Step. Optimization exit flag: %i', j, exitflag)
                        end
                    end
                    if isempty(z)
                        if options.warning
                            warning('DMU %i. Second Step. Optimization doesn''t return a result. Results set to NaN.', j)
                        end
                        z = nan(n + m + s, 1);                   
                    end

                    % Get results
                    lambda(j,:) = z(1:n);
                    slackX(j,:) = z(n + 1 : n + m);
                    slackY(j,:) = z(n + m + 1 : n + m + s);                
                    Eflag(j, 2) = exitflag;

                    % Compute efficient inputs and outputs
                    Xeff(j,:) = repmat(eff(j), 1, m) .* Xeval(j,:) - slackX(j,:);
                    Yeff(j,:) = Yeval(j,:) + slackY(j,:);

                end
                
            end
            
            % Compute efficient inputs and outputs
            % Xeff = repmat(eff, 1, m) .* Xeval - slackX;
            % Yeff = Yeval + slackY;

            
        case 'oo'
                        
            % For each DMU
            for j=1:neval
                
                % FIRST STEP:
                % Objective function (maximize)
                f = -[zeros(1,n), 1];
                
                % Constraints
                A = [ X', zeros(m,1);
                     -Y', Yeval(j,:)'];
                b = [Xeval(j,:)';
                     zeros(s,1)];
                Aeq = AeqRTS1;
                beq = beqRTS1;
                lb = zeros(1, n + 1);

                % Optimize
                [z, ~, exitflag, ~, dual] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                if isempty(z)
                    if options.warning
                        warning('DMU %i. First Step. Optimization doesn''t return a result. Efficiency set to NaN.', j)
                    end
                    z = nan(n + 1, 1);  
                    dual.ineqlin = nan(1, m + s);
                    dual.eqlin = nan(1, 1);
                end
                
                % Get efficiency
                phi = z(end);
                eff(j) = phi;
                Eflag(j, 1) = exitflag;
                
                % Get dual results
                dualineqlin(j, :) = dual.ineqlin;
                if ~isempty(beqRTS1)
                    dualeqlin(j, 1) = dual.eqlin;
                else
                    dualeqlin(j, 1) = NaN;
                end    
                
                % SECOND STEP
                
                if(options.secondstep)  && ~isnan(phi)
                                 
                    % Objective function
                    f = -[zeros(1, n), ones(1, m + s)];

                    % Constraints
                    Aeq = [ X', eye(m,m)  ,  zeros(m,s),  ;
                            Y', zeros(s,m), -eye(s,s)  ;                       
                           AeqRTS2];
                    beq = [Xeval(j,:)'
                           phi .* Yeval(j,:)';
                           beqRTS2];
                    lb = zeros(n + s + m, 1);

                    % Optimize
                    [z, ~, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], [], optimopts);
                    if exitflag ~= 1
                        if options.warning
                            warning('DMU %i. Second Step. Optimization exit flag: %i', j, exitflag)
                        end
                    end
                    if isempty(z)
                        if options.warning
                            warning('DMU %i. Second Step. Optimization doesn''t return a result. Results set to NaN.', j)
                        end
                        z = nan(n + m + s, 1);                   
                    end

                    % Get results
                    lambda(j,:) = z(1:n);
                    slackX(j,:) = z(n + 1 : n + m);
                    slackY(j,:) = z(n + m + 1 : n + m + s);                
                    Eflag(j, 2) = exitflag;

                    % Compute efficient inputs and outputs
                    Xeff(j,:) = Xeval(j,:) - slackX(j,:);
                    Yeff(j,:) = repmat(eff(j), 1, s) .* Yeval(j,:) + slackY(j,:);

                end
                    
            end
            
            % Compute efficient inputs and outputs
            % Xeff = Xeval - slackX;
            % Yeff = repmat(eff, 1, s) .* Yeval + slackY;

            
        case 'ddf'
                        
            % Get directions
            Gx = options.Gx;
            Gy = options.Gy;
                        
            if length(Gx) == 1
                Gx = repmat(Gx, size(X,1), size(X,2));
            elseif size(Gx, 1) == 1
                Gx = repmat(Gx, size(X,1), 1);
            end
            
            if length(Gy) == 1
                Gy = repmat(Gy, size(Y,1), size(Y,2));
            elseif size(Gy, 1) == 1
                Gy = repmat(Gy, size(Y,1), 1);
            end
            
            if isempty(Gx)
                Gx = Xeval;
            end
            
            if isempty(Gy)
                Gy = Yeval;
            end
                        
            % For each DMU
            for j=1:neval
                
                % FIRST STEP:
                % Objective function (maximize)
                f = -[zeros(1,n), 1];
                
                % Constraints
                A = [ X', Gx(j,:)';
                     -Y', Gy(j,:)'];
                b = [ Xeval(j,:)';
                     -Yeval(j,:)'];
                Aeq = AeqRTS1;
                beq = beqRTS1;
                %lb = zeros(1, n + 1);
                lb = [zeros(1, n), -inf];
                
                % Optimize
                [z, ~, exitflag, ~, dual] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                if isempty(z)
                    if options.warning
                        warning('DMU %i. First Step. Optimization doesn''t return a result. Efficiency set to NaN.', j)
                    end
                    z = nan(n + 1, 1);  
                    dual.ineqlin = nan(1, m + s);
                    dual.eqlin = nan(1, 1);
                end
                
                % Get efficiency
                beta = z(end);
                eff(j) = beta;
                Eflag(j, 1) = exitflag;
                
                % Get dual results
                dualineqlin(j, :) = dual.ineqlin;
                if ~isempty(beqRTS1)
                    dualeqlin(j, 1) = dual.eqlin;
                else
                    dualeqlin(j, 1) = NaN;
                end  
                
                % SECOND STEP
                
                if(options.secondstep) && ~isnan(beta)
                                 
                    % Objective function
                    f = -[zeros(1, n), ones(1, m + s)];

                    % Constraints
                    Aeq = [ X', eye(m,m)  ,  zeros(m,s),  ;
                            Y', zeros(s,m), -eye(s,s)  ;                       
                           AeqRTS2];
                    beq = [-beta .* Gx(j,:)' + Xeval(j,:)'
                            beta .* Gy(j,:)' + Yeval(j,:)';
                           beqRTS2];
                    lb = zeros(n + s + m, 1);

                    % Optimize
                    [z, ~, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], [], optimopts);
                    if exitflag ~= 1
                        if options.warning
                            warning('DMU %i. Second Step. Optimization exit flag: %i', j, exitflag)
                        end
                    end
                    if isempty(z)
                        if options.warning
                            warning('DMU %i. Second Step. Optimization doesn''t return a result. Results set to NaN.', j)
                        end
                        z = nan(n + m + s, 1);                   
                    end

                    % Get results
                    lambda(j,:) = z(1:n);
                    slackX(j,:) = z(n + 1 : n + m);
                    slackY(j,:) = z(n + m + 1 : n + m + s);                
                    Eflag(j, 2) = exitflag;

                    % Compute efficient inputs and outputs
                    Xeff(j,:) = Xeval(j,:) - repmat(eff(j), 1, m) .* Gx(j,:) - slackX(j,:);
                    Yeff(j,:) = Yeval(j,:) + repmat(eff(j), 1, s) .* Gy(j,:) + slackY(j,:);
                end
                
            end

            
    end   
  
    % Slacks structure
    slack.X = slackX;
    slack.Y = slackY;  
    
    % Dual structure
    dual.X = dualineqlin(:, 1:m);
    dual.Y = dualineqlin(:, m+1: m+s);
    if ~isempty(beqRTS2)
        dual.rts = dualeqlin(:);
    else
        dual.rts = nan(neval,1);
    end
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'radial', 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'dual', dual,...
        'exitflag', Eflag,...
        'dispstr', 'names/X/Y/eff/slack.X/slack.Y');


end

