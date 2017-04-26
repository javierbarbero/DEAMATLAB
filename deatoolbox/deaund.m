function [ out ] = deaund( X, Y, Yu, varargin)
%DEAUND Data envelopment analysis with undesirable outputs.
%   Computes data envelopment analysis model with undesirable outputs.
%
%   out = DEAUND(X, Y, Yu, Name, Value) computes data envelopment analysis 
%   model with inputs X, outputs Y, and undesirable outputs Yu. Model 
%   properties are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'names': DMU names.
%   - 'orient': orientation. Directional distane function with undesirable
%   outputs 'ddf' (Aparicio, Pastor and Zofio, 2013), default. Directional 
%   distance function with undesirable outputs 'ddf_cfg' (Chung, Fare and 
%   Grosskopf).
%
%   Advanced parameters:
%   - 'Xeval: inputs to evaluate if different from X.
%   - 'Yeval': outputs to evaluate if different from Y.
%
%   Example
%     
%      und = deaund(X, Y, Yu);
%
%   See also DEAOUT, DEA, DEAMALMLUEN
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 27, April, 2017
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y');
    end    
    
    if size(Yu,1) ~= size(Y,1)
        error('Number of rows in Yu must be equal to number of rows in Y');
    end
    
    % Get number of DMUs (n), inputs (m), outputs (s), and undesirable
    % outputs (r)
    [n, m] = size(X);
    s = size(Y,2);
    r = size(Yu,2); % Number of undesirable outputs
    
    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    
    % Orientation
    switch(options.orient)
        case {'none'}
            % Replace default 'none' orientation to 'ddf'
            orient = 'ddf';
        case {'ddf'}
            orient = 'ddf';
        case {'ddf_cfg'}
            orient = 'ddf_cfg';
        otherwise
            error('Orientation for the undesarible outputs model must be ddf');
    end 
    
    % Distance functions
    if ~isempty(options.Gx) || ~isempty(options.Gy)
        error('Distance functions Gx, Gy, and Gyu are automatically assigned in undesirable outputs dea model.');
    end
    
    % RETURNS TO SCALE
    rts = options.rts;
    switch(rts)
        case 'crs'
            AeqRTS1 = [];
            beqRTS1 = [];
            
            AeqRTS2 = [];
            beqRTS2 = [];
        otherwise
            error('DEA model with undesirable outputs only available for crs');
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
    
    if ~isempty(options.Yueval)
        Yueval = options.Yueval;
    else
        Yueval = Yu;
    end
    
    if size(Xeval,1) ~= size(Yeval,1)
        % Check size: rows
        error('Number of rows in Xref and Yref must be equal')
    end
    
    if size(Yueval,1) ~= size(Yeval,1)
        % Check size: rows
        error('Number of rows in Yuref and Yref must be equal')
    end
    
    if size(Xeval,2) ~= size(X,2)
        % Check columns Xref
        error('Number of columns in Xref and X must be equal')
    end
    
    if size(Yeval,2) ~= size(Y,2)
        % Check columns Yref
        error('Number of columns in Yref and Y must be equal')
    end
    
    if size(Yueval,2) ~= size(Yu,2)
        % Check columns Yref
        error('Number of columns in Yuref and Yu must be equal')
    end
    
    neval = size(Xeval,1);

    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;
    
    % Create variables to store results
    lambda = nan(neval, n);
    slackX = nan(neval, m);
    slackY = nan(neval, s);
    slackYu = nan(neval, r);
    eff = nan(neval, 1);
    Eflag = nan(neval, 2);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    Yueff = nan(neval, r);
    
    % OPTIMIZE: SOLVE LINEAR PROGRAMMING MODEL
    switch(orient)
        
        case 'ddf'
            % (Aparicio, Pastor and Zofio, 2013)
                        
            % Get directions
            %G = options.ddfG;
            %H = options.ddfH;
            Gx = zeros(n,m);
            Gy = Yeval;
            Gyu = Yueval;
                        
            if length(Gx) == 1
                Gx = repmat(Gx, size(X,1), size(X,2));
            end
            
            if length(Gy) == 1
                Gy = repmat(Gy, size(Y,1), size(Y,2));
            end
            
            if length(Gyu) == 1
                Gyu = repmat(Gyu, size(Y,1), size(Y,2));
            end
            
            maxYu = max(max([Yu; Yueval]));
 
           
            % For each DMU
            for j=1:neval
                
                % FIRST STEP:
                % Objective function (maximize)
                f = -[zeros(1,n), 1];
                
                % Constraints                
                A = [ X', Gx(j,:)';
                     -Y', Gy(j,:)';
                     Yu', Gyu(j,:)';
                     zeros(r,n), -Yueval(j,:)']; % Undesirible max condition
                b = [ Xeval(j,:)';
                     -Yeval(j,:)';
                     Yueval(j,:)';
                     maxYu - Yueval(j,:)';]; % Undesirible max condition
                 
                Aeq = AeqRTS1;
                beq = beqRTS1;
                %lb = zeros(1, n + 1);
                lb = [zeros(1, n), -inf];
    
                
                % Optimize
                [z, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);                
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                if isempty(z)
                    if options.warning
                        warning('DMU %i. First Step. Optimization doesn''t return a result in First Step. Efficiency set to NaN.', j) 
                    end
                    z = nan(n + 1, 1);   
                end
                
                % Get efficiency
                beta = z(end);
                eff(j) = beta;
                Eflag(j, 1) = exitflag;
                
                % SECOND STEP                
                if(options.secondstep) && ~isnan(beta)

                    % Objective function
                    f = -[zeros(1, n), ones(1, m + s + r)];

                    % Constraints
                    Aeq = [ X', eye(m,m)  ,  zeros(m,s), zeros(m,r) ;
                            Y', zeros(s,m), -eye(s,s)  , zeros(s,r) ;     
                           Yu', zeros(r,m),  zeros(r,s), eye(r,r);
                           AeqRTS2];
                    beq = [-beta .* Gx(j,:)' + Xeval(j,:)'
                            beta .* Gy(j,:)' + Yeval(j,:)';
                           -beta .* Gyu(j,:)'+ Yueval(j,:)';
                           beqRTS2];
                    lb = zeros(n + s + m + r, 1);
                    %lb = [zeros(1, n), -inf(1, s + m + r)];

                    % Optimize
                    z = linprog(f, [], [], Aeq, beq, lb, [], [], optimopts);
                    if exitflag ~= 1
                        if options.warning
                            warning('DMU %i. Second Step. Optimization exit flag: %i', j, exitflag)
                        end
                    end
                    if isempty(z)
                        if options.warning
                            warning('DMU %i. Second Step. Optimization doesn''t return a result. Results set to NaN.', j)
                        end
                        z = nan(n + m + s + r, 1);                   
                    end

                    % Get results
                    lambda(j,:) = z(1:n);
                    slackX(j,:) = z(n + 1 : n + m);
                    slackY(j,:) = z(n + m + 1 : n + m + s);   
                    slackYu(j,:) = z(n + m + s + 1 : n + m + s + r);
                    Eflag(j, 2) = exitflag;

                    % Compute efficient inputs and outputs
                    Xeff(j,:) = Xeval(j,:) - repmat(eff(j), 1, m) .* Gx(j,:) - slackX(j,:);
                    Yeff(j,:) = Yeval(j,:) + repmat(eff(j), 1, s) .* Gy(j,:) + slackY(j,:);
                    Yueff(j,:) = Yueval(j,:) - repmat(eff(j), 1, r) .* Gyu(j,:) - slackYu(j,:);
                end
                
            end
            
        case 'ddf_cfg'
            
            % (Chung, Fare and Grosskopf)            
            
            % For each DMU
            for j=1:neval
                % FIRST STEP:
                % Objective function (maximize)
                f = -[zeros(1,n), 1];

                % Constraints                                  
                A = [ X', zeros(m,1);                
                     -Y', Yeval(j,:)'];
                b = [Xeval(j,:)'; -Yeval(j,:)'; ];                
                Aeq = [Yu', Yueval(j,:)'];
                beq = Yueval(j,:)';           
                lb = [zeros(1, n), -inf];
  

                % Optimize
                [z, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);                
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                if isempty(z)
                    if options.warning
                        warning('DMU %i. First Step. Optimization doesn''t return a result in First Step. Efficiency set to NaN.', j)
                    end
                    z = nan(n + 1, 1);                    
                end                

                % Get efficiency
                beta = z(end);
                eff(j) = beta;
                Eflag(j, 1) = exitflag;
                
                % SECOND STEP
                % Not available
                lambda(j,:) = z(1:n);
                slackX(j,:) = nan(1, m);
                slackY(j,:) = nan(1, s);   
                slackYu(j,:) = nan(1, r);
                eff(j) = beta;
                Eflag(j, 2) = nan(1, 1);
                
                % Compute efficient inputs and outputs
                Xeff(j,:) = nan(1, m);
                Yeff(j,:) = nan(1, s);
                Yueff(j,:) = nan(1, r);
                
            end
            
            
    end
    
    % Slacks structure
    slack.X = slackX;
    slack.Y = slackY;
    slack.Yu = slackYu;
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'directional-undesirable', 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack, ...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/X/Y/Yu/eff/slack.X/slack.Y/slack.Yu',...
        'r', r, 'Yu', Yu, 'Yueff', Yueff);
 

end

