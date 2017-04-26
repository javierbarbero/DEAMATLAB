function [ out ] = deaadditsuper( X, Y, varargin )
%DEAADDITSUPER Data envelopment analysis super efficiency additive model
%   Computes data envelopment analysis super efficiency additive model 
%
%   out = DEAADDITSUPER(X, Y, Name, Value) computes data envelopment analysis 
%   super efficiency additive model with inputs X and outputs Y. Model 
%   properties are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
%   returns to sacle 'vrs'.
%   - 'rhoX': input slacks weights. Default is MIP: 1 ./ X.
%   - 'rhoY': output slacks weights. Default is MIP: 1 ./ Y.
%   - 'names': DMU names.
%
%   Example
%     
%      additsuper = deaadditsuper(X, Y, 'rts', 'vrs');
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
    
    % TODO: Additive model (other function)
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
    [n, m] = size(X);
    s = size(Y,2);
    neval = n;
    
    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    orient = options.orient;
    
    % Supperefficiency for Additive model
    if ~strcmp(orient, 'none')
        error('Additive super-efficiency model orientation msut be none.')
    end
    
    % Xeval, X and Yeval, Y must be equal in this function
    if ~isempty(options.Xeval) && size(options.Xeval) ~= size(X)
        error('Xeval and X must be equal')
    end
    
    if ~isempty(options.Yeval) && size(options.Yeval) ~= size(Y)
        error('Yeval and Y must be equal')
    end
    
    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;   
    
    % Create variable to store results
    lambda = nan(neval, n - 1);
    slackX = nan(neval, m);
    slackY = nan(neval, s);
    supereff = nan(n,1);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    Eflag = nan(neval, 1);
    
    % Returns ro scale
    rts = options.rts;
    switch(rts)
        case 'crs'
            %AeqRTS1 = [];
            %beqRTS1 = [];

            AeqRTS2super = [];
            beqRTS2super = [];
        case 'vrs'
            %AeqRTS1 = [ones(1,n), 0];
            %beqRTS1 = 1;

            AeqRTS2super = [ones(1,n - 1), zeros(1,m), zeros(1,s)];
            beqRTS2super = 1;
    end
    
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
    
    % For each DMU
    for j=1:n
        
        % ADDITIVE MODEL for each DMU
        tempdea = deaaddit(X, Y, varargin{:},...
                        'Xeval',X(j,:),...
                        'Yeval',Y(j,:));              
                
        % If it is efficeint, perform super additive model        
        if tempdea.eff < 1e-05
            
            % ADDITIVE SUPER-EFFICIENCY
            
            
            % Objective Function
            fsuper = [zeros(1,n - 1), rhoX(j,:) .* ones(1,m), rhoY(j,:) .* ones(1,s)];
            
            % Lower bound
            lbsuper = zeros(n + m + s - 1, 1);
            
            % Evaluate each DMU w.r.t all without including itself
            others = 1:n;
            others = others(others ~= j);
            
            % Constraints
            Asuper = [ X(others,:)',  -eye(m,m),  zeros(m,s);
                      -Y(others,:)', zeros(s,m), -eye(s,s)];
            bsuper = [X(j,:)'; -Y(j,:)'];
            [zsuper, ~, exitflag] = linprog(fsuper, Asuper, bsuper, AeqRTS2super, beqRTS2super, lbsuper, [], [], optimopts);
            if exitflag ~= 1
                if options.warning
                    warning('Optimization exit flag: %i', exitflag)
                end
            end
            if isempty(zsuper)
                if options.warning
                    warning('Optimization doesn''t return a result. Results set to NaN.')
                end
                zsuper = nan(n + m + s - 1, 1);                  
            end
            
            lambda(j,:) = zsuper(1: n-1);
            slackX(j,:) = zsuper(n: n + m - 1);
            slackY(j,:) = zsuper(n + m: n + m + s -1);
            supereff(j,:) = sum(rhoX(j,:) .* slackX(j,:)) + sum(rhoY(j,:) .* slackY(j,:));
            
            Xeff(j,:) = NaN;
            Yeff(j,:) = NaN;
            
            Eflag(j) = exitflag;

        else
            supereff(j) = NaN;
            lambda(j,:) = NaN;
            slackX(j,:) = NaN;
            slackY(j,:) = NaN;
            Xeff(j,:) = NaN;
            Yeff(j,:) = NaN;
            Eflag(j) = tempdea.exitflag;
        end
 
        
                
    end
    
    % Slacks structure
    slack.X = slackX;
    slack.Y = slackY;    

    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'additive-supereff', 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', supereff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/X/Y/slack.X/slack.Y/eff');

end

