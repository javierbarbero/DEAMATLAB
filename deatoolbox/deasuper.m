function [ out ] = deasuper( X, Y, varargin )
%DEASUPER Data envelopment analysis super efficiency radial and directional
%   Computes data envelopment analysis super efficiency radial and 
%   directional model
%
%   out = DEASUPER(X, Y, Name, Value) computes data envelopment analysis 
%   super efficiency model with inputs X and outputs Y. Model properties 
%   are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
%   directional distane function 'ddf'.
%   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
%   returns to sacle 'vrs'.
%   - 'Gx': input directions for 'ddf' orientation. Default is X.
%   - 'Gy': output directions for 'ddf' orientation. Default is Y.
%   - 'names': DMU names.
%
%   Example
%     
%      iosuper = deasuper(X, Y, 'orient', 'io');
%
%   See also DEAOUT, DEA, DEASCALE, DEAMALM, DEAADDIT, DEAADDITSUPER
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
    rts = options.rts;
    
    % Supperefficiency for Additive model
    if strcmp(orient, 'none')
        error('super-efficiency model must be oriented.')
    end
    
    % Xeval, X and Yeval, Y must be equal in this function
    if ~isempty(options.Xeval) && size(options.Xeval) ~= size(X)
        error('Xeval and X must be equal')
    end
    
    if ~isempty(options.Yeval) && size(options.Yeval) ~= size(Y)
        error('Yeval and Y must be equal')
    end
    
    % If DDF
    if strcmp(orient, 'ddf')
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
        
    end
    
    % Create variable to store results
    lambda = nan(neval, n - 1);
    slackX = nan(neval, m);
    slackY = nan(neval, s);
    supereff = nan(n,1);
    Xeff = nan(neval, m);
    Yeff = nan(neval, s);
    Eflag = nan(neval, 2);
    
    % For each DMU
    for j=1:n
        
        % Evaluate each DMU w.r.t all without including itself
        others = 1:n;
        others = others(others ~= j);
        
        if strcmp(orient, 'ddf')
            % DDF
            tempdea = dea(X(others, :), Y(others, :), varargin{:},...
                        'Xeval',X(j,:),...
                        'Yeval',Y(j,:),...
                        'Gx', Gx(j,:), 'Gy', Gy(j,:));   
        else
            tempdea = dea(X(others, :), Y(others, :), varargin{:},...
                        'Xeval',X(j,:),...
                        'Yeval',Y(j,:));   
        end
                    
        
           
        supereff(j) = tempdea.eff(1);
        lambda(j,:) = tempdea.lambda(1,:);
        slackX(j,:) = tempdea.slack.X(1,:);
        slackY(j,:) = tempdea.slack.Y(1,:);
        Xeff(j,:) = tempdea.Xeff(1,:);
        Yeff(j,:) = tempdea.Yeff(1,:);
        Eflag(j,:) = tempdea.exitflag;
                
    end
    
    % Slacks structure
    slack.X = slackX;
    slack.Y = slackY;    
    
    % Model
    if strcmp(orient, 'ddf')
        model = 'directional-supereff';
    else
        model = 'radial-supereff';
    end
    
    % SAVE results and input data
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', model, 'orient', orient, 'rts', rts,...
        'lambda', lambda, 'slack', slack,...
        'eff', supereff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/X/Y/eff/slack.X/slack.Y');


end

