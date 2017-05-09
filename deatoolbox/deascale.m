function [ out ] = deascale( X, Y, varargin )
%DEASCALE Data envelopment analysis scale efficiency
%   Computes data envelopment analysis scale efficiency for radial and 
%   directional model
%
%   out = DEASCALE(X, Y, Name, Value) computes data envelopment analysis 
%   scale efficiency model with inputs X and outputs Y. Model properties
%   are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
%   directional distane function 'ddf'.
%   - 'Gx': input directions for 'ddf' orientation. Default is X.
%   - 'Gy': output directions for 'ddf' orientation. Default is Y.
%   - 'names': DMU names.
%
%   Advanced parameters:
%   - 'Xeval: inputs to evaluate if different from X.
%   - 'Yeval': outputs to evaluate if different from Y.
%
%   Example
%     
%      io_scale = deascale(X, Y, 'orient', 'io');
%
%   See also DEAOUT, DEA, DEAMALM, DEAADDIT, DEASUPER
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 9, May, 2017
%

    % Get number of DMU's
    n = size(Y,1);

    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    
    % Xeval, X and Yeval, Y must be of the same size in this function
    if ~isempty(options.Xeval) 
        if size(options.Xeval) ~= size(X)
            error('Xeval and X must be of the same size')
        end
    end
    
    if ~isempty(options.Yeval) 
        if size(options.Yeval) ~= size(Y)
            error('Yeval and Y must be of the same size')
        end
    end

    % Compute CRS DEA
    crs = dea(X, Y, varargin{:}, 'rts', 'crs');
    crs_eff = crs.eff;
    
    % Compute VRS DEA
    vrs = dea(X, Y, varargin{:}, 'rts', 'vrs');
    vrs_eff = vrs.eff;
    
    % Divide
    switch(options.orient)
        case {'io','oo'}
            scaleeff = crs_eff ./ vrs_eff;
        case {'ddf'}
            scaleeff = crs_eff - vrs_eff;
    end
    
    % Efficiency
    eff.crs = crs_eff;
    eff.vrs = vrs_eff;
    eff.scale = scaleeff;
    
    % Extract some results
    neval = vrs.neval;
    s = vrs.s;
    m = vrs.m;
    slack.X = NaN;
    slack.Y = NaN;
    
    % Exit flags
    Eflag = nan(neval, 4);
    Eflag(:, 1:2) = crs.exitflag;
    Eflag(:, 3:4) = crs.exitflag;
    
    % Save results
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'radial', 'orient', options.orient, 'rts', 'scaleeff',...
        'lambda', NaN, 'slack', slack, ...
        'eff', eff, 'Xeff', NaN, 'Yeff', NaN,...
        'exitflag', Eflag,...
        'dispstr', 'names/eff.crs/eff.vrs/eff.scale');

end

