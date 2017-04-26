function [ S, SB, pvalue, critval ] = deatestrts( X, Y, varargin )
%DEATESTRTS Data envelopment analysis test of returns to scale
%   Computes data envelopment analysis test of returns to scale (RTS)
%
%   [S, SB, pvalue, critval] = DEATEST(X, Y, Name, Value) computes data 
%   envelopment analysis test of returns to scale with inputs X and outputs
%   Y. Model properties are specified using one or more Name ,Value pair 
%   arguments. The function returns the test statistic (S), the bootstrapped
%   statistic (SB), the p-value (pvalue) and the critical value (critval).
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
%   - 'names': DMU names.
%   - 'nreps': number of bootstrap replications. Default is 200.
%   - 'alpha': alpha value for confidence intervals. Default is 0.05.
%   - 'disp': set to 1 to display results on screen. Default is 0.
%
%   Example
%     
%      [S, SB, pvalue, critval] = deatestrts(X, Y, 'orient', 'io');
%
%   See also DEAOUT, DEA, DEASCALE, DEABOOT
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 16, March, 2016
%

   % Get number of DMU's
    n = size(Y,1);

    % Get DEA options
    options = getDEAoptions(n, varargin{:});
    nreps = options.nreps;
    alph = options.alpha;
    
    % Xeval, X and Yeval, Y must be equal in this function
    if ~isempty(options.Xeval) && size(options.Xeval) ~= size(X)
        error('Xeval and X must be equal')
    end
    
    if ~isempty(options.Yeval) && size(options.Yeval) ~= size(Y)
        error('Yeval and Y must be equal')
    end

    % Compute CRS DEA
    crs = dea(X, Y, varargin{:}, 'rts', 'crs');
    crs_eff = crs.eff;
    
    % Compute VRS DEA
    vrs = dea(X, Y, varargin{:}, 'rts', 'vrs');
    vrs_eff = vrs.eff;    
    
    % Observed S
    switch(options.orient)
        case {'io'}
            S = sum(crs_eff) ./ sum(vrs_eff);
        case {'oo'}
            S = sum(1 ./ crs_eff) ./ sum(1 ./ vrs_eff);            
        case {'ddf'}
            %scaleeff = crs_eff - vrs_eff;
            error('DEA rts test not available for ''ddf''');
    end

    % Bootstrap
    rset = rng(); % Get random number generator settings
    crsB = deaboot(X, Y, varargin{:}, 'rts', 'crs', 'nreps', nreps, 'alpha', alph);
    crsB = crsB.eff.Boot;
    rng(rset); % Set rando number generation settings to the same as when computing CRS
    vrsB = deaboot(X, Y, varargin{:}, 'rts', 'vrs', 'nreps', nreps, 'alpha', alph, 'effRef', crs_eff);
    vrsB = vrsB.eff.Boot;
    
    if strcmp(options.orient, 'oo')
        crsB = 1 ./ crsB;
        vrsB = 1 ./ vrsB;
    end
    
    % Bootstrapped statistic
    SB = sum(crsB) ./ sum(vrsB);
    
    % Compute p-value
    lower = sum(SB < S);
    pvalue = (lower + 1) ./ nreps;
    
    % Critical value
    SBsorted = sort(SB);
    critval = SBsorted(floor(alph * nreps));

    % Display test results
    if options.disp
        fprintf('_______________________________\n');
        fprintf('<strong>DEA Test of RTS</strong>\n\n');
        
        fprintf('H0: Globally CRS \n')
        fprintf('H1: VRS \n\n')
        
        fprintf('Bootstrap replications: %i \n', nreps);
        fprintf('Significance level: %4.2f \n \n', alph);
        
        fprintf('S statistic: %7.4f \n', S);
        fprintf('Critical value: %7.4f \n', critval);
        fprintf('p-value: %7.4f \n', pvalue)
    end
       

end

