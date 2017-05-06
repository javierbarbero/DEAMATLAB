function [ out ] = deamalmboot( X, Y, varargin )
%DEAMALMBOOT Data envelopment analysis Malmquist indices bootstrap
%   Computes data envelopment analysis Malmquist indices bootstrap
%
%   out = DEAMALMBOOT(X, Y, Name, Value) computes data envelopment analysis 
%   Malmquist indices bootstrap with inputs X and outputs Y. Model properties
%   are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
%   - 'names': DMU names.
%   - 'fixbaset': previous year 0 (default), first year 1.
%   - 'nreps': number of bootstrap replications. Default is 200.
%   - 'alpha': alpha value for confidence intervals. Default is 0.05.
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
%      iomalm = deamalmboot(X, Y, 'orient', 'io', 'nreps, 200);
%
%   See also DEAOUT, DEA, DEABOOT, DEAMALML
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
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
    [n, m, T] = size(X);
    s = size(Y,2);
    
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
        error('Malmquist index only available for ''crs'' returns to scale')
    end
    
    % Check orientation
    if strcmp(options.orient, 'ddf') || strcmp(options.orient, 'oo')
        error('Malmquist index bootstrap only for ''io'' orientation')
    end
    
    % Get number of Bootstrap replications and significance
    nreps = options.nreps;
    alph = options.alpha;
        
    % Create matrices to store results
    Mb = nan(n, T - 1);
    MTECb = nan(n, T - 1);
    MTCb = nan(n, T - 1);
    
    MB = nan(n, nreps, T - 1);
    MTECB = nan(n, nreps, T - 1);
    MTCB = nan(n, nreps, T - 1);
    EflagB = nan(n, nreps * 2, (T - 1)); % Only for T and T + 1
    
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
    
    % Original Malmquist indices
    tempmalm = deamalm(X, Y, varargin{:});
    Mo = tempmalm.eff.M;
    MTECo = tempmalm.eff.MTEC;
    MTCo = tempmalm.eff.MTC;
    
    % Bandwidth 
    % Silverman's (1986) suggestion for bivariate data
    h = (4 ./ (5 .* n)) .^ (1/6);
    
    % For each time period
    for t=1:T-1
        % Get base period
        if isempty(options.fixbaset) || options.fixbaset == 0
            tb = t;
        elseif options.fixbaset == 1
            tb = 1;
        end
        
        % A, B, and Delta matrix
        A = dea(X(:, :, tb), Y(:, :, tb), varargin{:}, 'secondstep', 0);
        A = A.eff;
        B = dea(X(:, :, t + 1), Y(:, :, t + 1), varargin{:}, 'secondstep', 0);
        B = B.eff;
        Delta = [   A,     B;
                2 - A,     B;
                2 - A, 2 - B;
                    A, 2 - B];
        
        % For each replication
        parfor i=1:nreps       
            % Draw random sample
            [DeltaStar, idx] = datasample(Delta, n);

            % delta means
            deltaMat = diag(mean(DeltaStar));
            
            % Covariance matrix
            Sigma = cov(A, B);
            SigmaR = Sigma;
            SigmaR(1, 2) = -SigmaR(1, 2);
            SigmaR(2, 1) = -SigmaR(2, 1);
            
            % Random numbers
            raMat = nan(n, 2);
            
            idxSigma = idx <= n | (idx > (2*n) & idx <= 3*n); % Drawn from [A B] or [2-A 2-B]
            raMat(idxSigma, :) = mvnrnd([0, 0], Sigma, sum(idxSigma));
            
            idxSigmaR = (idx > n & idx <= 2*n) | idx > 3*n; % Drawn from [2-A B] or [A - 2-B]
            raMat(idxSigmaR, :) = mvnrnd([0, 0], SigmaR, sum(idxSigmaR));


            % Gamma matrix   
            C = ones(n, 2);
            Gamma = (1 + h.^2).^(-1/2) .* ...
                (DeltaStar + h .* raMat - C * deltaMat) + ...
                (C * deltaMat);

            % Reflect values;
            GammaStar = Gamma;
            GammaStar(Gamma < 1) = 2 - Gamma(Gamma < 1);

            % Pseudosamples. Generate inefficient inputs or output.
            Xpseudo1 = X(:, :, tb) .* repmat(GammaStar(:, 1) ./ A, 1, m);
            Xpseudo2 = X(:, :, t + 1) .* repmat(GammaStar(:, 2) ./ B, 1, m);
            
            % MALMQUIST Index
            % Compute efficiency at base period
            temp_dea = dea(Xpseudo1, Y(:,:,tb), varargin{:}, 'secondstep', 0, 'Xeval', X(:, :, tb), 'Yeval', Y(:, :, tb) );
            tb_eff = temp_dea.eff;

            % Compute efficiency at time t + 1
            temp_dea = dea(Xpseudo2, Y(:,:,t + 1), varargin{:}, 'secondstep', 0, 'Xeval', X(:, :, t + 1), 'Yeval', Y(:, :, t + 1) );
            t1_eff = temp_dea.eff;

            % Evaluate each DMU at t + 1, with the others at base period                         
            temp_dea = dea(Xpseudo1, Y(:,:,tb), varargin{:},...
                        'Xeval', Xpseudo2,...
                        'Yeval', Y(:,:, t + 1), ...
                        'secondstep', 0);

            tbevalt1_eff = temp_dea.eff;
            
            % Additional calculatiosn for 'geomean' or 'comparison' period
            switch(options.period)
                case {'geomean','comparison'}
                    % Evaluate each DMU at t + 1, with the others at base period                         
                    temp_dea = dea(Xpseudo2, Y(:,:,t + 1), varargin{:},...
                            'Xeval',Xpseudo1,...
                            'Yeval', Y(:,:, tb));

                    t1evaltb_eff = temp_dea.eff;    
                case 'base'  
                    t1evaltb_eff = NaN;
            end            
            
            % Technical Efficiency
            MTECB(:, i, t) = t1_eff ./ tb_eff;

            % Technological Change
            switch(options.period)
                case 'geomean'
                    MTCB(:, i, t) = ((tbevalt1_eff ./ t1_eff) .* (tb_eff ./ t1evaltb_eff)).^(1/2);
                case 'base'
                    MTCB(:, i, t) = tbevalt1_eff ./ t1_eff;
                case 'comparison'
                    MTCB(:, i, t) = tb_eff ./ t1evaltb_eff ;
            end

            % Malmquist index
            MB(:, i, t) = MTECB(:, i, t) .* MTCB(:, i, t);
        
        end
        
        % Bootrstrap Technical Efficiency
        MTEC.bias(:, t) = mean(MTECB(:, :, t), 2) - MTECo(:, t);
        MTEC.b(:, t) = MTECo(:, t) - MTEC.bias(:, t);
        confInt = repelem(MTECo(:, t), 1, 2) + ...
            quantile(repmat(MTECo(:, t), 1, nreps) - MTECB(:, :, t), [0.5*alph, 1 - 0.5*alph], 2);
        MTEC.cL(:, t) = confInt(:, 1);
        MTEC.cU(:, t) = confInt(:, 2);
        
        % Technological Change
        MTC.bias(:, t) = mean(MTCB(:, :, t), 2) - MTCo(:, t);
        MTC.b(:, t) = MTCo(:, t) - MTC.bias(:, t);
        confInt = repelem(MTCo(:, t), 1, 2) + ...
            quantile(repmat(MTCo(:, t), 1, nreps) - MTCB(:, :, t), [0.5*alph, 1 - 0.5*alph], 2);
        MTC.cL(:, t) = confInt(:, 1);
        MTC.cU(:, t) = confInt(:, 2);        
        
        % Bootstrap Malmquist
        M.bias(:, t) = mean(MB(:, :, t), 2) - Mo(:, t);
        M.b(:, t) = Mo(:, t) - M.bias(:, t);
        confInt = repelem(Mo(:, t), 1, 2) + ...
            quantile(repmat(Mo(:, t), 1, nreps) - MB(:, :, t), [0.5*alph, 1 - 0.5*alph], 2);
        M.cL(:, t) = confInt(:, 1);
        M.cU(:, t) = confInt(:, 2);

    end
    
    % Store original malmquist
    M.o = Mo;
    MTEC.o = MTECo;
    MTC.o = MTCo;   
      
    % Store Malmquist results in the efficiency structure
    eff.M = M;
    eff.MTEC = MTEC;
    eff.MTC = MTC;
    eff.T = T;

    % Extract some results
    neval = NaN;
    lambda = NaN;
    slack.X = NaN;
    slack.Y = NaN;
    Xeff = NaN;
    Yeff = NaN;
    Eflag = NaN;
        
    % Save results
    out = deaout('n', n, 'neval', neval', 's', s, 'm', m,...
        'X', X, 'Y', Y, 'names', options.names,...
        'model', 'radial-malmquist-bootstrap', 'orient', options.orient, 'rts', options.rts,...
        'lambda', lambda, 'slack', slack, ...
        'eff', eff, 'Xeff', Xeff, 'Yeff', Yeff,...
        'exitflag', Eflag,...
        'dispstr', 'names/eff.M.o/eff.M.b/eff.M.cL/eff.M.cU',...
        'nreps', nreps, 'alpha', alph);
        
    out.period = options.period;
    out.fixbaset = options.fixbaset;    
    
    % Custom display texts
    out.disptext_text2 = 'Malmquist:';
    out.disptext_text4 = 'M = Malmquist. Mboot = Bootstrapped Malmquist. McLow = Lower confidence interval. McUpp: Upper confidence interval.';
    
    
end
