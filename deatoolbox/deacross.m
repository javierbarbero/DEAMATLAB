function [ out ] = deacross( X, Y, varargin )
%DEACROSS Data envelopment analysis cross efficiency model
%   Computes data envelopment analysis cross efficiency model using
%   Sexton's et al. (1986) model.
%
%   out = DEACROSS(X, Y, Name, Value) computes data envelopment analysis 
%   cross efficiency model with inputs X and outputs Y. Model properties are
%   specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io' (Default).
%   - 'model': 'linear' (Default).
%   - 'objective': minimization of the objective function, aggresive
%   approach ('aggresive'). Maximization of the objective function,
%   benevolent approach ('benevolent').
%   - 'mean': include the evaluated DMU when computing the cross efficiency 
%   score ('inclusive'); or excludes it ('exclusive'). Geometric mean 
%   including the evaluated DMU ('ginclusive'); or geometric mean excluding
%   it ('gexclusive').
%   - 'names': DMU names.
%
%
%   Example
%     
%      cross_agg = deacross(X, Y, 'objective', 'aggressive', 'mean', 'inclusive');
%      deadisp(cross_agg);
%
%      cross_ben = deacross(X, Y, 'objective', 'benevolent', 'mean', 'exclusive');
%      deadisp(cross_ben);
%
%   See also DEAOUT, DEA, DEABOOT
%
%   Copyright 2018 Inmaculada C. Alvarez, Javier Barbero, Jose L. Zofio
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 10, July, 2018
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end    
    
    % Get number of DMUs (n), inputs (m) and outputs (s)
    [n, m] = size(X);
    s = size(Y,2);

    % INPUT PARSER
    % Default optimization-options
    optimoptsdef = optimoptions('linprog','display','off', 'Algorithm','dual-simplex', 'TolFun', 1e-10, 'TolCon', 1e-7);
  
    % Parse Options
    p = inputParser;
    addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    
    % Generic options
    addPar(p,'names', cellstr(int2str((1:n)')),...
                @(x) iscellstr(x) && (length(x) == n) );
    addPar(p, 'optimopts', optimoptsdef, @(x) ~isempty(x));
    % Cross-eficiency options
    addPar(p,'orient','io',...
                @(x) any(validatestring(x,{'io','oo'})));
    addPar(p,'model','linear',...
                @(x) any(validatestring(x,{'linear'})));
    addPar(p,'objective','aggressive',...
                @(x) any(validatestring(x,{'aggressive','benevolent'})));
    addPar(p,'mean','inclusive',...
                @(x) any(validatestring(x,{'inclusive','exclusive', 'ginclusive', 'gexclusive'})));
    p.parse(varargin{:})
    options = p.Results;    
    
    % Correct names size (from row to column)
    if size(options.names, 2) > 1
        options.names = options.names';
    end
    
    %
    if strcmp(options.orient,'oo')
        error('Output oriented orientation not yet implemented');
    end
    
    % OPTIMIZATION OPTIONS:
    optimopts = options.optimopts;
    
    % MINIMIZE (AGGRESSIVE) / MAXIMIZE (BENEVOLENT) APPROACH
    switch(options.objective)
        case 'aggressive'
            appSign = 1;
        case 'benevolent'
            appSign = - 1;
    end

    % LINEAR APPROACH
    switch(options.model)
        case 'linear'
            % INPUT-OTINETED DEA
            eff = nan(n,1);
            v = nan(n,m);
            u = nan(n,s);            
            
            % For each DMU
            for j=1:n
                
                % Objective function
                f = -[zeros(1,m), Y(j,:)];
                
                % Constraints
                A = [-X, Y];
                b = [zeros(n,1)];
                Aeq = [X(j,:), zeros(1,s)];
                beq = 1;
                lb = zeros(1, m + s);
                                
                % Optimize
                [z, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                
                % Store efficiency and weights   
                v(j,:) = z(1:m);
                u(j,:) = z(m+1:m+s); 
                eff(j) = sum(z(m+1:m+s) .* Y(j,:)');
                
                
            end
            
            % AGGREGATE INPUTS AND OUTPUTs
            Xagg = sum(X) - X;
            Yagg = sum(Y) - Y;
            
            % AGGRESIVE / BENEVOLENT APPROACH
            v_ab = nan(n,m);
            u_ab = nan(n,s); 
            
            for j=1:n
                
                % Objective function
                f = appSign * [-Xagg(j,:), Yagg(j,:)];
                
                % Constraints
                A = [-X, Y];
                b = [zeros(n,1)];
                Aeq = [X(j,:), zeros(1,s);
                       -eff(j) * X(j,:), Y(j,:)];
                beq = [1;
                       0];
                lb = zeros(1, m + s);
                                
                % Optimize
                [z, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, [], [], optimopts);
                if exitflag ~= 1
                    if options.warning
                        warning('DMU %i. First Step. Optimization exit flag: %i', j, exitflag)
                    end
                end
                
                % Store weights
                v_ab(j,:) = z(1:m);
                u_ab(j,:) = z(m+1:m+s); 
                
            end
            
            % Peer-Appraisal
            PA = (u_ab * Y') ./ (v_ab * X');
            
            % Peer-Appraisal excluding self DMU
            PAex = PA;
            PAex = reshape(PAex, n*n, 1);
            PAex(1:n+1:n*n) = [];
            PAex = reshape(PAex, n-1, n);

            switch(options.mean)
                case 'inclusive'
                    crosseff = mean(PA)';
                case 'exclusive'
                    crosseff = mean(PAex)';
                case 'ginclusive'
                    crosseff = geomean(PA)';
                case 'gexclusive'
                    crosseff = geomean(PAex)';
            end
            
    end  
    
    % Generate output structure
    out.n = n;
    out.m = m;
    out.s = s;
    
    out.names = options.names;
    
    out.model = 'deacross-linear';
    out.orient = options.orient;
    out.rts = 'crs';
    
    % Eff-structure
    effs.eff      = eff;
    effs.crosseff = crosseff;
    effs.PA       = PA;
    effs.vab      = v_ab;
    effs.uab      = u_ab;
    
    % Return eff structure
    out.eff = effs;
    
    % Display string
    out.dispstr = 'names/eff.eff/eff.crosseff';
    
    % Custom display texts
    out.disptext_title = 'Data Envelopment Analysis (DEA)';
    out.disptext_text2 = 'Cross efficiency';
    switch(options.objective)
        case 'aggressive'
            out.disptext_text2 = [out.disptext_text2, ': Aggressive Approach'];
        case 'benevolent'
            out.disptext_text2 = [out.disptext_text2, ': Benevolent Approach'];
    end
    
    out.disptext_eff_eff = 'Eff';
    out.disptext_eff_crosseff = 'CrossEff';
    out.disptext_eff_PA = 'PA';
    out.disptext_eff_vab = 'v_ab';
    out.disptext_eff_uab = 'u_ab';
     
    out.disptext_text4 = 'CrossEff = Cross Efficiency';
    switch(options.mean)
        case 'inclusive'
            out.disptext_text4 = [out.disptext_text4, ' (inclusive mean)'];
        case 'exclusive'
            out.disptext_text4 = [out.disptext_text4, ' (exclusive mean)'];
        case 'ginclusive'
            out.disptext_text4 = [out.disptext_text4, ' (inclusive geometric mean)'];
        case 'gexclusive'
            out.disptext_text4 = [out.disptext_text4, ' (exclusive geometric mean)'];
    end
        
    
end
