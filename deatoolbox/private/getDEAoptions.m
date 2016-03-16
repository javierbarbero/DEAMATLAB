function [ options ] = getDEAoptions( n, varargin )
%GETDEAOPTIONS Private function
%   Private function
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 16, March, 2016
%
    
    % Default optimization-options
    optimoptsdef = optimoptions('linprog','display','off', 'Algorithm','dual-simplex');

    % Parse Options
    p = inputParser;
    if verLessThan('matlab', '8.2')
        addPar = @(v1,v2,v3,v4) addParamValue(v1,v2,v3,v4);
    else
        addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    end
    
    % Generic options
    addPar(p,'names', cellstr(int2str((1:n)')),...
                @(x) iscellstr(x) && (length(x) == n) );
    addPar(p, 'optimopts', optimoptsdef, @(x) isstruct(x));
    addPar(p, 'disp', 0, @(x) ismember(x, [0, 1]) );
    % Radial models
    addPar(p,'orient','none',...
                @(x) any(validatestring(x,{'io','oo','ddf','none'})));
    addPar(p,'rts','crs',...
                @(x) any(validatestring(x, {'crs','vrs'})));
    addPar(p, 'Gx', [], @(x) isnumeric(x) );
    addPar(p, 'Gy', [], @(x) isnumeric(x) );   
    addPar(p, 'secondstep', 1, @(x) isnumeric(x) );
    % Special options (used by other functions)
    addPar(p, 'Xeval', [], @(x) isnumeric(x)  );
    addPar(p, 'Yeval', [], @(x) isnumeric(x)  );
    % Additive
    addPar(p, 'rhoX', [], @(x) isnumeric(x) );
    addPar(p, 'rhoY', [], @(x) isnumeric(x) );
    % Malmquist Index
    addPar(p, 'fixbaset', [], @(x) ismember(x, [0, 1]));
    addPar(p, 'geomean', 1, @(x) ismember(x, [0, 1]));
    % Allocative Efficiency
    addPar(p, 'Xprice', [], @(x) isnumeric(x));
    addPar(p, 'Yprice', [], @(x) isnumeric(x));
    % Undesarible outputs
    addPar(p, 'Yueval', [], @(x) isnumeric(x));
    % Bootstrap
    addPar(p, 'nreps', 200, @(x) isnumeric(x) & (x > 0) );
    addPar(p, 'alpha', 0.05, @(x) isnumeric(x) & (x > 0));
    addPar(p, 'effRef', [], @(x) isnumeric(x));
    
    p.parse(varargin{:})
    options = p.Results;
    
    % Correct names size (from row to column)
    if size(options.names, 2) > 1
        options.names = options.names';
    end
    

end

