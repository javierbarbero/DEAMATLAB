function [ out ] = deaout( varargin )
%DEAOUT Generates a default deaout structure to store dea results
%   Generates a default deaout structure to store dea results
%
%   Information fields:
%   - n: number of DMU's.
%   - neval: number of evaluated DMU's.
%   - m: number of inputs.
%   - s: number of outputs.
%   - r: number of undesirable outputs.
%   - model: dea model.
%   - orient: orientation.
%   - rts: returns to scale.
%
%   Common fields:
%   - names: DMU names.
%   - X: inputs.
%   - Y: outputs.
%   - eff: efficiency measure
%   - slack.X: input slacks.
%   - slack.Y: ouput slacks.
%   - lambda: computed lambda'.
%   - Xeff: efficient X's.
%   - Yeff: efficient Y's.
%   - dual.X: input shadow prices.
%   - dual.Y: output shadow prices.
%   - dual.rts: RTS dual.
%   - exitflag: exit flags of the optimization.
%
%   Scale efficiency models:
%   - eff.crs: CRS efficiency.
%   - eff.vrs: VRS efficiency.
%   - eff.scale: Scale efficiency.
%
%   Malmquist index:
%   - eff.M: Malmquist index.
%   - eff.MTEC: Technical efficiency change.
%   - eff.MTC: Technical change.
%
%   Allocative efficiency model:
%   - W: Cost.
%   - P: Price.
%   - eff.C: Cost efficiency.
%   - eff.R: Revenue efficiency.
%   - eff.P: Profit efficiency.
%   - eff.A: Allocative efficiency.
%   - eff.T: Technical efficiency.
%
%   Undesirable outputs model;
%   - Yu: Undesirable outputs.
%   - slack.Yu: Undesirable outputs slacks.
%
%   Malmquist-Luenberger index:
%   - eff.ML: Malmquist-Luenberger index.
%   - eff.MLTEC: Technical efficiency change.
%   - eff.MLTC: Technical change.
%
%   DEA Bootstrap:
%   - eff.o: Origianl efficiency.
%   - eff.b: Bootstraped efficiency.
%   - eff.c: Efficiency confidence interval.
%
%   Example
%     
%      model = deaout();
%
%   See also DEA
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 1, September, 2016
%   
    
    
    % Parse options
    p = inputParser;
    if verLessThan('matlab', '8.2')
        addPar = @(v1,v2,v3,v4) addParamValue(v1,v2,v3,v4);
    else
        addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    end
    % Generic options
    addPar(p, 'n', NaN, @(x) isnumeric(x));
    addPar(p, 'neval', NaN, @(x) isnumeric(x));
    addPar(p, 's', NaN, @(x) isnumeric(x));
    addPar(p, 'm', NaN, @(x) isnumeric(x));
    addPar(p, 'X', NaN, @(x) isnumeric(x));
    addPar(p, 'Y', NaN, @(x) isnumeric(x));
    addPar(p,'names', [],...
                @(x) iscellstr(x));
    addPar(p,'model','radial',...
                @(x) any(validatestring(x,{'radial','radial-supereff','radial-malmquist',...
                'directional','directional-supereff','directional-undesirable','directional-malmquist-luenberger',...
                'additive','additive-supereff','additive-profit',...
                'allocative-cost','allocative-revenue','allocative-profit',...
                'radial-bootstrap','radial-malmquist-bootstrap'})));  
    addPar(p,'orient','none',...
                @(x) any(validatestring(x,{'io','oo','ddf','none','ddf_cfg'})));
    addPar(p,'rts','crs',...
                @(x) any(validatestring(x, {'crs','vrs','scaleeff'})));
    addPar(p, 'lambda', NaN, @(x) isnumeric(x));
    addPar(p, 'slack', NaN, @(x) isstruct(x));
    addPar(p, 'eff', NaN, @(x) isnumeric(x) || isstruct(x));
    addPar(p, 'Xeff', NaN, @(x) isnumeric(x));
    addPar(p, 'Yeff', NaN, @(x) isnumeric(x));
    addPar(p, 'dual', NaN, @(x) isstruct(x));
    addPar(p, 'dispstr', 'names/X/Y/slackX/slackY/eff', @(x) ischar(x));
    addPar(p, 'exitflag', NaN, @(x) isnumeric(x));
    
    % Allocative
    addPar(p, 'Xprice', NaN, @(x) isnumeric(x));
    addPar(p, 'Yprice', NaN, @(x) isnumeric(x));
    
    % Undesirable output
    addPar(p, 'r', NaN, @(x) isnumeric(x));
    addPar(p, 'Yu', NaN, @(x) isnumeric(x));
    addPar(p, 'Yueff', NaN, @(x) isnumeric(x));
    
    % Bootstrap
    addPar(p, 'nreps', NaN, @(x) isnumeric(x));
    addPar(p, 'alpha', NaN, @(x) isnumeric(x));
    
    p.parse(varargin{:})
    out = p.Results;

end

