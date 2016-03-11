function [ name, format ] = getDEAformat( fieldname, orient )
%GETDEAFORMAT Private function
%   Computes data envelopment analysis radial and directional model
%
%   out = DEA(X, Y, Name, Value) computes data envelopment analysis model
%   with inputs X and outputs Y. Model properties are specified using 
%   one or more Name ,Value pair arguments.
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
%   LAST UPDATE: 11, March, 2016
%

    if nargin < 2
        orient = 'none';
    end
    
    fmtNumber = '%7.4f';

    switch fieldname
        % Common fields
        case 'names'
            name = 'DMU';
            format = '%s';
        case 'X'
            name = 'X';
            format = fmtNumber;
        case 'Y'
            name = 'Y';
            format = fmtNumber;
        case 'eff'
            switch orient
                case 'io'
                    name = 'Theta';
                case 'oo'
                    name = 'Phi';
                case 'ddf'
                    name = 'Beta';
                otherwise
                    name = 'Eff';
            end
            format = fmtNumber'; 
        case 'slack.X'
            name = 'slackX';
            format = fmtNumber;
        case 'slack.Y'
            name = 'slackY';
            format = fmtNumber;                   
        case 'lambda'
            name = 'lambda';
            format = fmtNumber;
        case 'Xeff'
            name = 'Xeff';
            format = fmtNumber;
        case 'Yeff'
            name = 'Yeff';
            format = fmtNumber;
        case 'exitflag'
            name = 'EFlag';
            format = '%i';
        % Scale efficiency
        case 'eff.crs'
            name = 'CRS';
            format = fmtNumber;
        case 'eff.vrs'
            name = 'VRS';
            format = fmtNumber;
        case 'eff.scale'
            name = 'ScaleEff';
            format = fmtNumber;
        % Malmquist
        case 'eff.M'
            name = 'M';
            format = fmtNumber;
        case 'eff.MTEC'
            name = 'MTEC';
            format = fmtNumber;
        case 'eff.MTC';
            name = 'MTC';
            format = fmtNumber;
        % Allocative efficiency
        case 'Xprice'
            name = 'Xprice';
            format = fmtNumber;
        case 'Yprice'
            name = 'Yprice';
            format = fmtNumber;
        case 'eff.C'
            name = 'CostEff';
            format = fmtNumber;
        case 'eff.R'
            name = 'RevEff';
            format = fmtNumber;
        case 'eff.P'
            name = 'ProfEff';
            format = fmtNumber;
        case 'eff.A'
            name = 'AllocEff';
            format = fmtNumber;
        case 'eff.T'
            name = 'TechEff';
            format = fmtNumber;            
        % Undesirable outputs
        case 'Yu'
            name = 'Yu';
            format = fmtNumber;
        case 'slack.Yu'
            name = 'slackYu';
            format = fmtNumber;
        case 'Yueff'
            name = 'Yueff';
            format = fmtNumber;
        % Malmquist-Luenberger
        case 'eff.ML'
            name = 'ML';
            format = fmtNumber;
        case 'eff.MLTEC'
            name = 'MLTEC';
            format = fmtNumber;
        case 'eff.MLTC'
            name = 'MLTC';
            format = fmtNumber;
        % Bootstrap
        case 'eff.o'
            name = 'eff';
            format = fmtNumber;
        case 'eff.b'
            name = 'effboot';
            format = fmtNumber;
        case 'eff.c'
            name = 'effCI';
            format = fmtNumber;
        % Malmquist Bootstrap
        case 'eff.M.o'
            name = 'M';
            format = fmtNumber;
        case 'eff.M.b'
            name = 'Mboot';
            format = fmtNumber;
        case 'eff.M.cL'
            name = 'McL';
            format = fmtNumber;
        case 'eff.M.cU'
            name = 'McU';
            format = fmtNumber;
        case 'eff.MTEC.o'
            name = 'MTEC';
            format = fmtNumber;
        case 'eff.MTEC.b'
            name = 'MTECboot';
            format = fmtNumber;
        case 'eff.MTEC.cL'
            name = 'MTECcL';
            format = fmtNumber;
        case 'eff.MTEC.cU'
            name = 'MTECcU';
            format = fmtNumber;
        case 'eff.MTC.o'
            name = 'MTC';
            format = fmtNumber;
        case 'eff.MTC.b'
            name = 'MTCboot';
            format = fmtNumber;
        case 'eff.MTC.cL'
            name = 'MTCcL';
            format = fmtNumber;
        case 'eff.MTC.cU'
            name = 'MTCcU';
            format = fmtNumber;
        otherwise
            error('Field %s not found', fieldname)
    end   
    


end

