function [ name, format ] = getDEAformat( fieldname, orient )
%GETDEAFORMAT Private function
%   Private function
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 6, May, 2017
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
        case 'dual.X'
            name = 'dualX';
            format = fmtNumber;
        case 'dual.Y'
            name = 'dualY';
            format = fmtNumber;
        case 'dual.rts'
            name = 'dualRTS';
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
            name = 'McLow';
            format = fmtNumber;
        case 'eff.M.cU'
            name = 'McUpp';
            format = fmtNumber;
        case 'eff.MTEC.o'
            name = 'MTEC';
            format = fmtNumber;
        case 'eff.MTEC.b'
            name = 'MTECboot';
            format = fmtNumber;
        case 'eff.MTEC.cL'
            name = 'MTECcLow';
            format = fmtNumber;
        case 'eff.MTEC.cU'
            name = 'MTECcUpp';
            format = fmtNumber;
        case 'eff.MTC.o'
            name = 'MTC';
            format = fmtNumber;
        case 'eff.MTC.b'
            name = 'MTCboot';
            format = fmtNumber;
        case 'eff.MTC.cL'
            name = 'MTCcLow';
            format = fmtNumber;
        case 'eff.MTC.cU'
            name = 'MTCcUpp';
            format = fmtNumber;
        otherwise
            name = [];
            format = fmtNumber;
    end   
    


end

