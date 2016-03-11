function [  ] = deadisp( out, dispstr )
%DEADISP Display data envelopment analysis results
%   Display data envelopment analysis results stored in a 'deaout'
%   structure.
%   DEADISP( out ) Display data envelopment analysis results.
%   DEADISP( out, dispstr ) Display results using the specified 'dispstr'.
%
%   Example
%       
%       io = dea(X, Y, 'orient', 'io');
%       deadisp(io);
%
%       deadisp(io, 'names/lambda/eff');
%
%   See also DEAOUT, DEA2TABLE
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 11, March, 2016
%
    
    if nargin < 2
        dispstr = out.dispstr;
    end

    % HEADER
    fprintf('_______________________________\n');
    fprintf('<strong>Data Envelopment Analysis (DEA)</strong>\n\n');
    
    % DMU and number of inputs and outputs information
    fprintf('DMUs: %i ', out.n);
    if out.neval ~= out.n && ~isnan(out.neval)
        fprintf('(%i evaluated)', out.neval)
    end
    fprintf('\n');
    fprintf('Inputs: %i     Outputs: %i ', out.m, out.s);
    if ~isnan(out.r)
        fprintf('    Undesirable: %i ', out.r);
    end
    fprintf('\n');
    
    % Model
    fprintf('Model: %s ', out.model);
    if strcmp(out.model,'allocative')
        fprintf('(%s)', out.modelalloc);
    end
    fprintf('\n');
    
    % Orientation
    fprintf('Orientation: %s ', out.orient);
    switch(out.orient)
        case 'io'
            fprintf('(Input oriented)');
        case 'oo'
            fprintf('(Output oriented)');
        case 'ddf'
            fprintf('(Directional distance function)');
    end
    fprintf('\n');
    
    % Returns to scale
    fprintf('Returns to scale: %s ', out.rts)
    switch(out.rts)
        case 'crs'
            fprintf('(Constant)')
        case 'vrs'
            fprintf('(Variable)')
        case 'scaleeff'
            fprintf('(Scale efficiency)')
    end
    fprintf('\n');
    
    % Bootstrap and significance
    if ~isnan(out.nreps)
        fprintf('Bootstrap replications: %i \n', out.nreps);
    end
    if ~isnan(out.alpha)
        fprintf('Significance level: %4.2f \n', out.alpha);
    end
    
    fprintf('\n');
        
    % Malmquist
    switch(out.model)
        case {'radial-malmquist','directional-malmquist-luenberger','radial-malmquist-bootstrap'}
        
            switch(out.model)
                case {'radial-malmquist'}
                    disp('Malmquist:');
                case {'directional-malmquist-luenberger'}
                    disp('Malmquist-Luenberger:');
            end

            if out.eff.fixbaset == 1
                disp('Base period is period 1')
            else
                disp('Base period is previous period');
            end
            disp(' ');
    end
    
    %SDisp = setSDispExp(out.orient);
    %SDisp.disp(out, dispstr);
    
    % Display table
    dispstr = strsplit(dispstr, '/');
    tabAll = [];
    for i=1:length(dispstr)          
            % Get param name
            paramstr = char(dispstr(i));
            
            % Get name and format
            [name, format] = getDEAformat(paramstr, out.orient);
            
            % Get data
            dat = eval(sprintf('out.%s', paramstr));
            if ~iscell(dat)
                % Convert to cell if not cell
                dat = num2cell(dat);                
            end
            
            % Number of columns
            ncols = size(dat, 2);
                
            % For each column in the data
            for j=1:ncols
                
                
                % Get Body
                bodyj = cellfun(@(x) sprintf(format, x), dat(:, j),'Unif',false);

                % Header
                if ncols > 1
                    % If more than 1 columns add number to name
                    namej = [name, num2str(j)];
                else
                    namej = name;
                end
                
                headerj = cellstr(namej);

                % All all together
                allThis_c = [headerj; ' '; bodyj];

                % Convert to char
                tabThis = char(allThis_c);

                % Right align
                tabThis = strjust(tabThis, 'right');

                % Add blanks to left and separator to right
                tabThis(:, 2:end + 1) = tabThis(:,:); % Move one space to right
                tabThis(:, 1) = ' '; % Add blank
                tabThis(:, end + 1) = '|'; % Add separator to right

                % Add to table
                tabAll = [tabAll tabThis];
                
            end
            
            % Replace second row with separator
            tabAll(2, :) = '-';            
            
    end
    
    disp(repelem('-', size(tabAll, 2))); % Upper Line
    disp(tabAll); % Table
    disp(repelem('-', size(tabAll, 2))); % Lower Line
    
    
    % More infor in Malmquist or Malmquist-Luenberger
    if strcmp(out.model, 'radial-malmquist')
        disp('M = Malmquist. MTEC = Technical Efficiency Change. MTC = Technical Change.')
    end
    if strcmp(out.model, 'directional-malmquist-luenberger')
        disp('ML: Malmquist-Luenberger. MLTEC: Technical Efficiency Change. MLTC: Technical Change.')
    end

    
end

