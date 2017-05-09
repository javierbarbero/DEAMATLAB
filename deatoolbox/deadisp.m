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
%   LAST UPDATE: 9, May, 2017
%

    % Check if input is a structure
    if ~isstruct(out)
        error('Input is not a structure');
    end
    
    % If not custom dispstr is specified
    if nargin < 2
        % Check if the structure has a dispstr
        if isfield(out, 'dispstr')
            dispstr = out.dispstr;
        else
            error('Input structure does not have a display string, ''dispstr'', field')
        end
    end

    % TITLE
    fprintf('_______________________________\n');
    if isfield(out, 'disptext_title')
        % Display specified title
        fprintf('<strong>%s</strong>\n\n',out.disptext_title);
    else
        fprintf('<strong>Data Envelopment Analysis (DEA)</strong>\n\n');
    end
    
    % TEXT 1: Before model information
    if isfield(out, 'disptext_text1')
        fprintf('%s\n', out.disptext_text1)
    end
    
    % MODEL INFORMATION    
    % DMU and number of inputs and outputs information
    if isfield(out, 'n')
        fprintf('DMUs: %i ', out.n);
        if isfield(out, 'neval')
            if out.neval ~= out.n && ~isnan(out.neval)
                fprintf('(%i evaluated)', out.neval)
            end
        end
        fprintf('\n');
    end
    
    if isfield(out, 'm') && isfield(out, 's') 
        fprintf('Inputs: %i     Outputs: %i ', out.m, out.s);
        if isfield(out, 'r')
            if ~isnan(out.r)
                fprintf('    Undesirable: %i ', out.r);
            end
        end
        fprintf('\n');
    end
    
    % Model
    if isfield(out, 'model')
        fprintf('Model: %s ', out.model);
        fprintf('\n');
    end
    
    % Orientation
    if isfield(out, 'orient')
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
    end
    
    % Returns to scale
    if isfield(out, 'rts')
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
    end
            
    % Bootstrap and significance
    if isfield(out, 'nreps')
        if ~isnan(out.nreps)
            fprintf('Bootstrap replications: %i \n', out.nreps);
        end
    end
    if isfield(out, 'alpha')
        if ~isnan(out.alpha)
            fprintf('Significance level: %4.2f \n', out.alpha);
        end
    end
    
    fprintf('\n');
    
    % TEXT 2: After model information
    if isfield(out, 'disptext_text2')
        fprintf('%s\n', out.disptext_text2)
    end
        
    % Period (for temporal models)
    if isfield(out, 'period')
        switch(out.period)
            case 'base'
                disp('Reference period is base period');
            case 'comparison'
                disp('Reference period is comparison period');
            case 'geomean'
                disp('Geometric mean is computed');
        end
    end
    
    % Fixbase t (for temporal models)
    if isfield(out, 'fixbaset')
        if out.fixbaset == 1
            disp('Base period is period 1')
        else
            disp('Base period is previous period');
        end
        disp(' ');            
    end
    
    % TEXT 3: Before table
    if isfield(out, 'disptext_text3')
        fprintf('%s\n', out.disptext_text3)
    end  
    
    % TABLE
    dispstr = strsplit(dispstr, '/');
    tabAll = [];
    for i=1:length(dispstr)          
            % Get param name
            paramstr = char(dispstr(i));
            
            % Get name and format
            [name, format] = getDEAformat(paramstr, out.orient);
            % If no name in output structure
            if isempty(name)
                disptext_field = sprintf('disptext_%s', strrep(paramstr,'.','_'));
                if isfield(out, disptext_field)
                    % If custom name exists in the output structure use it
                    name = eval(sprintf('out.%s',disptext_field));
                else
                    % If not, display paramstr name without eff.
                    name = strrep(paramstr, 'eff.', '');
                end
            end
            
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

    % TEXT 4: After table
    if isfield(out, 'disptext_text4')
        fprintf('%s\n', out.disptext_text4)
    end
    
    
end

