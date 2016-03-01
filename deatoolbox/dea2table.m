function [ T ] = dea2table( out, dispstr )
%DEA2TABLE Convert 'deaout' results into a table object
%   Convert data envelopment analysis results in a 'deaout' structure into
%   a MATLAB table object.
%
%   T = DEA2TABLE( out ) Converts 'deaout' structure into a MATLAB table
%   object.
%   T = DEA2TABLE( out, dispstr ) Converts 'deaout' structure into a MATLAB
%   table object using the specified 'dispstr' structure.
%
%   Example
%       
%       io = dea2table(X, Y, 'orient', 'io');
%       T = dea2table(io);
%
%       T2 = dea2table(io, 'names/lambda/eff');
%
%   See also DEAOUT, DEADISP
%
%   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.deatoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 1, March, 2016
%

    if nargin < 2
        dispstr = out.dispstr;
    end    
    
    % Convert to table
    dispstr = strsplit(dispstr, '/');
    
    % Create empty table
    T = table();
    
    % Build Header and Body
    for i=1:length(dispstr)  
            % Get param name
            paramstr = char(dispstr(i));
            
            % Get data
            dat = eval(sprintf('out.%s', paramstr));
            
            % Append to Table
            T = [T table(dat)];
            
            % Set variable name
            [colname, ~] = getDEAformat(paramstr, out.orient);
            T.Properties.VariableNames(size(T,2)) = cellstr(colname);
            
    end
    

end

