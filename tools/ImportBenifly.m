function filedata = ImportBenifly(filename, dataLines)
%% ImportBenifly: Import data from a text file, return the data as a table.

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Time","Head","LWing","RWing","Abdomen","WBF"];
opts.VariableTypes = ["double","double", "double", "double","double","double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
filedata = readtable(filename, opts);

end