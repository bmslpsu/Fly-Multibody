function S = CatStructFields(dim, varargin)
%% CatStructFields: concatenates fields in stuctures
%
%   INPUT:
%       dim : dimension to concatenate
%   
%   OUPUT:
%       S : concatenated structure
%
%   Usage: S = CatStructFields(dim, struct1, struct2, ..., structn)
%

if length(varargin) > 1
    F = cellfun(@fieldnames, varargin, 'UniformOutput', false); % field names for all structures
    assert(isequal(F{:}), ...  % make sure field names are the same for all structures
        'All structures must have the same field names.')
    T = [varargin{:}]; % concatenate sructures into stucture array
    S = struct(); % new structure
    F = F{1}; % get field names from 1st stucture
    for f = 1:numel(F) % each field
        % Concatenate field
        S.(F{f}) = cat(dim, T.(F{f}));

        % If structure field is the same for all stuctures 
        % just use the first structure
        if isnumeric(S.(F{f}))
            if isreal(S.(F{f}))
                dS = diff(S.(F{f}),1,dim);
                nan_flag = all(~isnan(dS));
                dS(isnan(dS)) = 0;
                if all(all(~dS) & nan_flag)
                    S.(F{f}) = T.(F{f});
                end
            end
        elseif isstring(S.(F{f}))
            %S.(F{f}) = cat(dim, T.(F{f}));
        end
    end
else % no concatenation needed if one input
    S = varargin{1};
end

end
