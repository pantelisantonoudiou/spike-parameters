function index = get_index(cell_array, match)
% index = get_index(cell_array, match)
% find index from 1D cell array
%
% example
% --------
%
% a = {1, '2', 3 ,50, 120, []}
% index = get_index(a, 120)

% preallocate index to empty
index = [];

% convert input to string
match = num2str(match);

for i = 1:length(cell_array)
    
    % convert cell array element to string
    search_str = num2str(cell_array{i});
    
    % check if element is identical to match
    if strcmp(search_str, match) == true
        index = i;
        return
    end

end