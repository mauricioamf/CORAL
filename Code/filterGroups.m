function filteredGroups = filterGroups(groups)
    % Remove groups with a size of 1x1
    groupNames = fieldnames(groups);
    toRemove = false(size(groupNames));
    
    for i = 1:numel(groupNames)
        toRemove(i) = numel(groups.(groupNames{i})) == 1;
    end
    
    % Remove fields with 1x1 size
    filteredGroups = rmfield(groups, groupNames(toRemove));
end