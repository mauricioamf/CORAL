function ecModel = applyRatioConstraints(ecModel, enzymeTable)

% Set ratio constraint (E_i * K_i - E_m * K_m = 0)

% coeffRatio = cell(length(enzymeTable.kcats),1);
% for c = 1:length(coeffRatio)
%     if find(endsWith(enzymeTable.enzNames{c}, '_1')) == 1
%         coeffRatio{c,1} = -enzymeTable.kcats(c,1);
%     else
%         coeffRatio{c,1} = enzymeTable.kcats(c,1);
%     end
% end
% coeffRatio = cell2mat(coeffRatio);
% 
% modelpFBA.c(enzymeTable.rxnUpIdx) = coeffRatio;

% enzNames = cell(length(enzymeTable.enzNames),1);
enzNames = {};
for e = 1:length(enzymeTable.enzNames)
    enzNames{e,1} = strtrim(['usage_' enzymeTable.enzNames{e,1}]);
end

% Define the pattern using regular expression
pattern = '^usage_prot_(\w+)_(\d+)$';

% Initialize structures to store grouped strings and numbers
groupedEnzs = struct();
groupedKcats = struct();

% Loop through each string in the cell array
for i = 1:numel(enzNames)
    % Match the string against the pattern
    match = regexp(enzNames{i}, pattern, 'tokens', 'once');

    % If there is a match, group the string and number
    if ~isempty(match)
        % Extract the matched components
        id = match{1};
        number = match{2};

        % Create a group name based on the id
        groupName = ['Group_', id];

        % Check if the group already exists in the structure
        if isfield(groupedEnzs, groupName)
            % If it exists, append the string and number to the existing group
            groupedEnzs.(groupName) = [groupedEnzs.(groupName), enzNames{i}];
            groupedKcats.(groupName) = [groupedKcats.(groupName), num2cell(enzymeTable.enzCoef(i))];
        else
            % If it doesn't exist, create a new group
            groupedEnzs.(groupName) = enzNames(i);
            groupedKcats.(groupName) = num2cell(enzymeTable.enzCoef(i));
        end
    end
end

% Filter out groups with a size of 1x1
groupedEnzs = filterGroups(groupedEnzs);
groupedKcats = filterGroups(groupedKcats);

% Convert the structures to cell arrays
groupedEnzs = struct2cell(groupedEnzs);
groupedKcats = struct2cell(groupedKcats);

% Transpose cells
for cell = 1:length(groupedKcats)
    groupedEnzs{cell} = groupedEnzs{cell}';
    groupedKcats{cell} = groupedKcats{cell}';
end

% Loop over reaction groups
% rxnsToAdd = {};
rxnsToAdd = struct();
rxnsToAdd.natEnz = {};
rxnsToAdd.undEnz = {};
rxnsToAdd.numi = {};
rxnsToAdd.numj = {};

for i = 1:length(groupedEnzs)
    currentGroup = groupedEnzs{i};

    % Find reactions in the current group
    currentGroupReactions = findRxnIDs(ecModel, currentGroup);

    % Identify the index_i (the first reaction in the group)
    % index_i = currentGroupReactions(1);

    % Loop over the rest of the reactions in the group
    % for j = 2:length(currentGroupReactions)
    for j = 1:(length(currentGroupReactions) - 1)
        % Identify the index_j (subsequent reactions in the group)
        % index_j = currentGroupReactions(j);
        index_i = currentGroupReactions(j);
        index_j = currentGroupReactions(j + 1);

        % Specify the ratio values from your cell array
        % num_i = groupedKcats{i}{1}; % Assuming you have a cell array for each group

        % The index_j component of ratioValues may depend on the specific structure of your data
        % num_j = groupedKcats{i}{j};
        
        if index_i ~= index_j
            num_i = groupedKcats{i}{j};     % Assuming you have a cell array for each group
            num_j = groupedKcats{i}{j + 1}; % Assuming you have a cell array for each group

            rxnsToAdd.natEnz = [rxnsToAdd.natEnz currentGroup{1}];
            rxnsToAdd.undEnz = [rxnsToAdd.undEnz currentGroup{j}];
            rxnsToAdd.numi = [rxnsToAdd.numi num_i];
            rxnsToAdd.numj = [rxnsToAdd.numj num_j];

        end

        % Add the ratio constraint to the model
        % ecModel = addRatioReaction(ecModel, {currentGroup{1} currentGroup{j}}, [num_j, num_i]);

        % rxnsToAdd = [rxnsToAdd currentGroup{1} currentGroup{j} num_j num_i];

        % rxnsToAdd.natEnz = [rxnsToAdd.natEnz currentGroup{1}];
        % rxnsToAdd.undEnz = [rxnsToAdd.undEnz currentGroup{j}];
        % rxnsToAdd.numi = [rxnsToAdd.numi num_i];
        % rxnsToAdd.numj = [rxnsToAdd.numj num_j];
    end
end

for e = 1:length(rxnsToAdd.natEnz)
    rxnsToAdd.natEnz(e) = erase(rxnsToAdd.natEnz(e), 'usage_');
    rxnsToAdd.undEnz(e) = erase(rxnsToAdd.undEnz(e), 'usage_');
end

rxnsToAdd.mets = {};
rxnsToAdd.stoichCoeffs = {};

for met = 2:length(rxnsToAdd.natEnz)
    rxnsToAdd.mets{met,1} = [rxnsToAdd.undEnz(met) rxnsToAdd.undEnz(met-1)];
    rxnsToAdd.stoichCoeffs{met,1} = [rxnsToAdd.numj{met}, rxnsToAdd.numi{met}*-1];
end

rxnsToAdd.lb = ones(length(rxnsToAdd.natEnz),1) * -1000;
rxnsToAdd.ub = ones(length(rxnsToAdd.natEnz),1) * 1000;

if isempty(rxnsToAdd.mets{1,1}) == 1
    rxnsToAdd.mets(1) = [];
    rxnsToAdd.stoichCoeffs(1) = [];
    rxnsToAdd.lb(1) = [];
    rxnsToAdd.ub(1) = [];
end

% for i = 1:length(rxnsToAdd.natEnz)
%     % Extract numbers from the second column
%     numFromSecondString = rxnsToAdd.undEnz{i}(end-3:end);
% 
%     % Create the third string using the specified pattern
%     thirdString = sprintf('Ratio_%s_%s', rxnsToAdd.natEnz{i}, numFromSecondString);
% 
%     % Store the third string in the new cell array
%     rxnsToAdd.rxns{i} = thirdString;
% end

for i = 1:length(rxnsToAdd.mets)
    sub = rxnsToAdd.mets{i,1}{1,2}(end-3:end);
    prod = rxnsToAdd.mets{i,1}{1,1}(end-3:end);

    rxnName = sprintf('Ratio_%s_%s_%s', rxnsToAdd.natEnz{i}(1:end-5), sub, prod);
    
    rxnsToAdd.rxns{i} = rxnName;

end

rxnsToAdd.rxns = rxnsToAdd.rxns';
rxnsToAdd.rxnNames = rxnsToAdd.rxns;

indexes11 = find(cellfun(@(x) endsWith(x, '_0001_0001'), rxnsToAdd.rxns));
rxnsToAdd.rxns(indexes11) = [];
rxnsToAdd.rxnNames(indexes11) = [];
rxnsToAdd.mets(indexes11) = [];
rxnsToAdd.stoichCoeffs(indexes11) = [];
rxnsToAdd.lb(indexes11) = [];
rxnsToAdd.ub(indexes11) = [];

indexesDif = [];

% Loop through the nested cells
for i = 1:length(rxnsToAdd.mets)
    % Extract the two strings from the current nested cell
    strings = rxnsToAdd.mets{i};
    
    % Check if the patterns of the two strings are different
    if ~isequal(regexprep(strings{1}, '_\d+', ''), regexprep(strings{2}, '_\d+', ''))
        % If patterns are different, add the index to the array
        indexesDif = [indexesDif, i];
    end
end
indexesDif = indexesDif';

rxnsToAdd.rxns(indexesDif) = [];
rxnsToAdd.rxnNames(indexesDif) = [];
rxnsToAdd.mets(indexesDif) = [];
rxnsToAdd.stoichCoeffs(indexesDif) = [];
rxnsToAdd.lb(indexesDif) = [];
rxnsToAdd.ub(indexesDif) = [];

ecModel = addRxns(ecModel, rxnsToAdd);

end