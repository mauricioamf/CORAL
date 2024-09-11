function ecModelEnzymeSubpools = getEnzymeSubpools(ecModel)
%
%   Divides promiscuous enzymes into sequentially renamed individual 
%   subenzymes and adds pseudo-reactions for enzyme subpool usage. Here. a 
%   certain enzyme Ei that catalyses a number r of reactions has its 
%   .mets and .rxns IDs renamed to Ei_(1:r). For example, the enzyme E1 
%   catalyse 3 reactions, so E1 is renamed to E1_1, E1_2 and E1_3. For 
%   each subenzyme, a enzyme usage pseudo-reaction is added (e.g. 
%   usage_prot_E1_1 for E1_1). The pseudo-reaction usage_prot_E1 that 
%   draws from prot_pool is renamed to pool_E1. The added enzyme 
%   usage pseudo-reactions for E1_1, E1_2 and E1_3 then draws from 
%   pool_E1. In this regard, all enzymes subpools should draw from the same
%   pool. This function changes only the structure of the reactions, and
%   does not affect the predictions in any way. Thus, it generates the 
%   exact same predictions as the original ecModel when using basic FBA.
%   This restructuring is necessary in order to apply the constraints of
%   ratios in the CORAL framework.
%
%   Usage
%       ecModelEnzymeSubpools = getEnzymeSubpools(ecModel)
%
%   Parameters
%       ecModel                       (struct) a GECKO3 ecModel structure
%
%   Outputs
%       ecModelEnzymeSubpools         (struct) a GECKO3 ecModel structure 
%                                     with subenzymes
%
% .. Author:
%       - Mauricio Ferreira           2023.10.02

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

%% Retrieve IDs and coefficients of mets and rxns associated with kcats
disp("Gathering enzyme data from model... (1/5)")

rxnsInd = num2cell(1:1:length(ecModel.rxns))';

[substrateMets, productMets, substrateCoeffs, productCoeffs] = getRxnData(ecModel, rxnsInd);

%% Retrives rxn IDs, GPRs and join substrates and products into tables
combinedTable = getRxnTable(ecModel, rxnsInd, substrateMets, productMets, substrateCoeffs, productCoeffs);

% Remove reactions with more than one enzyme. These reactions should have
% all been removed proviously, but one is never too careful
combinedTableSingles = combinedTable(~contains(combinedTable.genes, "and"), :);
clear combinedTable

%% Reorder enzymes by gene and then MW/kcat coefficient
% Convert doubles to cells
[row, col] = size(combinedTableSingles.stoichCoeffs); % Get the size of the cell array
for i = 1:row
    for j = 1:col
        if isnumeric(combinedTableSingles.stoichCoeffs{i, j}) % Check if the element is numeric
            combinedTableSingles.enzCoeffs{i, j} = num2cell(combinedTableSingles.stoichCoeffs{i, j}); % Convert to cell
        end
    end
end

% Extract MW/kcat coefficient from converted cells
filteredCoeffs = cell(numel(combinedTableSingles.enzCoeffs), 1);

for i = 1:length(combinedTableSingles.enzCoeffs)
    currentCell = combinedTableSingles.enzCoeffs{i};
    for j = 1:length(currentCell)
        if cell2mat(currentCell(j)) ~= 1 && cell2mat(currentCell(j)) ~= -1 && cell2mat(currentCell(j)) ~= 2 && cell2mat(currentCell(j)) ~= -2
            filteredCoeffs{i,1} = cell2mat(currentCell(j));
        end
    end 
end

combinedTableSingles.enzCoeffs = filteredCoeffs;

isScalar = zeros(numel(combinedTableSingles.enzCoeffs), 1);

for i = 1:numel(combinedTableSingles.enzCoeffs)
    isScalar(i) = isnumeric(combinedTableSingles.enzCoeffs{i}) && isscalar(combinedTableSingles.enzCoeffs{i});
end

indexesNotScalar = find(isScalar == 0);
combinedTableSingles(indexesNotScalar, :) = [];

% Extract enzyme pseudo-metabolites to another cell
enzymeMets = cell(length(combinedTableSingles.mets), 1);

for i = 1:length(combinedTableSingles.mets)
    currentCell = combinedTableSingles.mets{i};
    for j = 1:length(currentCell)
        enzymeMets{i,1} = char(currentCell(startsWith(currentCell, "prot_"), 1));
    end
end

combinedTableSingles.enzymeMets = enzymeMets;

% Reorder table first by gene, then by MW/kcat coeffs
combinedTableSingles = sortrows(combinedTableSingles,["genes" "enzCoeffs"], ["ascend" "descend"]);

% Remove rows without enzymeMets
toDel = [];
for dd = 1:length(combinedTableSingles.enzymeMets)
    if isempty(combinedTableSingles.enzymeMets{dd})
        toDel = [toDel dd];
    end
end
toDel = toDel';
combinedTableSingles(toDel,:) = [];

%% Rename enzyme pseudo-metabolites sequentially
disp("Constructing enzyme subpools... (2/5)")

enzCounts = struct();

% Loop through the cell array
for i = 1:length(combinedTableSingles.enzymeMets)
    entry = combinedTableSingles.enzymeMets{i};
    
    % Check if the entry exists in the entryCounts struct
    if isfield(enzCounts, entry)
        % If it exists, increment the count
        enzCounts.(entry) = enzCounts.(entry) + 1;
    else
        % If it doesn't exist, initialize the count to 1
        enzCounts.(entry) = 1;
    end
    
    % Rename the entry in the cell array
    combinedTableSingles.enzymeMets{i} = [entry, '_000', num2str(enzCounts.(entry))];
end

% Some subpools have 6 digits, this can cause problems later
for sbp = 1:length(combinedTableSingles.enzymeMets)
    % Check if the string ends with a 5-digit number
    if ~isempty(regexp(combinedTableSingles.enzymeMets{sbp}, '\d{6}$', 'once'))
        sixDigitNumber = regexp(combinedTableSingles.enzymeMets{sbp}, '\d{6}$', 'match', 'once');
        trimmedNumber = sixDigitNumber(2:end);
        modifiedString = regexprep(combinedTableSingles.enzymeMets{sbp}, '\d{6}$', trimmedNumber);
        combinedTableSingles.enzymeMets{sbp} = modifiedString;
    else
        combinedTableSingles.enzymeMets{sbp} = combinedTableSingles.enzymeMets{sbp};
    end
end

% Some subpools have 5 digits, this can cause problems later
for sbp = 1:length(combinedTableSingles.enzymeMets)
    % Check if the string ends with a 5-digit number
    if ~isempty(regexp(combinedTableSingles.enzymeMets{sbp}, '\d{5}$', 'once'))
        fiveDigitNumber = regexp(combinedTableSingles.enzymeMets{sbp}, '\d{5}$', 'match', 'once');
        trimmedNumber = fiveDigitNumber(2:end);
        modifiedString = regexprep(combinedTableSingles.enzymeMets{sbp}, '\d{5}$', trimmedNumber);
        combinedTableSingles.enzymeMets{sbp} = modifiedString;
    else
        combinedTableSingles.enzymeMets{sbp} = combinedTableSingles.enzymeMets{sbp};
    end
end

% Replace enzyme names in mets field
for i = 1:length(combinedTableSingles.mets)
    currentCell = combinedTableSingles.mets{i};
    for j = 1:length(currentCell)
        if convertCharsToStrings(currentCell{j}) == convertCharsToStrings(char(currentCell(startsWith(currentCell, "prot_"), 1)))
            combinedTableSingles.mets{i}{j} = combinedTableSingles.enzymeMets{i};
        end
    end
end

% Convert table to struct
equations = struct();
equations.rxns = combinedTableSingles.rxns;
equations.mets = combinedTableSingles.mets;
equations.stoichCoeffs = combinedTableSingles.stoichCoeffs;

% Metabolite struct for addMets function
metsToAdd = struct();
metsToAdd.mets = combinedTableSingles.enzymeMets;
metsToAdd.compartments = cell(length(combinedTableSingles.enzymeMets), 1);
for i = 1:length(metsToAdd.compartments)
    metsToAdd.compartments{i} = "c";
end

%% Apply changes in enzyme usage to model
disp("Renaming preexisting enzyme pseudo-metabolites... (3/5)")

% Rename current usage_prot_ pseudoreactions
enzymeIds = find(~cellfun('isempty',strfind(ecModel.rxns,'usage_prot')));
enzymeIds(end) = []; % Remove usage_prot_standard
enzymesToChange = ecModel.rxns(enzymeIds);
enzymesToChange = cellfun(@(x) strrep(x, 'usage_prot_', 'pool_'), enzymesToChange, 'UniformOutput', false);
ecModel.rxns(enzymeIds) = enzymesToChange;
ecModel.rxnNames(enzymeIds) = enzymesToChange;

metEnzymeIds = find(~cellfun('isempty',strfind(ecModel.mets,'prot_')));
metEnzymeIds(end-1:end) = []; % Remove usage_prot_standard
metsToChange = ecModel.mets(metEnzymeIds);
toDel = [];
for mets = 1:length(metsToChange)
    if ~startsWith(metsToChange{mets}, 'prot')
        toDel = [toDel mets];
    end
end
metsToChange(toDel) = [];
metEnzymeIds(toDel) = [];
metsToChange = cellfun(@(x) strrep(x, 'prot_', 'pool_'), metsToChange, 'UniformOutput', false);
ecModel.mets(metEnzymeIds) = metsToChange;
ecModel.metNames(metEnzymeIds) = metsToChange;

% Add metabolites
disp("Adding enzyme subpools pseudo-metabolites to model... (4/5)")
ecModel = addMets(ecModel, metsToAdd);

% equationProts = cell(size(equations.mets));
% for eqs = 1:length(equationProts)
%     currentEq = equations.mets{eqs};
%     equationProts{eqs,1} = currentEq{find(startsWith(currentEq, 'prot_'))};
% end

% Replace reactions with new metabolites
ecModel = changeRxns(ecModel, equations.rxns, equations, 1);
ecModel = sortIdentifiers(ecModel);

% Save checkpoint
% save('checkpoint2');
% load('checkpoint2.mat');

%% Add subpool pseudo-reactions for newly-created enzymes
disp("Creating enzyme subpools pseudo-reactions for newly-created pools... (5/5)")

enzymesToAdd = combinedTableSingles.enzymeMets;

% Generate .mets field
enzUsageSubs = {};
currentBaseName = '';
currentSubArray = {};

for i = 1:length(enzymesToAdd)
    str = enzymesToAdd{i};
    baseName = strtok(strrep(str, 'prot_', ''), '_');
    if ~strcmp(baseName, currentBaseName)
        if ~isempty(currentSubArray)
            currentSubArray{end+1} = ['prot_' currentBaseName];
            enzUsageSubs{end+1} = currentSubArray;
        end
        
        currentBaseName = baseName;
        currentSubArray = {str};
    else
        currentSubArray{end+1} = str;
    end
end
if ~isempty(currentSubArray)
    currentSubArray{end+1} = ['prot_' currentBaseName];
    enzUsageSubs{end+1} = currentSubArray;
end

enzUsageSubs = enzUsageSubs';

for n = 1:length(enzUsageSubs)
    currentEnz = enzUsageSubs{n};
    for m = 1:length(currentEnz)
        enzUsageSubs{n,1}{end} = strrep(enzUsageSubs{n,1}{end}, 'prot_', 'pool_');
    end
end

usageSubpool = {};
for i = 1:numel(enzUsageSubs)
    subarray = enzUsageSubs{i};  % Get the current subarray
    lastString = subarray{end};  % Get the last string
    
    % Create a new subarray with pairs of strings and the last string
    for j = 1:numel(subarray)-1
        newSubarray = {subarray{j}, lastString};
        usageSubpool = [usageSubpool; newSubarray];  % Append the new subarray to the output
    end
end

usageSubpoolCell = cell(length(usageSubpool),1);
for i = 1:size(usageSubpool, 1)
    % Combine the elements of the current row into a subarray
    usageSubpoolCell{i} = {usageSubpool{i, 1}, usageSubpool{i, 2}};
end

% Generate .rxns field
enzBaseNames = usageSubpool(:,1);
enzBaseNames = cellfun(@(x) strrep(x, 'prot_', 'usage_prot_'), enzBaseNames, 'UniformOutput', false);

% Generate .StoichCoeffs field
poolCoeffs = cell(length(enzBaseNames),1);

for m = 1:length(usageSubpoolCell)
    currentEnz = usageSubpoolCell{m};
    coeffs = zeros(1, length(currentEnz));
    for n = 1:length(currentEnz)
        if n < length(currentEnz)
            coeffs(n) = -1;
        elseif n == length(currentEnz)
            coeffs(n) = 1;
        end
    end
    poolCoeffs{m} = coeffs;
end

% Generate LBs and UBs
enzLBs = zeros(length(poolCoeffs), 1);
for k = 1:length(enzLBs)
    enzLBs(k) = -1000;
end

enzUBs = zeros(length(poolCoeffs), 1);

% Generate struct
rxnsToAdd = struct();
rxnsToAdd.rxns = enzBaseNames;
rxnsToAdd.rxnNames = enzBaseNames;
rxnsToAdd.mets = usageSubpoolCell;
rxnsToAdd.stoichCoeffs = poolCoeffs;
rxnsToAdd.lb = enzLBs;
rxnsToAdd.ub = enzUBs;

% Add subpool reactions to model
ecModel = addRxns(ecModel, rxnsToAdd);
ecModel = sortIdentifiers(ecModel);

% save('checkpoint3');
% load('checkpoint3.mat');

%% Final output
ecModelEnzymeSubpools = ecModel;
disp("Done!")
end