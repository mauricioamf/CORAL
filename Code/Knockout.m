%% Cleaning the workspace and the command window and adjusting settings
tic
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
clear;clc

%% Predict usage of underground reations
% load('../Models/eciML1515u_v2_stage2_DLKcat.mat');

% load('../Models/eciML1515u_CORAL_BRENDA.mat');
load('../Models/eciML1515u_CORAL_DLKcat.mat');

%% Define constraints and other parameters
enzymeTable = getEnzymeTable(model);

biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
biomassID = find(contains(model.rxns, biomassRxnID));
% vBio = 0.1;

% M9 media without carbon
exchangeIndex = find(contains(model.rxnNames, "exchange"));
exchangeIDs = model.rxns(exchangeIndex);
exchangeIDs(end) = [];

M9_components = ["EX_pi_e", "EX_co2_e", "EX_fe3_e", "EX_h_e", ...
    "EX_mn2_e", "EX_fe2_e", "EX_zn2_e", "EX_mg2_e", ...
    "EX_ca2_e", "EX_ni2_e", "EX_cu2_e", "EX_sel_e", ...
    "EX_cobalt2_e", "EX_h2o_e", "EX_mobd_e", "EX_so4_e", ...
    "EX_nh4_e", "EX_k_e", "EX_na1_e", "EX_cl_e", ...
    "EX_o2_e", "EX_tungs_e", "EX_slnt_e"];

model = changeRxnBounds(model, exchangeIDs, 0, 'l');
model = changeRxnBounds(model, M9_components, -1000, 'l');

% Block underground reactions
% undIndex = find(contains(model.rxns, 'u0'));
% undIDs = model.rxns(undIndex);
% 
% model = changeRxnBounds(model, undIDs, 0, 'b');

% Anaerobic growth
% model = changeRxnBounds(model, 'EX_o2_e', 0, 'l');

% Carbon sources
model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'b'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', -10, 'l'); %glycerol
% model = changeRxnBounds(model, 'EX_xyl__D_e', -10, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc__L_e', -10, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab__L_e', -10, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_fru_e', -10, 'l'); %fructose

% Desired product
% model = changeRxnBounds(model, 'EX_alltn_e', 0.05, 'l'); %Adenosylcobalamin

%% Set ratio constraint (E_i * K_i - E_m * K_m = 0)
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

enzNames = cell(length(enzymeTable.enzNames),1);
for e = 1:length(enzNames)
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
    currentGroupReactions = findRxnIDs(model, currentGroup);

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
        % model = addRatioReaction(model, {currentGroup{1} currentGroup{j}}, [num_j, num_i]);

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

for i = 1:length(rxnsToAdd.natEnz)
    % Extract numbers from the second column
    numFromSecondString = str2double(rxnsToAdd.undEnz{i}(end));

    % Create the third string using the specified pattern
    thirdString = sprintf('Ratio_%s_%d', rxnsToAdd.natEnz{i}, numFromSecondString);

    % Store the third string in the new cell array
    rxnsToAdd.rxns{i} = thirdString;
end

rxnsToAdd.rxns = rxnsToAdd.rxns';
rxnsToAdd.rxnNames = rxnsToAdd.rxns;

indexes11 = find(cellfun(@(x) endsWith(x, '_1_1'), rxnsToAdd.rxns));
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

model = addRxns(model, rxnsToAdd);

%% Construct and solve the problem
modelpFBA = model;
modelpFBA = buildRxnGeneMat(modelpFBA);

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));
enzUsageIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'usage_prot_')));
enzUsageIds(end) = [];

% Fix biomass
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio*(1-0.05), 'l');
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio*(1+0.05), 'u');

% Set lower and upper bounds
modelpFBA.lb(modelpFBA.lb==-Inf) = -1000;
modelpFBA.ub(modelpFBA.ub==Inf) = 1000;

% Set the objective
[~,nRxns] = size(modelpFBA.S);

modelpFBA.c = zeros(nRxns,1);
modelpFBA.c(biomassID) = 1;
% modelpFBA.c(subpoolIds) = 1;
% modelpFBA.c(enzUsageIds) = 1;
% modelpFBA.c(enzUsageIds) = enzymeTable.kcats.*modelpFBA.c(enzUsageIds);

modelpFBA.osense = -1;

% Solve the problem
SolutionEKcat = solveCobraLP(modelpFBA);
% SolutionEKcat = optimizeCbModel(modelpFBA);

% Get the solutions
SolutionEKcat.x = SolutionEKcat.full;
% SolutionEKcat.full = SolutionEKcat.x;

printFluxes(modelpFBA, SolutionEKcat.x, true);
% printFluxes(modelpFBA, SolutionEKcat.x, false);

%% Retrieve the native subenzyme IDs
% undRxns = startsWith(enzymeTable.enzRxnNames, 'u0');
% KOtable = enzymeTable(undRxns, :);
% KOtable = enzymeTable;

% numRows = size(KOtable, 1);
% columnIndex = find(strcmp(KOtable.Properties.VariableNames, 'enzNames'));
% nativeEnzs = cell(length(numRows),1);
% 
% for i = 1:numRows
%     currentEnz = KOtable{i, columnIndex}{1};
%     underscoreIndices = strfind(currentEnz, '_');
% 
%     if numel(underscoreIndices) >= 2
%         modifiedString = currentEnz(1:underscoreIndices(2)-1);
%         modifiedString = ['usage_' modifiedString '_2'];
%         nativeEnzs{i, 1} = modifiedString;
%     end
% end
% nativeEnzs = unique(nativeEnzs);

nativeEnzs = find(endsWith(enzymeTable.enzNames, '_1'));
nativeEnzs = enzymeTable.enzNames(nativeEnzs);
for e = 1:length(nativeEnzs)
    nativeEnzs{e,1} = ['usage_' nativeEnzs{e,1}];
end

%% Knockout a subenzyme_1 and solve a second problem using the predicted subpool usage as constraints
modelKO = model;
modelKO = buildRxnGeneMat(modelKO);
% enzymeTable = getEnzymeTable(modelKO);

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));
poolID = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'prot_pool_exchange')));

% modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID)*(1-0.25), 'l');
% modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID), 'u');
 
% modelKO.lb(subpoolIds) = SolutionEKcat.x(subpoolIds)*(1+0.05);
% modelKO.ub(subpoolIds) = SolutionEKcat.x(subpoolIds)*(1-0.05);
% modelKO.lb(poolID) = SolutionEKcat.x(poolID)*(1+0.05);
% modelKO.ub(poolID) = SolutionEKcat.x(poolID)*(1-0.05);

modelKO.c = zeros(nRxns,1);
modelKO.c(biomassID) = 1;
% modelKO.c(subpoolIds) = 1;
% modelKO.c(enzUsageIds) = enzymeTable.kcats.*modelKO.c(enzUsageIds);

% modelKO.c(enzymeTable.rxnUpIdx) = -coeffRatio;

modelKO.osense = -1;

% Initialize a cell array to store the results
% KOresults = cell(length(nativeEnzs), 4);
KOresults = {};

% Loop through each reaction to knockout
for i = 1:length(nativeEnzs)
    % Clone the original model to avoid modifying it
    model_knockout = modelKO;
    
    % Knockout the reaction
    % disp(['Knocking out subenzyme ' nativeEnzs{i}]);
    rxn_index = findRxnIDs(model_knockout, nativeEnzs{i});
    model_knockout.lb(rxn_index) = 0;
    model_knockout.ub(rxn_index) = 0;
    
    % Run FBA
    SolutionKO = solveCobraLP(model_knockout);
    % SolutionKO = optimizeCbModel(model_knockout);
    % SolutionKO.full = SolutionKO.x;
    
    % Check if the model grows or not
    if SolutionKO.stat ~= 0
        KOresults{i, 1} = nativeEnzs{i};
        KOresults{i, 2} = 'Feasible';
        KOresults{i, 3} = SolutionKO.full(biomassID);
        KOresults{i, 4} = SolutionKO.full;
    else
        KOresults{i, 1} = nativeEnzs{i};
        KOresults{i, 2} = 'Infeasible';
    end
end

% SolutionKO = solveCobraLP(modelKO);

% SolutionKO.x = SolutionKO.full;

% printFluxes(modelKO, SolutionKO.x, true);
% printFluxes(modelKO, SolutionKO.x, false);

toc

%% Retrieve solutions for all subenzymes across knockouts
for e = 1:length(nativeEnzs)
    currentEnz = nativeEnzs(e);
    currentEnz = currentEnz{1}(1:19-1);
    
    subEnzs = startsWith(model.rxnNames, currentEnz);
    subEnzs = model.rxnNames(subEnzs, 1);

    % subEnzUsage = cell(length(subEnzs), 2);
    subEnzUsage = {};

    for s = 1:length(subEnzs)
        subEnzID = find(endsWith(model.rxnNames, subEnzs(s)));
        subEnzUsage{s, 1} = SolutionEKcat.full(subEnzID);
        if KOresults{e, 2} == "Feasible"
            subEnzUsage{s, 2} = KOresults{e, 4}(subEnzID);
        else
            subEnzUsage{s, 2} = 0;
        end
    end

    KOresults{e, 5} = subEnzUsage;

end

%% Filter KOresults to keep only solutions where E_1 is used in SolutionEKcat
KOresultsE1used = KOresults;

for k = 1:length(KOresultsE1used)
    if KOresultsE1used{k, 5}{1,1} == 0
        KOresultsE1used{k, 5} = [];
    end
end

notUsed = any(cellfun(@isempty, KOresultsE1used), 2);
KOresultsE1used(notUsed,:) = [];

% Filter out zero growth mutants
zeroGrowth = [];
for k = 1:length(KOresultsE1used)
    if KOresultsE1used{k, 3} == 0
        zeroGrowth = [zeroGrowth k];
    end
end
KOresultsE1used(zeroGrowth,:) = [];

% Filter out knockout results with size 2x1
tofiltersize = [];
size12 = [1,2];
for k = 1:length(KOresultsE1used)
    if size(KOresultsE1used{k, 5}) == size12
        tofiltersize = [tofiltersize k];
    end
end
KOresultsE1used(tofiltersize,:) = [];

toc

%% Perform double KOs based on filtered KOresults
% DKOresults = cell(length(KOresultsE1used), 1);
DKOresults = {};

for ko = 1:length(KOresultsE1used)

    modelDKO = modelKO;
    
    % disp(['Knocking out subenzyme ' KOresultsE1used{ko,1} ' and ...']);
    rxn_index = findRxnIDs(modelDKO, KOresultsE1used{ko,1});
    
    modelDKO.lb(rxn_index) = 0;
    modelDKO.ub(rxn_index) = 0;

    for dko = 1:length(nativeEnzs)
        modelDKOsub = modelDKO;

        % disp(['... subenzyme ' nativeEnzs{dko}]);
        rxn_index2 = findRxnIDs(modelDKOsub, nativeEnzs{dko});
        modelDKOsub.lb(rxn_index2) = 0;
        modelDKOsub.ub(rxn_index2) = 0;
        
        % Run FBA
        SolutionDKO = solveCobraLP(modelDKOsub);
        % SolutionDKO = optimizeCbModel(modelDKOsub);
        % SolutionDKO.full = SolutionDKO.x;
        
        % Check if the model grows or not
        if SolutionDKO.stat ~= 0
            DKOresults{ko}{dko, 1} = nativeEnzs{dko};
            DKOresults{ko}{dko, 2} = 'Feasible';
            DKOresults{ko}{dko, 3} = SolutionDKO.full(biomassID);
            DKOresults{ko}{dko, 4} = SolutionDKO.full;
        else
            DKOresults{ko}{dko, 1} = nativeEnzs{dko};
            DKOresults{ko}{dko, 2} = 'Infeasible';
        end
    end

end

DKOresults = DKOresults';


%% Retrieve solutions for all subenzymes across double knockouts
for e = 1:length(nativeEnzs)
    currentEnz = nativeEnzs(e);
    currentEnz = currentEnz{1}(1:19-1);
    
    subEnzs = startsWith(model.rxnNames, currentEnz);
    subEnzs = model.rxnNames(subEnzs, 1);

    % DKOsubEnzUsage = cell(length(subEnzs), 2);
    DKOsubEnzUsage = {};

    for s = 1:length(subEnzs)
        subEnzID = find(endsWith(model.rxnNames, subEnzs(s)));
        DKOsubEnzUsage{s, 1} = SolutionEKcat.full(subEnzID);
        
        for k = 1:length(DKOresults) 
            if DKOresults{k,1}{e, 2} == "Feasible"
                DKOsubEnzUsage{s, 2} = DKOresults{k,1}{e, 4}(subEnzID);
            else
                DKOsubEnzUsage{s, 2} = 0;
            end
            DKOresults{k,1}{e, 5} = DKOsubEnzUsage;
        end

    end    

end

%% Filter double KO results to keep only solutions where E_1 is used in SolutionEKcat
DKOresultsE1used = DKOresults;

for k = 1:length(DKOresultsE1used)

    for dk = 1:length(DKOresultsE1used{k,1})
        if DKOresultsE1used{k,1}{dk,5}{1,1} == 0
            DKOresultsE1used{k,1}{dk,5} = [];
        end
    end

    notUsed = any(cellfun(@isempty, DKOresultsE1used{k,1}), 2);
    DKOresultsE1used{k,1}(notUsed,:) = [];

end

notUsed = any(cellfun(@isempty, KOresultsE1used), 2);
KOresultsE1used(notUsed,:) = [];

%%
toc