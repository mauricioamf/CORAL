%% Cleaning the workspace and the command window and adjusting settings
changeCobraSolver('gurobi', 'LP');
changeCobraSolverParams('LP', 'feasTol', 1e-9);
clear;clc

%% Predict usage of underground reations
% load('../Models/eciML1515u_CORAL_DLKcat.mat');
% load('../Models/eciML1515u_CORAL_TurNuP.mat');

%% Define constraints and other parameters
enzymeTable = getEnzymeTable(model);

biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
biomassID = find(contains(model.rxns, biomassRxnID));

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'l'); % D-glucose
% vBio = 0.11;

% model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'b'); % glucose

% model = changeRxnBounds(model, 'EX_glyc_e', -1.5, 'l'); % glycerol
% vBio = 0.08;

% model = changeRxnBounds(model, 'EX_xyl__D_e', -1.5, 'l'); % xylose
% vBio = 0.09;

% model = changeRxnBounds(model, 'EX_fuc__L_e', -1, 'l'); % fucose
% vBio = 0.07;

% model = changeRxnBounds(model, 'EX_arab__L_e', -1.5, 'l'); % arabinose
% vBio = 0.09;

% model = changeRxnBounds(model, 'EX_fru_e', -1.4, 'l'); % fructose
% vBio = 0.1;

filename = "SKO_redist_Delta_Ratio_glucose.csv";

SolutionVbio = optimizeCbModel(model);
vBio = SolutionVbio.f;
printFluxes(model, SolutionVbio.x, true);

%% Construct and solve the problem
modelpFBA = model; 

modelpFBA = buildRxnGeneMat(modelpFBA);
modelpFBA = applyRatioConstraints(modelpFBA, enzymeTable);

poolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'pool_')));
enzUsageIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'usage_prot_')));
enzUsageIds(end) = [];

% Fix biomass
modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio, 'l');
modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio, 'u');

% Set lower and upper bounds
modelpFBA.lb(modelpFBA.lb==-Inf) = -1000;
modelpFBA.ub(modelpFBA.ub==Inf) = 1000;

% Set the objective
[~,nRxns] = size(modelpFBA.S);

modelpFBA.c = zeros(nRxns,1);
% modelpFBA.c(biomassID) = 1;
modelpFBA.c(poolIds) = 1;
% modelpFBA.c(enzUsageIds) = 1;
% modelpFBA.c(enzUsageIds) = enzymeTable.kcats.*modelpFBA.c(enzUsageIds);

modelpFBA.osense = 1;

% Solve the problem
SolutionEKcat = solveCobraLP(modelpFBA);
% SolutionEKcat = optimizeCbModel(modelpFBA);

% Get the solutions
SolutionEKcat.x = SolutionEKcat.full;
% SolutionEKcat.full = SolutionEKcat.x;

printFluxes(modelpFBA, SolutionEKcat.x, true);
% printFluxes(modelpFBA, SolutionEKcat.x, false);

%% Retrieve the native subenzyme IDs
nativeEnzs = find(endsWith(enzymeTable.enzNames, ['_0001' ...
    '']));
nativeEnzs = enzymeTable.enzNames(nativeEnzs);
for e = 1:length(nativeEnzs)
    nativeEnzs{e,1} = ['usage_' nativeEnzs{e,1}];
end

%% Knockout a subenzyme_1 and solve a second problem using the predicted subpool usage as constraints
tic

modelKO = model;
modelKO = buildRxnGeneMat(modelKO);
modelKO = applyRatioConstraints(modelKO, enzymeTable);

poolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'pool_')));
totalPoolID = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'prot_pool_exchange')));

modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID)*(1-0.05), 'l');
modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID), 'u');

modelKO.lb(poolIds) = SolutionEKcat.x(poolIds)*(1+0.05);
modelKO.ub(poolIds) = SolutionEKcat.x(poolIds)*(1-0.05);
modelKO.lb(totalPoolID) = SolutionEKcat.x(totalPoolID)*(1+0.05);
modelKO.ub(totalPoolID) = SolutionEKcat.x(totalPoolID)*(1-0.05);

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
    % SolutionKO = solveCobraLP(model_knockout);
    SolutionKO = optimizeCbModel(model_knockout);
    SolutionKO.full = SolutionKO.x;
    
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

% Filter out infeasibilities
infeasible = [];
for k = 1:length(KOresultsE1used)
    if KOresultsE1used{k, 2} == "Infeasible"
        infeasible = [infeasible k];
    end
end
KOresultsE1used(infeasible,:) = [];

% % Filter out knockout results with size 2x1
% tofiltersize = [];
% size12 = [1,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size12
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];

% % Filter out knockout results with size 2x2
% tofiltersize = [];
% size22 = [2,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size22
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];

% % Filter out knockout results with size 3x2
% tofiltersize = [];
% size32 = [3,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size32
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];

% % Filter out knockout results with size 4x2
% tofiltersize = [];
% size42 = [4,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size42
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];

% Keep only desired solutions 
for k = 1:length(KOresultsE1used)
    if KOresultsE1used{k, 5}{1,1} == 0
        KOresultsE1used{k, 5} = [];
    end
end

notUsed = any(cellfun(@isempty, KOresultsE1used), 2);
KOresultsE1used(notUsed,:) = [];

% % Filter out zero growth mutants
% zeroGrowth = [];
% for k = 1:length(KOresultsE1used)
%     if KOresultsE1used{k, 3} == 0
%         zeroGrowth = [zeroGrowth k];
%     end
% end
% KOresultsE1used(zeroGrowth,:) = [];

%%
save("SKO_redist_Delta_Ratio_glucose.mat", "-v7.3")

%% Get subpool deltas
subpools = {};
enzymeWT = [];
enzymeKO = [];
enzymeDelta = [];

for e = 1:length(KOresultsE1used(:,1))
    currentEnz = KOresultsE1used{e,5};

    for s = 1:length(currentEnz)
        enzymeWT = [enzymeWT currentEnz{s,1}];
        enzymeKO = [enzymeKO currentEnz{s,2}];
        enzymeDelta = [enzymeDelta currentEnz{s,1} - currentEnz{s,2}];
        subpools = [subpools KOresultsE1used{e,1}];
    end

end

enzymeWT = enzymeWT';
enzymeKO = enzymeKO';
enzymeDelta = enzymeDelta';

% enzymeDelta = enzymeDelta(enzymeDelta ~= 0);
% enzymeDelta = enzymeDelta(~cellfun(@isempty, enzymeDelta) & ~cellfun(@(x) isnumeric(x) && x==0, enzymeDelta));

subpools = subpools';

counterMap = containers.Map();

for i = 1:numel(subpools)
    currentString = subpools{i};

    if isKey(counterMap, currentString)
        counter = counterMap(currentString) + 1;
        counterMap(currentString) = counter;

        newString = sprintf('%s_%d', currentString(1:end-2), counter);

        subpools{i} = newString;

    else
        counterMap(currentString) = 1;
    end

end

deltaTable = table();
deltaTable.Subpools = subpools;
deltaTable.UsageWT = enzymeWT;
deltaTable.UsageKO = enzymeKO;
deltaTable.Deltas = enzymeDelta;

deltaTable = deltaTable(deltaTable.Deltas ~= 0, :);

writetable(deltaTable, filename, 'Delimiter','\t')

%% Get subpool delta for all subpools

% subpools = {};
% enzymeDelta = [];
% 
% for e = 1:length(KOresultsE1used)
%     currentEnzKO = KOresultsE1used{e,4}(9345:16604,:);
%     currentEnzEKcat = SolutionEKcat.full(9345:16604,:);
% 
%     for s = 1:length(currentEnzKO)
%         enzymeDelta = [enzymeDelta currentEnzKO(s,1) - currentEnzEKcat(s,1)];
%         subpools = [subpools KOresultsE1used{e,1}];
%     end
% 
% end
% 
% enzymeDelta = enzymeDelta';
% 
% % enzymeDelta = enzymeDelta(enzymeDelta ~= 0);
% % enzymeDelta = enzymeDelta(~cellfun(@isempty, enzymeDelta) & ~cellfun(@(x) isnumeric(x) && x==0, enzymeDelta));
% 
% subpools = subpools';
% 
% counterMap = containers.Map();
% 
% for i = 1:numel(subpools)
%     currentString = subpools{i};
% 
%     if isKey(counterMap, currentString)
%         counter = counterMap(currentString) + 1;
%         counterMap(currentString) = counter;
% 
%         newString = sprintf('%s_%d', currentString(1:end-2), counter);
% 
%         subpools{i} = newString;
% 
%     else
%         counterMap(currentString) = 1;
%     end
% 
% end
% 
% deltaTable = table();
% deltaTable.Subpools = subpools;
% deltaTable.Deltas = enzymeDelta;
% 
% deltaTable = deltaTable(deltaTable.Deltas ~= 0, :);
% 
% % writetable(deltaTable, filename, 'Delimiter','\t')
% 
% KOsolutionTable = table();
% KOsolutionTable.EKcat = SolutionEKcat.full(9345:16604,:);
% 
% solMatrix = [];
% 
% for sol = 1:length(KOresultsE1used)
%     currentSol = KOresultsE1used{sol,4}(9345:16604,:);
%     solMatrix = [solMatrix currentSol];
% end
% 
% colNames = KOresultsE1used(:,1);
% solTable = array2table(solMatrix, 'VariableNames', colNames);
% 
% KOsolutionTable = [KOsolutionTable solTable];
% 
% rxnNames = model.rxnNames(9345:16604);
% 
% KOsolutionTable.Properties.RowNames = rxnNames;
% 
% KOsolutionTableDelta = KOsolutionTable;
% 
% for col = 2:width(KOsolutionTableDelta)
%     currentCol = KOsolutionTableDelta(:,col);
%     KOsolutionTableDelta(:,col) = currentCol - KOsolutionTable.EKcat;
% end
% 
% writetable(KOsolutionTableDelta, "SKOsolutionTableDelta_glucose.csv", 'Delimiter','\t')