%% Cleaning the workspace and the command window and adjusting settings
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
clear;clc

%% Predict usage of underground reations
% load('../Models/eciML1515u_v2_stage2_DLKcat.mat');

% load('../Models/eciML1515u_CORAL_BRENDA.mat');
% load('../Models/eciML1515u_CORAL_DLKcat.mat');

load('../Models/eciML1515u_TurNuP_CORAL.mat');

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'b'); %D-glucose
vBio = 0.11;

% model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'b'); %D-glucose

% model = changeRxnBounds(model, 'EX_glyc_e', -1.5, 'l'); %glycerol
% vBio = 0.08;

% model = changeRxnBounds(model, 'EX_xyl__D_e', -1.5, 'l'); %D-xylose
% vBio = 0.09;

% model = changeRxnBounds(model, 'EX_fuc__L_e', -1, 'l'); %L-fucose
% vBio = 0.07;

% model = changeRxnBounds(model, 'EX_arab__L_e', -1.5, 'l'); %L-arabinose
% vBio = 0.09;

% model = changeRxnBounds(model, 'EX_fru_e', -1.4, 'l'); %fructose
% vBio = 0.1;

%% Construct and solve the problem
modelpFBA = model;
modelpFBA = buildRxnGeneMat(modelpFBA);
modelpFBA = applyRatioConstraints(modelpFBA, enzymeTable);

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));
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
modelpFBA.c(subpoolIds) = 1;
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

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));
poolID = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'prot_pool_exchange')));

modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID)*(1-0.9), 'l');
modelKO = changeRxnBounds(modelKO, biomassRxnID, SolutionEKcat.x(biomassID), 'u');
 
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
% 
% % Filter out knockout results with size 2x2
% tofiltersize = [];
% size22 = [2,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size22
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];
% 
% % Filter out knockout results with size 3x2
% tofiltersize = [];
% size32 = [3,2];
% for k = 1:length(KOresultsE1used)
%     if size(KOresultsE1used{k, 5}) == size32
%         tofiltersize = [tofiltersize k];
%     end
% end
% KOresultsE1used(tofiltersize,:) = [];
% 
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

% Filter out zero growth mutants
% zeroGrowth = [];
% for k = 1:length(KOresultsE1used)
%     if KOresultsE1used{k, 3} == 0
%         zeroGrowth = [zeroGrowth k];
%     end
% end
% KOresultsE1used(zeroGrowth,:) = [];


%% Perform double KOs based on filtered KOresults
tic
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
toc

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

% Filter out infeasible results
for k = 1:length(DKOresultsE1used)
    
    infeasible = [];
    for dk = 1:length(DKOresultsE1used{k,1})

        if DKOresultsE1used{k,1}{dk,2} == "Infeasible"
            infeasible = [infeasible dk];
        end
        
    end
    DKOresultsE1used{k,1}(infeasible,:) = [];
end

% % Filter out knockout results with size 2x1
% for k = 1:length(DKOresultsE1used)
% 
%     tofiltersize = [];
%     size12 = [1,2];
% 
%     for dk = 1:length(DKOresultsE1used{k,1})
% 
%         if size(DKOresultsE1used{k,1}{dk,5}) == size12
%             tofiltersize = [tofiltersize dk];
%         end
% 
%     end
%     DKOresultsE1used{k,1}(tofiltersize,:) = [];
% 
% end
% 
% % Filter out knockout results with size 2x2
% for k = 1:length(DKOresultsE1used)
% 
%     tofiltersize = [];
%     size22 = [2,2];
% 
%     for dk = 1:length(DKOresultsE1used{k,1})
% 
%         if size(DKOresultsE1used{k,1}{dk,5}) == size22
%             tofiltersize = [tofiltersize dk];
%         end
% 
%     end
%     DKOresultsE1used{k,1}(tofiltersize,:) = [];
% 
% end
% 
% % Filter out knockout results with size 3x2
% for k = 1:length(DKOresultsE1used)
% 
%     tofiltersize = [];
%     size32 = [3,2];
% 
%     for dk = 1:length(DKOresultsE1used{k,1})
% 
%         if size(DKOresultsE1used{k,1}{dk,5}) == size32
%             tofiltersize = [tofiltersize dk];
%         end
% 
%     end
%     DKOresultsE1used{k,1}(tofiltersize,:) = [];
% 
% end
% 
% % Filter out knockout results with size 4x2
% for k = 1:length(DKOresultsE1used)
% 
%     tofiltersize = [];
%     size42 = [4,2];
% 
%     for dk = 1:length(DKOresultsE1used{k,1})
% 
%         if size(DKOresultsE1used{k,1}{dk,5}) == size42
%             tofiltersize = [tofiltersize dk];
%         end
% 
%     end
%     DKOresultsE1used{k,1}(tofiltersize,:) = [];
% 
% end

% Only keep desired solutions
for k = 1:length(DKOresultsE1used)

    for dk = 1:length(DKOresultsE1used{k,1})
        if DKOresultsE1used{k,1}{dk,5}{1,1} == 0
            DKOresultsE1used{k,1}{dk,5} = [];
        end
    end

    notUsed = any(cellfun(@isempty, DKOresultsE1used{k,1}), 2);
    DKOresultsE1used{k,1}(notUsed,:) = [];

end

% notUsed = any(cellfun(@isempty, KOresultsE1used), 2);
% KOresultsE1used(notUsed,:) = [];

%% Export the results
exportTable = table();
pattern = '^usage_prot_(\w+)_(\d+)$';

for e = 1:length(KOresultsE1used)
    currentFirstKO = KOresultsE1used{e,1};
    currentFirstMu = KOresultsE1used{e,3};

    match = regexp(currentFirstKO, pattern, 'tokens', 'once');
    id = match{1};
    currentFirstSubpoolID = ['subpool_', id];
    currentFirstSubpoolNumID = findRxnIDs(model_knockout, currentFirstSubpoolID);
    currentFirstSubpoolNum = KOresultsE1used{e,4}(currentFirstSubpoolNumID);

    currentDKOCell = DKOresultsE1used{e,1};

    for d = 1:length(currentDKOCell)
        currentDKO = currentDKOCell{d,1};
        currentDMu = currentDKOCell{d,3};

        match = regexp(currentDKO, pattern, 'tokens', 'once');
        id = match{1};
        currentDKOSubpoolID = ['subpool_', id];
        currentDKOSubpoolNum = findRxnIDs(model_knockout, currentFirstSubpoolID);
        currentDKOSubpoolNum = currentDKOCell{d,4}(currentDKOSubpoolNum);

        exportRow{1,1} = currentFirstKO;
        exportRow{1,2} = currentFirstMu;
        exportRow{1,3} = currentFirstSubpoolID;
        exportRow{1,4} = abs(currentFirstSubpoolNum);
        exportRow{1,5} = abs(SolutionEKcat.full(currentFirstSubpoolNumID));
        exportRow{1,6} = currentDKO;
        exportRow{1,7} = currentDMu;
        exportRow{1,8} = currentDKOSubpoolID;
        exportRow{1,9} = abs(currentDKOSubpoolNum);

        exportTable = [exportTable; exportRow];

    end

end

exportTable = renamevars(exportTable, ["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9"],...
    ["FirstKO","FirstMu","FirstSubpoolID","FirstSubpoolNum","EsKcatNum","DKO","DKOMu","DKOSubpoolID","DKOSubpoolNum"]);

filename = "DKO_ratio_glucose_E1_TurNuP.csv";
writetable(exportTable, filename, 'Delimiter','\t')
% save("DKO_ratio_arabinose_E1.mat", "-v7.3")
toc

msgbox('Your code has finished running!', 'Notification');