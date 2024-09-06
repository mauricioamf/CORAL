%% Cleaning the workspace and the command window and adjusting settings
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
% initCobraToolbox(false);
clear;clc

%% Predict usage of underground reations
load('../Models/iML1515u_v2.mat');
model = iML1515_underground;
clear iML1515_underground
model = ravenCobraWrapper(model);

%% Define constraints and other parameters
biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
biomassID = find(contains(model.rxns, biomassRxnID));
maxBio = 0.6;
minBio = 0.1;

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
model = changeRxnBounds(model, 'EX_glc__D_e', -2, 'b'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', -10, 'l'); %glycerol
% model = changeRxnBounds(model, 'EX_xyl__D_e', -10, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc__L_e', -10, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab__L_e', -10, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_fru_e', -10, 'l'); %fructose

% Desired product
% model = changeRxnBounds(model, 'EX_alltn_e', 0.05, 'l'); %Adenosylcobalamin

%% Construct and solve the problem
modelpFBA = model;
modelpFBA = buildRxnGeneMat(modelpFBA);

% Fix biomass
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, minBio, 'l');
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, maxBio, 'u');

% Set lower and upper bounds
modelpFBA.lb(modelpFBA.lb==-Inf) = -1000;
modelpFBA.ub(modelpFBA.ub==Inf) = 1000;

% Set the objective
[~,nRxns] = size(modelpFBA.S);

modelpFBA.c = zeros(nRxns,1);
% modelpFBA.c = ones(nRxns,1);
modelpFBA.c(biomassID) = 1;

modelpFBA.osense = -1;

% Solve the problem
% SolutionEKcat = solveCobraLP(modelpFBA);
SolutionEKcat = optimizeCbModel(modelpFBA);

% Get the solutions
% SolutionEKcat.x = SolutionEKcat.full;
SolutionEKcat.full = SolutionEKcat.x;

printFluxes(modelpFBA, SolutionEKcat.x, true);
% printFluxes(modelpFBA, SolutionEKcat.x, false);

%% Perform single and double KOs
% tic
% [grRatioSing, grRateKOSing, grRateWTSing, hasEffectSing, delRxnsSing, fluxSolutionSing] = singleGeneDeletion(model);
% [grRatioDble, grRateKO, grRateWT] = doubleGeneDeletion(model);
% toc
% save("DKO_GEM.mat", "-v7.3")

%%
% newCellArray = cell(0, 4);
% 
% for i = 1:1530
%     for j = 1:1530
%         % Extract data for single and double KO
%         gene1 = model.genes{i};
%         singleKORate = grRateKOSing(i);
%         gene2 = model.genes{j};
%         doubleKORate = grRateKO(i, j);
% 
%         % Create a new row for the new cell array
%         newRow = {gene1, singleKORate, gene2, doubleKORate};
% 
%         % Concatenate the new row to the existing cell array
%         newCellArray = [newCellArray; newRow];
%     end
% end
% toc

%% Knock-out a single gene
tic
KOresults = {};

% Loop through each reaction to knockout
for i = 1:length(modelpFBA.genes)
    % Clone the original model to avoid modifying it
    model_knockout = modelpFBA;
    
    % Knockout the reaction
    % disp(['Knocking out subenzyme ' nativeEnzs{i}]);
    model_knockout = deleteModelGenes(model_knockout, model_knockout.genes{i});
    
    % Run FBA
    SolutionKO = solveCobraLP(model_knockout);
    % SolutionKO = optimizeCbModel(model_knockout);
    % SolutionKO.full = SolutionKO.x;
    
    % Check if the model grows or not
    if SolutionKO.stat ~= 0
        KOresults{i, 1} = model_knockout.geneUniprotID{i};
        KOresults{i, 2} = 'Feasible';
        KOresults{i, 3} = SolutionKO.full(biomassID);
        KOresults{i, 4} = SolutionKO.full;
    else
        KOresults{i, 1} = model_knockout.geneUniprotID{i};
        KOresults{i, 2} = 'Infeasible';
    end
end
toc

%% Filter out infeasibilities
infeasible = [];
for k = 1:length(KOresults)
    if KOresults{k, 2} == "Infeasible"
        infeasible = [infeasible k];
    end
end
KOresults(infeasible,:) = [];

%% Perform double KOs based on filtered KOresults
tic
% DKOresults = cell(length(KOresultsE1used), 1);
DKOresults = {};

for ko = 766:1530

    modelDKO = modelpFBA;
    
    % modelDKO = deleteModelGenes(modelDKO, model_knockout.genes{ko});

    genePair = {};
    genePair{1,1} = modelDKO.genes{ko};

    for dko = 1:length(model.genes)
        modelDKOsub = modelDKO;

        % disp(['... subenzyme ' nativeEnzs{dko}]);
        genePair{2,1} = modelDKO.genes{dko};
        modelDKOsub = deleteModelGenes(modelDKOsub, genePair);
        
        % Run FBA
        SolutionDKO = solveCobraLP(modelDKOsub);
        % SolutionDKO = optimizeCbModel(modelDKOsub);
        % SolutionDKO.full = SolutionDKO.x;
        
        % Check if the model grows or not
        if SolutionDKO.stat ~= 0
            DKOresults{ko}{dko, 1} = model_knockout.geneUniprotID{dko};
            DKOresults{ko}{dko, 2} = 'Feasible';
            DKOresults{ko}{dko, 3} = SolutionDKO.full(biomassID);
            DKOresults{ko}{dko, 4} = SolutionDKO.full;
        else
            DKOresults{ko}{dko, 1} = model_knockout.geneUniprotID{dko};
            DKOresults{ko}{dko, 2} = 'Infeasible';
        end
    end

end

DKOresults = DKOresults';
toc

%% Split the results cell array
KOresultsQuarter1 = KOresults(1:383, :);
KOresultsQuarter2 = KOresults(384:765, :);
KOresultsQuarter3 = KOresults(766:1147, :);
KOresultsQuarter4 = KOresults(1148:length(KOresults), :);

DKOresultsQuarter1 = DKOresults(1:383, :);
DKOresultsQuarter2 = DKOresults(384:765, :);
DKOresultsQuarter3 = DKOresults(766:1147, :);
DKOresultsQuarter4 = DKOresults(1148:length(KOresults), :);

%% Export results
% pattern = '^usage_prot_(\w+)_(\d+)$';
clear ko modelDKO genePair dko modelDKOsub SolutionDKO exportTable currentFirstKO currentFirstMu currentDKOCell exportRow

exportTable = table();

for e = 1:length(KOresultsQuarter4)

    % if e == 14, continue, end

    currentFirstKO = KOresultsQuarter4{e,1};
    currentFirstMu = KOresultsQuarter4{e,3};

    % match = regexp(currentFirstKO, pattern, 'tokens', 'once');
    % id = match{1};
    % currentFirstSubpoolID = ['subpool_', id];
    % currentFirstSubpoolNumID = findRxnIDs(model_knockout, currentFirstSubpoolID);
    % currentFirstSubpoolNum = KOresultsTopE1used{e,4}(currentFirstSubpoolNumID);

    currentDKOCell = DKOresultsQuarter4{e,1};

    for d = 1:length(currentDKOCell)
        currentDKO = currentDKOCell{d,1};
        currentDMu = currentDKOCell{d,3};

        % match = regexp(currentDKO, pattern, 'tokens', 'once');
        % id = match{1};
        % currentDKOSubpoolID = ['subpool_', id];
        % currentDKOSubpoolNum = findRxnIDs(odel_knockout, currentFirstSubpoolID);
        % currentDKOSubpoolNum = currentDKOCell{d,4}(currentDKOSubpoolNum);

        exportRow{1,1} = currentFirstKO;
        exportRow{1,2} = currentFirstMu;
        exportRow{1,3} = currentDKO;
        if isempty(currentDMu) == true
            exportRow{1,4} = 0;
        else
            exportRow{1,4} = currentDMu;
        end

        exportTable = [exportTable; exportRow];

    end

end

exportTable = renamevars(exportTable, ["Var1","Var2","Var3","Var4"],...
    ["FirstKO","FirstMu","DKO","DKOMu"]);

filename = "DKO_GEM_quarter4.csv";
writetable(exportTable, filename, 'Delimiter','\t')

% save("DKO_GEM_Half2.mat", "-v7.3")

%%
toc