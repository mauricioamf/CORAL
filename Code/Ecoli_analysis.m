%% Cleaning the workspace and the command window and adjusting settings
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Predict usage of underground reations
load('../Models/eciML1515u_CORAL_DLKcat_Ratio.mat');

%% Define constraints and other parameters
% enzymeTable = getEnzymeTable(model);

biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
biomassID = find(contains(model.rxns, biomassRxnID));
vBio = 0.23;

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
% model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'b'); %D-glucose

model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'b'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', -1000, 'l'); %glycerol
model = changeRxnBounds(model, 'EX_xyl__D_e', -1000, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc__L_e', -1000, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab__L_e', -1000, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_fru_e', -10, 'l'); %fructose

% Desired product
% model = changeRxnBounds(model, 'EX_alltn_e', 0.05, 'l'); %Adenosylcobalamin

%% How many underground reaction does pFBA/EKcat predicts?
% SolutionEKcat = solveKcatE(model, biomassRxnID, vBio);

%% Solve CORAL LP
% [~, SolutionTheta] = solveCoralLP(model, biomassRxnID, vBio, 1, 1);

%% Construct and solve the problem
modelpFBA = model;
modelpFBA = buildRxnGeneMat(modelpFBA);
% modelpFBA = applyRatioConstraints(modelpFBA, enzymeTable);

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));
enzUsageIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'usage_prot_')));
enzUsageIds(end) = [];

% Fix biomass
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio, 'l');
% modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio, 'u');

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
% SolutionEKcat = solveCobraLP(modelpFBA);
SolutionEKcat = optimizeCbModel(modelpFBA);

% Get the solutions
% SolutionEKcat.x = SolutionEKcat.full;
% SolutionEKcat.full = SolutionEKcat.x;

printFluxes(modelpFBA, SolutionEKcat.x, true);
% printFluxes(modelpFBA, SolutionEKcat.x, false);