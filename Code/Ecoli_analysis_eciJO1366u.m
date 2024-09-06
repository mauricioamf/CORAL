%% Cleaning the workspace and the command window and adjusting settings
clear;clc
changeCobraSolver('gurobi', 'LP');

%% Predict usage of underground reations
load('../Models/eciJO1366u_BRENDA.mat');
% load('../Models/eciJO1366u_DLKcat.mat');

%% Define constraints and other parameters
biomassRxnID = 'Ec_biomass_iJO1366_core_53p95M';
vBio = 0.06;

% M9 media without carbon
exchangeIndex = find(contains(model.rxnNames, "exchange"));
exchangeIDs = model.rxns(exchangeIndex);
exchangeIDs(end) = [];

model = changeRxnBounds(model, exchangeIDs, 0, 'l');

M9_components = ["EX_pi_e", "EX_co2_e", "EX_fe3_e", "EX_h_e", ...
    "EX_mn2_e", "EX_fe2_e", "EX_zn2_e", "EX_mg2_e", ...
    "EX_ca2_e", "EX_ni2_e", "EX_cu2_e", "EX_sel_e", ...
    "EX_cobalt2_e", "EX_h2o_e", "EX_mobd_e", "EX_so4_e", ...
    "EX_nh4_e", "EX_k_e", "EX_na1_e", "EX_cl_e", ...
    "EX_o2_e", "EX_tungs_e", "EX_slnt_e"];

model = changeRxnBounds(model, M9_components, -1000, 'l');

% Anaerobic growth
% model = changeRxnBounds(model, 'EX_o2_e', 0, 'l');

% Carbon sources
model = changeRxnBounds(model, 'EX_glc_e', -10, 'l'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', -10, 'l'); %glycerol
% model = changeRxnBounds(model, 'EX_xyl_D_e', -10, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc_L_e', -10, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab_L_e', -10, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_ac_e', -10, 'l'); %acetate

% Desired product
% model = changeRxnBounds(model, 'EX_alltn_e', 0.01, 'l'); %Adenosylcobalamin

%% How many underground reaction does pFBA/EKcat predicts?
SolutionEKcat = solveKcatE(model, biomassRxnID, vBio);

%% Solve CORAL LP
[~, SolutionTheta] = solveCoralLP(model, biomassRxnID, vBio, 1, 1);
