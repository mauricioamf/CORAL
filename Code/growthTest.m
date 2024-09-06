%% Cleaning the workspace and the command window and adjusting settings
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Predict usage of underground reations
% load('../Models/eciML1515u_DLKcat_CORAL.mat');
% load('../Models/eciML1515u_CORAL_TurNuP.mat');

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'b'); %D-glucose

% model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'b'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', -1.5, 'l'); %glycerol
% model = changeRxnBounds(model, 'EX_xyl__D_e', -1.5, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc__L_e', -1, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab__L_e', -1.5, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_fru_e', -1.4, 'l'); %fructose

%%
solution = optimizeCbModel(model);

printFluxes(model, solution.x, true);
% printFluxes(model, solution.x, false);