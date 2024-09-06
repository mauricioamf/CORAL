%
% Minimizes Etot - sum Ei*MW
%
%% Cleaning the workspace and the command window
clear;clc

%% Configurations for MATLAB
% initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
%run("~/Softwares/RAVEN/installation/checkInstallation.m");
warning("off");
tic
%% FVA for models without underground reactions

fprintf('\n' + "Starting FVA for models without underground reactions");

% load('../Models/eciML1515u_CORAL_DLKcat.mat');
load('../Models/eciML1515u_TurNuP_CORAL.mat');

model = buildRxnGeneMat(model);
enzymeTable = getEnzymeTable(model);
model = applyRatioConstraints(model, enzymeTable);

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'l'); %D-glucose

% Block underground reactions
undIndex = find(contains(model.rxns, 'u0'));
undIDs = model.rxns(undIndex);

model = changeRxnBounds(model, undIDs, 0, 'b');

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model.rxns,'usage_prot_'));
toRemove1 = model.rxnNames(toRemove1);
toRemove2 = find(contains(model.rxns,'pool_'));
toRemove2 = model.rxnNames(toRemove2);

rxnNamesList = model.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);
rxnNamesList = setdiff(rxnNamesList, toRemove2);

[minFlux, maxFlux] = fluxVariability(model, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'enzymes' 'ranges' 'minUsage' 'maxUsage'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_eciML1515u_CORAL_TurNuP_noBio_noUnd.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');
toc
%% FVA for models with underground reactions

clear

fprintf('\n' + "Starting FVA for models with underground reactions");

% load('../Models/eciML1515u_CORAL_DLKcat.mat');
load('../Models/eciML1515u_TurNuP_CORAL.mat');

model = buildRxnGeneMat(model);
enzymeTable = getEnzymeTable(model);
model = applyRatioConstraints(model, enzymeTable);

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'l'); %D-glucose

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model.rxns,'usage_prot_'));
toRemove1 = model.rxnNames(toRemove1);
toRemove2 = find(contains(model.rxns,'pool_'));
toRemove2 = model.rxnNames(toRemove2);

rxnNamesList = model.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);
rxnNamesList = setdiff(rxnNamesList, toRemove2);

[minFlux, maxFlux] = fluxVariability(model, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'enzymes' 'ranges' 'minUsage' 'maxUsage'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_eciML1515u_CORAL_TurNuP_noBio_Und_Ratio.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');

%% FVA for models without underground reactions and fixing biomass

clear

fprintf('\n' + "Starting FVA for models without underground reactions and fixing biomass");

% load('../Models/eciML1515u_CORAL_DLKcat.mat');
load('../Models/eciML1515u_TurNuP_CORAL.mat');

model = buildRxnGeneMat(model);
enzymeTable = getEnzymeTable(model);
model = applyRatioConstraints(model, enzymeTable);

% Constraints and other parameters
biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
vBio = 0.09;

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'l'); %D-glucose
model = changeRxnBounds(model, biomassRxnID, vBio, 'b');

% Block underground reactions
undIndex = find(contains(model.rxns, 'u0'));
undIDs = model.rxns(undIndex);

model = changeRxnBounds(model, undIDs, 0, 'b');

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model.rxns,'usage_prot_'));
toRemove1 = model.rxnNames(toRemove1);
toRemove2 = find(contains(model.rxns,'pool_'));
toRemove2 = model.rxnNames(toRemove2);

rxnNamesList = model.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);
rxnNamesList = setdiff(rxnNamesList, toRemove2);

[minFlux, maxFlux] = fluxVariability(model, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'enzymes' 'ranges' 'minUsage' 'maxUsage'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_eciML1515u_CORAL_TurNuP_Bio_noUnd_Ratio.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');

%% FVA for models with underground reactions and fixing biomass

clear

fprintf('\n' + "Starting FVA for models with underground reactions and fixing biomass");

% load('../Models/eciML1515u_CORAL_DLKcat.mat');
load('../Models/eciML1515u_TurNuP_CORAL.mat');

model = buildRxnGeneMat(model);
enzymeTable = getEnzymeTable(model);
model = applyRatioConstraints(model, enzymeTable);

% Constraints and other parameters
biomassRxnID = 'BIOMASS_Ec_iML1515_core_75p37M';
vBio = 0.09;

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
model = changeRxnBounds(model, 'EX_glc__D_e', -1.5, 'l'); %D-glucose
model = changeRxnBounds(model, biomassRxnID, vBio, 'b');

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model.rxns,'usage_prot_'));
toRemove1 = model.rxnNames(toRemove1);
toRemove2 = find(contains(model.rxns,'pool_'));
toRemove2 = model.rxnNames(toRemove2);

rxnNamesList = model.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);
rxnNamesList = setdiff(rxnNamesList, toRemove2);

[minFlux, maxFlux] = fluxVariability(model, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'enzymes' 'ranges' 'minUsage' 'maxUsage'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_eciML1515u_CORAL_TurNuP_Bio_Und_Ratio.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');

%%
toc
