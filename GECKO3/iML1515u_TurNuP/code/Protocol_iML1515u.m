%%
% initCobraToolbox(false);
cd('~/CORAL/GECKO3/iML1515u_TurNuP/code/');
clc;clear

% delete ../data/DLKcat.tsv ../data/kegg.tsv ../data/smilesDB.tsv ../data/uniprot.tsv

%% STAGE 1: Expansion from a starting metabolic model to an ecModel structure
adapterLocation = '~/CORAL/GECKO3/iML1515u_TurNuP/iML1515uAdapter.m';
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Load conventional GEM
model = loadConventionalGEM();

% Correct metabolites
correctMets = cell(size(model.mets));
correctMetNames = cell(size(model.metNames));

for i = 1:numel(model.mets)
    correctMets{i} = strtrim(model.mets{i});
    correctMetNames{i} = strtrim(model.metNames{i});
end

model.mets = correctMets;
model.metNames = correctMetNames;

% Prepare ecModel
% We will make a full GECKO ecModel
[ecModel, noUniprot] = makeEcModel(model,false);

% Annotate with complex data
complexInfo = getComplexData();
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% Fix subsystems
ecModel.subSystems = getSubSystem_str(ecModel);

% Store model in YAML format
% saveEcModel(ecModel,'eciML1515u_TurNuP_stage1.yml');

%% STAGE 2: integration of kcat into the ecModel structure
% ecModel = loadEcModel('eciML1515u_TurNuP_stage1.yml'); 

% Gather EC numbers
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
ecModel         = getECfromDatabase(ecModel);

% Gather kcat values from BRENDA
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% Gather metabolite SMILES
[ecModel, noSmiles] = findMetSmiles(ecModel);

% Prepare DLKcat input file
% writeDLKcatInput(ecModel,[],[],[],[],true);

% Run DLKcat
% runDLKcat();

% ecModel.metNames = metNames_new;

% Load DLKcat output
kcatList_DLKcat = readDLKcatOutput(ecModel);

% Combine kcat from BRENDA and DLKcat
% kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% Take kcatList and populate edModel.ec.kcat
% ecModel_BRENDA  = selectKcatValue(ecModel, kcatList_fuzzy);
ecModel  = selectKcatValue(ecModel, kcatList_DLKcat);
% ecModel_Merged  = selectKcatValue(ecModel, kcatList_merged);

%  Apply TurNuP kcat values
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% Get kcat values across isozymes
% ecModel_BRENDA = getKcatAcrossIsozymes(ecModel_BRENDA);
ecModel = getKcatAcrossIsozymes(ecModel);
% ecModel_Merged = getKcatAcrossIsozymes(ecModel_Merged);

% Get standard kcat
% ecModel_BRENDA.subSystems = getSubSystem_str(ecModel_BRENDA);
ecModel.subSystems = getSubSystem_str(ecModel);
% ecModel_Merged.subSystems = getSubSystem_str(ecModel_Merged);

% [ecModel_BRENDA, rxnsMissingGPR_BRENDA, standardMW_BRENDA, standardKcat_BRENDA] = getStandardKcat(ecModel_BRENDA);
[ecModel, rxnsMissingGPR_DLKcat, standardMW_DLKcat, standardKcat_DLKcat] = getStandardKcat(ecModel);
% [ecModel_Merged, rxnsMissingGPR_Merged, standardMW_Merged, standardKcat_Merged] = getStandardKcat(ecModel_Merged);

% Apply kcat constraints from ecModel.ec.kcat to ecModel.S
% ecModel_BRENDA = applyKcatConstraints(ecModel_BRENDA);
ecModel = applyKcatConstraints(ecModel);
% ecModel_Merged = applyKcatConstraints(ecModel_Merged);

% Set upper bound of protein pool
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

% ecModel_BRENDA = setProtPoolSize(ecModel_BRENDA,Ptot,f,sigma);
ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);
% ecModel_Merged = setProtPoolSize(ecModel_Merged,Ptot,f,sigma);

% Store model in YAML format
% saveEcModel(ecModel_BRENDA,'eciML1515u_v2_stage2_BRENDA.yml');
saveEcModel(ecModel,'eciML1515u_TurNuP_stage2.yml');
% saveEcModel(ecModel_Merged,'eciML1515u_v2_stage2_Merged.yml');

%% STAGE 3: model tuning
% ecModel = loadEcModel('eciML1515u_TurNuP_stage2.yml', ModelAdapter); 

% Test maximum growth rate
ecModel = setParam(ecModel,'lb','EX_glc__D_e',-1000);

% And set growth maximization as the objective function.
ecModel = setParam(ecModel,'obj','BIOMASS_Ec_iML1515_core_75p37M',1);

% Run FBA.
sol = solveLP(ecModel,1);

fprintf('Growth rate that is reached: %f h-1\n', abs(sol.f));
%printFluxes(ecModel,sol.x,true);

% Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

% Inspect the tunedKcats structure in table format.
struct2table(tunedKcats)

% This functional ecModel will also be kept in the GECKO GitHub.
saveEcModel(ecModel,'eciML1515u_TurNuP_stage3.yml');

saveEcModel(ecModel,'eciML1515u_TurNuP.yml');