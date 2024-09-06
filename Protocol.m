%
% Protocol for running CORAL Toolbox. To be expanded.
%
%%
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP');
cd /home/mferreira/Doutorado/4.Underground_metabolism/underground_metabolism/Code/

%% Cleaning the workspace and the command window and adjusting settings
clear;clc

%% Loading the enzyme-constrained model and other data
% model = readYAMLmodel('../Models/eciML1515_underground_stage2_BRENDA.yml');
% load('../Models/eciML1515_underground_stage2_BRENDA.mat');

% model = readYAMLmodel('../Models/eciML1515_underground_stage2_DLKcat.yml');
% load('../Models/eciML1515_underground_stage2_DLKcat.mat');

% model = readYAMLmodel('../Models/eciML1515u_v2_stage2_BRENDA.yml');
% load('../Models/eciML1515u_v2_stage2_BRENDA.mat');

% model = readYAMLmodel('../Models/eciML1515u_v2_stage2_DLKcat.yml');
load('../Models/eciML1515u_v2_stage2_DLKcat.mat');
% load('../Models/eciML1515u_v2_stage3_DLKcat.mat');

% load('../Models/eciML1515u_TurNuP.mat');

adapterLocation = '/home/mferreira/Doutorado/4.Underground_metabolism/underground_metabolism/GECKO3/iML1515u/iML1515uAdapter.m';
% adapterLocation = '/home/mferreira/Doutorado/4.Underground_metabolism/underground_metabolism/GECKO3/iML1515u_TurNuP/code/iML1515uAdapter.m';
modelAdapter = ModelAdapterManager.setDefault(adapterLocation);

%% STEP 1 Simplify GPR rules by splitting reactions catalysed by enzyme complexes into multiple reactions
model = expandEnzComplex(model);

% Manual corrections to model
% model.grRules{965,1} = 'b0650';
model = rmfield(model, "rules");

fixedGrRules = cell(size(model.grRules));
for i = 1:numel(model.grRules)
    fixedGrRules{i} = strrep(model.grRules{i}, '( ', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '(( ', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '(', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '((', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ')', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '))', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ' )', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ' ))', '');
end
model.grRules = strtrim(fixedGrRules);

fixedGenes = cell(size(model.genes));
for i = 1:numel(model.genes)
    fixedGenes{i} = strrep(model.genes{i}, '( ', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '(( ', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '(', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '((', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ')', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '))', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ' )', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ' ))', '');
end
model.genes = strtrim(fixedGenes);

% Save checkpoint
% save('checkpoint1');
% load('checkpoint1.mat');

%% STEP 2 Restructure the model to include enzyme usage of promiscuous enzymes as enzyme subpools
model = getEnzymeSubpools(model);

% Save checkpoint
% save('checkpoint2');
% load('checkpoint2.mat');

%% STEP 3 Organize new information into the new .und field
model = CreateUndField(model, modelAdapter);

%% Test if model can produce growth
% testEcoliModels(model, 2);
% testEcoliModels_iJO1366(model, 2);

%% Export model
% writeYAMLmodel(model, '../Models/eciML1515u_BRENDA.yml');
% save('../Models/eciML1515u_BRENDA.mat','model');

% writeYAMLmodel(model, '../Models/eciML1515u_DLKcat.yml');
% save('../Models/eciML1515u_DLKcat_CORAL_test.mat','model');

save('../Models/eciML1515u_DLKcat_CORAL.mat','model');

%%
toc