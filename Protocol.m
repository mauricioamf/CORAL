%% Protocol for generating a pcGEM with the CORAL Toolbox formulation. 
%
% This script contain step-by-step instructions on how to run CORAL for
% your pcGEM to restructure its enzyme usage. It replicates the Escherichia
% coli pcGEM used for all analyses in the accompanying manuscript. Please
% adjust all variables, paths and data for your organism as needed.
% 
% This restructuring was tested only with a GECKO 3 version of the 
% eciML1515 model. Regardless, it should work for any pcGEM reconstructed 
% using GECKO 3 up to version 3.1.3. Please open an Issue at the GitHub
% page in case it does not work for your pcGEM.
%
% Your pcGEM might require curation and evaluation before and/or after the
% restructuring using CORAL. 
%
% The CORAL Toolbox might receive updates after the manuscript is 
% published. Refer to version 1.0.0 if you are looking into reproducing the
% results in the manuscript.
%
% The functions used for restructuring the pcGEM are:
%
%   - expandEnzComplex.m
%   - getEnzymeSubpools.m
%   - CreateUndField.m

%% Preparation for the pcGEM restructuring
% Make sure all dependencies are installed and working properly.
% Initialize COBRA
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');

% Check that RAVEN and GECKO are working.
checkInstallation;
GECKOInstaller.install

% Load the enzyme-constrained model.
load('./Models/eciML1515u_v2_stage2_DLKcat.mat');

% Load the modelAdapter. It should be the same as used for reconstructing
% the pcGEM. Use the absolute path of the modelAdapter directory.
adapterLocation = '/absolute_path/CORAL/GECKO3/iML1515u/code/iML1515uAdapter.m';
% Open the modelAdapter file and adjust the path on the obj.params.path 
% variable before proceeding.
modelAdapter = ModelAdapterManager.setDefault(adapterLocation);

%% STAGE 1: Simplify GPR rules
% STEP 1: Change directory to Code folder.
currentPath = pwd;
cd ./Code/

% STEP 2: Split reactions catalysed by enzyme complexes into multiple 
% reactions. This step does not require additional inputs from the user.
model = expandEnzComplex(model);

% (optional) STEP 2.1: The eciML1515 model requires manual corrections to GPR rules to proceed 
% to the next stage. Make sure your model also has the correct format for
% GPR rules and other fields.
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

% STEP 3: Save checkpoint.
save('CORAL_checkpoint1');
% load('CORAL_checkpoint1.mat');

%% STAGE 2: Restructure the model
% STEP 4: Include enzyme usage of promiscuous enzymes as enzyme subpools.
% This step does not require additional inputs from the user.
model = getEnzymeSubpools(model);

% STEP 5: Save checkpoint.
save('CORAL_checkpoint2');
% load('CORAL_checkpoint2.mat');

% STEP 6: Organize new information into the new .und field
model = CreateUndField(model, modelAdapter);

%% STAGE 3: Export model
% STEP 7: Export to YAML format.
writeYAMLmodel(model, '../Models/eciML1515u_CORAL.yml');
% STEP 8: Export to MAT format.
save('../Models/eciML1515u_CORAL.mat','model');

%% Return to starting directory
cd(currentPath)