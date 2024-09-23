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
%   - fixGPRrules.m
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

% (optional) STEP 2: The eciML1515 model requires corrections to GPR 
% rules to proceed to the next stage. Make sure your model has the 
% correct format for GPR rules and other fields. This includes but is not 
% limited to leading and trailing spaces, non-conventional characters, 
% parenthesis or brackets. 
% 
% Potential sources for errors also include duplicated genes, misspelled
% genes, fields with inconsistent sizes, and non-conventional formatting of 
% GPR rules.
% 
% The function fixGPRrules does not extensively check for irregularities in 
% the model.grRules field and some manual inspections might be necessary.
model = fixGPRrules(model);

% STEP 3: Split reactions catalysed by enzyme complexes into multiple 
% reactions. This step does not require additional inputs from the user.
model = expandEnzComplex(model);

% STEP 4: Save checkpoint.
save('CORAL_checkpoint1');
% load('CORAL_checkpoint1.mat');

%% STAGE 2: Restructure the model
% STEP 5: Include enzyme usage of promiscuous enzymes as enzyme subpools.
% This step does not require additional inputs from the user.
model = getEnzymeSubpools(model);

% STEP 6: Save checkpoint.
save('CORAL_checkpoint2');
% load('CORAL_checkpoint2.mat');

% STEP 7: Organize new information into the new .und field
model = CreateUndField(model, modelAdapter);

%% STAGE 3: Export model
% STEP 8: Export to YAML format.
writeYAMLmodel(model, '../Models/eciML1515u_CORAL.yml');
% STEP 9: Export to MAT format.
save('../Models/eciML1515u_CORAL.mat','model');

%% Return to starting directory
cd(currentPath)

%% STAGE 4: What's next?
% Your CORAL-restructured pcGEM should be ready to use. It should be
% compatible with all COBRA and RAVEN functions, but not with GECKO 3
% functions, since they depend on the exact match between the conventional
% model fields (.rxns, .mets, etc) and the model.ec fields.
