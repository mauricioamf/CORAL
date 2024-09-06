%% Cleaning the workspace and the command window
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Loading the models and other data
% iJO1366
load('../Models/iJO1366.mat') 
iJO1366 = model;
clear model

% iJO1366 with underground reactions
load('../Models/iJO1366_underground.mat')                  
iJO1366_und = model;
clear model

% iML1515
load('../Models/iML1515.mat')         
iML1515 = model;
clear model

% iML1515 with adjustments
load('../Models/iML1515_adjusted.mat')
iML1515_adj = model;
clear model

%% Display the number of metabolites, reactions, and genes in each model

models = ["iJO1366", "iJO1366_und", "iML1515", "iML1515_adj"];

metabolites = zeros(4,1);
reactions = zeros(4,1);
genes = zeros(4,1);
connectivities = zeros(4,1);
uniquePathways = zeros(4,1);

for i = 1:numel(models)
    metabolites(i,1) = length(eval(models{i}).mets);
    reactions(i,1) = length(eval(models{i}).rxns);
    genes(i,1) = length(eval(models{i}).genes);
    connectivities(i,1) = sum(sum(eval(models{i}).S ~= 0));
end

models_table = rows2vars(table(metabolites, reactions, genes, connectivities, 'RowNames', models'));
models_table.Properties.VariableNames(1) = "ModelStructure";
models_table.ModelStructure = char(models_table.ModelStructure);
models_table.ModelStructure(1,1:11) = "Metabolites";
models_table.ModelStructure(2,1:9) = "Reactions";
models_table.ModelStructure(3,1:5) = "Genes";
models_table.ModelStructure(4,1:14) = "Connectivities";

disp(models_table)

%% Get the set of unique reactions in each model
% Between iML1515 versions
newiML1515_adj  = setdiff(iML1515_adj.rxns, iML1515.rxns);
removed_iML1515 = setdiff(iML1515.rxns, iML1515_adj.rxns);

% Between iML1515 and iJO1366
uniqueToiJO1366 = setdiff(iJO1366.rxns, iML1515.rxns);
uniqueToiML1515 = setdiff(iML1515.rxns, iJO1366.rxns);

% Between iML1515_adj and iJO1366_und
uniqueToiJO1366_und = setdiff(iJO1366_und.rxns, iML1515_adj.rxns);
uniqueToiML1515_adj = setdiff(iML1515_adj.rxns, iJO1366_und.rxns);