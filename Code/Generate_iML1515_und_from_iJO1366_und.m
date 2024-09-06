%% Merge models
%% Load models
load("../Models/iML1515.mat");
iML1515  = model;
clear model

load("../Models/iJO1366.mat");
iJO = model;
clear model

load("../Models/iJO1366_underground.mat");
iJO_underground = model;
clear model

%% Discard reactions already in iML1515
iJO = ravenCobraWrapper(iJO);
iJO_underground = ravenCobraWrapper(iJO_underground);
iML1515 = ravenCobraWrapper(iML1515);

[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(iML1515);
iML1515.grRules = grRules;

clear grRules rxnGeneMat indexes2check

[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(iJO);
iJO.grRules = grRules;

clear grRules rxnGeneMat indexes2check

% From this list, through manual curation define the following new grRules
fid         = fopen('../Models/change_grRules.txt');
loadedData  = textscan(fid,'%s %s','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = loadedData{2};

iJO = changeGrRules(iJO,rxns,grRules,true);

[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(iJO);
iJO.grRules = grRules;

clear grRules rxnGeneMat indexes2check

[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(iJO_underground);
iJO_underground.grRules = grRules;

clear grRules rxnGeneMat indexes2check

%From this list, through manual curation define the following new grRules
fid         = fopen('../Models/change_grRules.txt');
loadedData  = textscan(fid,'%s %s','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = loadedData{2};

iJO_underground = changeGrRules(iJO_underground,rxns,grRules,true);

[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(iJO_underground);
iJO_underground.grRules = grRules;

clear grRules rxnGeneMat indexes2check

%% Filter underground reactions
iJO_underground = removeReactions(iJO_underground,contains(iJO_underground.rxns,iJO.rxns),true,true,true);
 
disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(iJO_underground.genes)) ' / ' ...
    num2str(length(iJO_underground.rxns)) ' / ' ...
    num2str(length(iJO_underground.mets))])

iJO_underground = removeReactions(iJO_underground,contains(iJO_underground.rxns,iML1515.rxns),true,true,true);
 
disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(iJO_underground.genes)) ' / ' ...
    num2str(length(iJO_underground.rxns)) ' / ' ...
    num2str(length(iJO_underground.mets))])

toExclude = iJO_underground.rxns(find(~startsWith(iJO_underground.rxns, 'u0')));

iJO_underground = removeReactions(iJO_underground,toExclude,true,true,true);

%% Add underground reactions to iML1515
modelComb = mergeModels({iML1515,iJO_underground},"metNames",true);
% iML1515_underground = contractModel(modelComb);
iML1515_underground = modelComb;

[Unq, isUnique] = checkCobraModelUnique(iML1515_underground);

[modelOut, removedRxnInd, keptRxnInd] = checkDuplicateRxn(iML1515_underground);

iML1515_underground = modelOut;

iML1515_underground.subSystems = getSubSystem_str(iML1515_underground);

% iML1515_underground = removeMetabolites(iML1515_underground,'murein5px4p_p_iRN1366u');

% iML1515_underground = removeGenes(iML1515_underground,'b0650');

% iML1515_underground = ravenCobraWrapper(iML1515_underground);

writeYAMLmodel(iML1515_underground,'iML1515u.yml');

model = readYAMLmodel('iML1515_underground.yml');

checkModelStruct(model)

save('iML1515u_v2.mat','iML1515_underground');
%load('iML1515_underground.mat','iML1515_underground');

exportToExcelFormat(iML1515_underground,'iML1515u_v2.xlsx');


%% Compare models

load("iML1515.mat");
iML1515  = model;
clear model

load("iJO1366.mat");
iJO1366 = model;
clear model

load("iJO1366_underground.mat");
iJO1366_und = model;
clear model

load("iML1515_underground.mat");
iML1515_und = iML1515_underground;
clear iML1515_underground

% Display the number of metabolites, reactions, and genes in each model

models = ["iJO1366", "iJO1366_und", "iML1515", "iML1515_und"];

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


%% Test iML1515_und
%% Checking fluxes without adjusting parameters 

model = iML1515_und;

solution = optimizeCbModel(model, 'max');
try
    printFluxes(model, solution.x, true);
catch
    disp("No feasible solution")
end
clear solution

%%