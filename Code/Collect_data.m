%
% LP
%
%% Cleaning the workspace and the command window
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Loading the enzyme-constrained model and other data
% Conventional ecmodel
% load('../Models/eciML1515_stage3_BRENDA.mat');

% Underground ecmodel
model = readYAMLmodel('../Models/eciML1515_underground_stage2_BRENDA.yml');
filename = 'eciML1515_underground_stage2_BRENDA.csv';

%% Collect data from the model

EcKcat = model.ec.kcat;
EcRxns = model.ec.rxns;

EcRxnNames{length(EcRxns),1} = {};

for i = 1:length(EcRxns)
    rxnIndex = findRxnIDs(model, EcRxns{i});
    if ~isempty(rxnIndex)
        EcRxnNames{i} = model.rxnNames{rxnIndex};
    else
        EcRxnNames{i} = 'Reaction ID not found';
    end
end

EcRxnIDsAll = findRxnIDs(model, EcRxns);

model = buildRxnGeneMat(model);

EcGene{length(EcRxnIDsAll),1} = {};

for i = 1:length(EcRxnIDsAll)
    EcGene{i,1} = model.genes(find(model.rxnGeneMat(EcRxnIDsAll(i),:)));
end

EcRxnIDsEc = transpose(1:length(EcRxns));

EcEnz{length(EcRxnIDsEc),1} = {};

for i = 1:length(EcRxnIDsEc)
    EcEnz{i,1} = model.ec.enzymes(find(model.ec.rxnEnzMat(EcRxnIDsEc(i),:)));
end

EcTable = {};
EcTable(:,1) = num2cell(EcRxnIDsAll);
EcTable(:,2) = num2cell(EcRxnIDsEc);
EcTable(:,3) = EcRxns;
EcTable(:,4) = EcRxnNames;
EcTable(:,5) = EcGene;
EcTable(:,6) = EcEnz;
EcTable(:,7) = num2cell(EcKcat);
EcTable = cell2table(EcTable);
EcTable.Properties.VariableNames = {'ModelRxnIDs' ...
                                    'EcRxnIDs' ...
                                    'Reactions' ...
                                    'ReactionNames' ...
                                    'Genes' ...
                                    'Enzymes' ...
                                    'Kcats'};

mask = cellfun(@(x) numel(x) > 1, EcTable.Genes);
EcTable = EcTable(~mask, :);

mask2 = strcmp(EcTable.Enzymes, 'standard');
EcTable = EcTable(~mask2, :);

for i = 1:numel(EcTable.Genes)
    currentCell = EcTable.Genes{i};
    EcTable.Genes{i} = strjoin(currentCell, ' ');
    EcTable.Enzymes{i} = strjoin(currentCell, ' ');
end

%%
writetable(EcTable, filename, 'Delimiter','\t');