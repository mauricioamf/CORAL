%
% Optimises Z = v_bio - lambda * sum(theta_i)
%
%% Cleaning the workspace and the command window
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Loading the enzyme-constrained model and other data
% Underground ecmodel
model = readYAMLmodel('../Models/eciML1515_underground_stage2_DLKcat.yml');
filename = 'KcatRatios_DLKcat.csv';

%% Retrieve kcat values
rxns = model.ec.rxns;
updateRxns = true(numel(rxns),1);
KcatAll=zeros(numel(updateRxns)*10,5);
updateRxns=find(updateRxns);
kcatFirst=0;

for i=1:numel(updateRxns)
    j=updateRxns(i);
    enzymes   = find(model.ec.rxnEnzMat(j,:));
    kcatLast  = kcatFirst+numel(enzymes);
    kcatFirst = kcatFirst+1;
    KcatAll(kcatFirst:kcatLast,1) = j;
    KcatAll(kcatFirst:kcatLast,2) = enzymes;
    KcatAll(kcatFirst:kcatLast,3) = model.ec.rxnEnzMat(j,enzymes);
    KcatAll(kcatFirst:kcatLast,4) = model.ec.kcat(j);
    KcatAll(kcatFirst:kcatLast,5) = model.ec.mw(enzymes);
    kcatFirst = kcatLast;
end

KcatAll(kcatLast+1:end,:)=[];

%%
% Column index for grouping
groupColumnIndex = 2; % Grouping based on the second column

% Get unique values from the grouping column
uniqueValues = unique(KcatAll(:, groupColumnIndex));

% Create a struct to hold grouped data
EnzymeAll = struct();

% Iterate through unique values
for i = 1:length(uniqueValues)
    value = uniqueValues(i);
    
    % Find rows with the same value in the grouping column
    matchingRows = KcatAll(KcatAll(:, groupColumnIndex) == value, :);
    
    % Create a field in the struct
    fieldname = ['Enzyme_' num2str(value)];
    EnzymeAll.(fieldname) = matchingRows;
end

%%
enzymes = fieldnames(EnzymeAll);

KcatRatio = {};

for i = 1:numel(enzymes)
    if ( isnumeric(EnzymeAll.(enzymes{i})) )
        currentField = EnzymeAll.(enzymes{i});
        currentField = sortrows(currentField,4,"descend");
        
        numEntries = size(currentField, 1);
        ratios = zeros(numEntries - 1, 1);
        
        for j = 1:(numEntries - 1)
            kcat_n = currentField(j,4);
            kcat_np = currentField(j+1,4);
            ratio_nnp = kcat_n / kcat_np;
            ratios(j) = ratio_nnp;
        end
    
        KcatRatio{i,1} = ratios;
    end
end

%%
writecell(KcatRatio, filename, 'Delimiter','\t');