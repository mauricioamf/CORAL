function combinedTable = getTable(ecModel, rxnsInd, substrateMets, productMets, substrateCoeffs, productCoeffs)
% getTable
%
%   Generated a table containing data from all reactions catalysed by
%   promiscuous enzymes. The table is then used with the
%   'getSubEnzymePools' function to obtain the subenzymes. This function is
%   not intended to be called outside of 'getSubEnzymePools' or for any 
%   other usage.
%
%   Usage
%       combinedTable = getTable(ecModel, rxnsInd, substrateMets, productMets, substrateCoeffs, productCoeffs)
%
%   Parameters
%       ecModel                       a GECKO3 ecModel structure
%       rxnsInd                       cell containing the indexes of all
%                                     reactions
%       substrateMets                 cell containing the substrate mets of
%                                     all rxns
%       productMets                   cell containing the product mets of
%                                     all rxns
%       substrateCoeffs               cell containing the stoichiometric
%                                     coefficient of all substrates
%       productCoeffs                 cell containing the stoichiometric
%                                     coefficient of all products
%
%   Outputs
%       combinedTable                 a table combining all inputs
%
% .. Author:
%       - Mauricio Ferreira &         2023.09.29
%         Eduardo Almeida


% Retrieve rxn IDs associated reactions containing kcat values
ecModel = buildRxnGeneMat(ecModel);

ecRxns = cell(numel(rxnsInd), 1);

for i = 1:numel(rxnsInd)
    ecRxns{i,1} = ecModel.rxns{rxnsInd{i}};
end

% Retrieve genes associated reactions containing kcat values
ecGenes = cell(length(rxnsInd), 1);

for i = 1:numel(rxnsInd)
    ecGenes{i,1} = ecModel.genes(find(ecModel.rxnGeneMat(rxnsInd{i},:)));
end

%% Join cell arrays
% Join substrates and products
combinedMets = cell(numel(rxnsInd), 1);

for i = 1:numel(rxnsInd)
    currentSubstrateMets = substrateMets{i};
    currentProductMets = productMets{i};
    combinedMets{i} = [currentSubstrateMets; currentProductMets];
end

combinedCoeffs = cell(numel(rxnsInd), 1);

for i = 1:numel(rxnsInd)
    combinedCoeffs{i} = [substrateCoeffs{i}, productCoeffs{i}];
end

% Combine all data into a table
combinedTable = table();
combinedTable.rxns = ecRxns;
combinedTable.mets = combinedMets;
combinedTable.stoichCoeffs = combinedCoeffs;
combinedTable.genes = ecGenes;

% Format genes to cell array
for i = 1:numel(combinedTable.genes)
    currentCell = combinedTable.genes{i};
    % GPRs of GECKO3-generated pcGEMs are either single enzymes or enzyme complexes ("AND")
    combinedTable.genes{i} = strjoin(currentCell, ' and '); 
end

% Remove rows with empty genes
combinedTable = combinedTable(cellfun(@(x) numel(x) > 1, combinedTable.genes), :);

end