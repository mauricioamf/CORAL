function [substrateMets, productMets, substrateCoeffs, productCoeffs] = getRxnData(ecModel, rxnsInd)
%
%   Retrieve substrate and product information from the stoichiometric 
%   matrix. This function is not intended to be called outside of 
%   'expandEnzComplex', 'getSubEnzymePools', or for any other usage.
%
%   Usage
%       [substrateMets, productMets, substrateCoeffs, productCoeffs] = getRxnData(ecModel, rxnsInd)
%
%   Parameters
%       ecModel                       (struct) a GECKO3 ecModel structure
%       rxnsInd                       (cell) cell containing the indexes of
%                                     all reactions
%
%   Outputs
%       substrateMets                 (cell) cell containing the substrate
%                                     mets of all rxns
%       productMets                   (cell) cell containing the product 
%                                     mets of all rxns
%       substrateCoeffs               (cell) cell containing the stoichiometric
%                                     coefficient of all substrates
%       productCoeffs                 (cell) cell containing the stoichiometric
%                                     coefficient of all products
%
% .. Author:
%       - Mauricio Ferreira         2023.09.29

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

% Loop over stoichiometric matrix to get metabolite IDs and coefficients
substrateIDs = cell(numel(rxnsInd), 1);
substrateCoeffs = cell(numel(rxnsInd), 1);
productIDs = cell(numel(rxnsInd), 1);
productCoeffs = cell(numel(rxnsInd), 1);

for i=1:numel(rxnsInd)
    substrateIDs{i,1}    = find(ecModel.S(:,rxnsInd{i}) < 0);
    substrateCoeffs{i,1} = ecModel.S(substrateIDs{i,1}(:), rxnsInd{i,1})';
    productIDs{i,1}      = find(ecModel.S(:,rxnsInd{i}) > 0);
    productCoeffs{i,1}   = ecModel.S(productIDs{i,1}(:), rxnsInd{i,1})';
end

% Refactor sparse double to double
substrateCoeffsFull = cell(numel(substrateCoeffs), 1);

for i = 1:numel(substrateCoeffs)
    substrateCoeffsFull{i} = full(substrateCoeffs{i});
end

substrateCoeffs = substrateCoeffsFull;

% Refactor sparse double to double
ProductCoeffsFull = cell(numel(productCoeffs), 1);

for i = 1:numel(productCoeffs)
    ProductCoeffsFull{i} = full(productCoeffs{i});
end

productCoeffs = ProductCoeffsFull;

% Get mets inside substrateIDs and productIDs
substrateMets = cell(numel(rxnsInd), 1);
productMets = cell(numel(rxnsInd), 1);

for i = 1:numel(rxnsInd)
    currentSubstrateIDs = substrateIDs{i};
    currentSubstrateNames = cell(numel(currentSubstrateIDs), 1);
    
    currentProductIDs = productIDs{i};
    currentProductNames = cell(numel(currentProductIDs), 1);


    for j = 1:numel(currentSubstrateIDs)
        currentSubstrateID = currentSubstrateIDs(j);
        currentSubstrateNames{j} = ecModel.mets{currentSubstrateID};
    end

    for k = 1:numel(currentProductIDs)
        currentProductID = currentProductIDs(k);
        currentProductNames{k} = ecModel.mets{currentProductID};
    end
    
    substrateMets{i} = currentSubstrateNames;
    productMets{i} = currentProductNames;
end
end