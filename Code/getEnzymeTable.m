function enzymeTable = getEnzymeTable(ecModel)
%
%   Generates a table containing data from all enzymes and reactions
%   asssociated with enzymes. This function is not intended to be called 
%   outside of 'solveKcatE', 'solveCoralLP' or for any other usage.
%
%   Usage
%       enzymeTable = getEnzymeTable(ecModel)
%
%   Parameters
%       ecModel                       (struct) a GECKO3 ecModel structure
%
%   Outputs
%       enzymeTable                   (table) a table contaning enzyme
%                                     usage information
%
% .. Author:
%       - Mauricio Ferreira           2023.09.29

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

enzIdx = find(startsWith(ecModel.metNames,'prot_'));
enzIdx(end-1:end) = [];
enzCoef = cell(length(enzIdx),1);
rxnUpIdx = [];
enzRxns = [];

for i=1:numel(enzIdx)
    % pIdx(i,1) = find(strcmpi(ecModel.metNames,join(['prot_' char(proteins(i)) '_'],"")));
    rxnIdx = find(ecModel.S(enzIdx(i),:) < 0);
    rxnUpIdx = [rxnUpIdx rxnIdx(end)];
    rxnIdx(end) = [];
    enzRxns = [enzRxns rxnIdx];
    enzCoef{i,1} = abs(ecModel.S(enzIdx(i),rxnIdx));
end

rxnUpIdx = rxnUpIdx';
enzRxns = enzRxns';
enzRxnNames = ecModel.rxnNames(enzRxns);
enzNames = ecModel.metNames(enzIdx);

enzUniprot = erase(enzNames, 'prot_');
pattern = '(_\d+)$';
enzUniprot_temp = cellfun(@(str) regexprep(str, pattern, ''), enzUniprot, 'UniformOutput', false);
enzUniprot = enzUniprot_temp;

enzymeTable = table();
enzymeTable.enzRxns = enzRxns;
enzymeTable.enzRxnNames = enzRxnNames;
enzymeTable.enzIdx = enzIdx;
enzymeTable.rxnUpIdx = rxnUpIdx;
enzymeTable.enzUniprot = enzUniprot;
enzymeTable.enzNames = enzNames;
enzymeTable.enzCoef = enzCoef;

preCoralEcFields = table();
preCoralEcFields.enzUniprot = ecModel.ec.enzymes;
preCoralEcFields.MW = ecModel.ec.mw;

enzymeTable = outerjoin(enzymeTable, preCoralEcFields);
enzymeTable.Properties.VariableNames{5} = 'enzUniprot';
enzymeTable.enzUniprot_preCoralEcFields = [];
enzymeTable(end,:) = [];

enzymeTable(any(ismissing(enzymeTable),2), :) = [];

numericArray = zeros(numel(enzymeTable.enzCoef),1);

for k = 1:numel(numericArray)
    numericArray(k,1) = enzymeTable.enzCoef{k};
end

enzymeTable.enzCoef = numericArray;

MWs = enzymeTable.MW;

kcats = MWs ./ numericArray;

enzymeTable.kcats = kcats;

% clear i pattern rxnIdx enzCoef rxnUpIdx enzIdx enzNames enzRxnNames enzRxns enzUniprot enzUniprot_temp kcats MWs numericArray preCoralEcFields

end