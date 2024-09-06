%%
% initCobraToolbox(false);
cd ~/CORAL/Code/
clc;clear

%%
load('~/CORALModels/eciML1515u_DLKcat_CORAL.mat');

%%
customKcats = readtable("../Data/customKcats_sidesubpools.csv");

enzymeTable = getEnzymeTable(model);
enzymeTable.Properties.VariableNames(6) = "subpool";
enzymeTable = removevars(enzymeTable, "kcats");
enzymeTable = removevars(enzymeTable, "enzCoef");

customKcats = outerjoin(customKcats, enzymeTable);
rowsWithEmptyCells = any(ismissing(customKcats), 2);
customKcats(rowsWithEmptyCells, :) = [];
customKcats.enzCoef = customKcats.MW ./ customKcats.kcat;

%%
rxnIds = customKcats.rxns;
kcat = customKcats.kcat;
subpool = customKcats.subpool_customKcats;
enzCoef = customKcats.enzCoef;
rxn = customKcats.enzRxns;
met = customKcats.enzIdx;

model = setKcatForReactions(model,rxnIds,kcat);
% model = applyKcatConstraints(model);

for i=1:length(rxn)
    model.S(met(i), rxn(i)) = -enzCoef(i);
end

%%
% enzIdx = find(contains(model.metNames,subpool));
% % enzIdx(end-1:end) = [];
% % enzCoef = cell(length(enzIdx),1);
% rxnUpIdx = [];
% enzRxns = [];
% 
% for i=1:numel(enzIdx)
%     % pIdx(i,1) = find(strcmpi(model.metNames,join(['prot_' char(proteins(i)) '_'],"")));
%     rxnIdx = find(model.S(enzIdx(i),:) < 0);
%     rxnUpIdx = [rxnUpIdx rxnIdx(end)];
%     rxnIdx(end) = [];
%     enzRxns = [enzRxns rxnIdx];
%     % enzCoef{i,1} = abs(model.S(enzIdx(i),rxnIdx));
%     model.S(enzIdx(i),rxnIdx) = enzCoef{i,1};
% end

%%
save('../Models/eciML1515u_CORAL_TurNuP.mat','model');