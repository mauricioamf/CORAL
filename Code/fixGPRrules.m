function ecModel = fixGPRrules(ecModel)

% Remove .rules field. It can cause problems if it exists
% alongside model.grRules.
if isfield(ecModel,'rules')
    ecModel = rmfield(ecModel, "rules");
end

% Fix symbols and trim spaces
fixedGrRules = cell(size(ecModel.grRules));
for i = 1:numel(ecModel.grRules)
    fixedGrRules{i} = strrep(ecModel.grRules{i}, '( ', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '(( ', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '(', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '((', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ')', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, '))', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ' )', '');
    fixedGrRules{i} = strrep(fixedGrRules{i}, ' ))', '');
end
ecModel.grRules = strtrim(fixedGrRules);

fixedGenes = cell(size(ecModel.genes));
for i = 1:numel(ecModel.genes)
    fixedGenes{i} = strrep(ecModel.genes{i}, '( ', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '(( ', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '(', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '((', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ')', '');
    fixedGenes{i} = strrep(fixedGenes{i}, '))', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ' )', '');
    fixedGenes{i} = strrep(fixedGenes{i}, ' ))', '');
end
ecModel.genes = strtrim(fixedGenes);

end