function FVAtable = Fluxes_FVA(model)
% enzymeUsage_FVA
%   Perf
%
%   Usage: FVAtable = enzymeUsage_FVA(model,enzymes)
%
%   Ivan Domenzain.     Last edited 2020-05-27

% if nargin<2
%     enzymes = model.enzymes;
% end
%Get parsimonious protein usages
tempModel  = model;
% pool_indx  = find(strcmpi(model.rxns,'prot_pool_exchange'));
% prot_indxs = find(contains(model.rxnNames,'prot_'));
% prot_indxs = prot_indxs(1:(end-1));
% tempModel = setParam(tempModel, 'obj',pool_indx,-1);
sol       = solveLP(tempModel,1);
%initialize variables
ranges    = [];
minUsgs   = [];
maxUsgs   = [];
pFlux = [];

% Perform Flux Variability Analysis on enzyme usage
geneList = model.genes;

% Initialize a cell array to store the results
reactionIDs = cell(length(geneList), 1); % 1 column for associated reactions

for i = 1:length(geneList)
    gene = geneList{i};
    
    % Find the index of the gene in the gene list of the model
    geneIndex = find(strcmp(model.genes, gene));
    
    if ~isempty(geneIndex)
        % Get reaction IDs associated with the gene
        reactionIDs{i} = model.rxns(model.rxnGeneMat(:, geneIndex) > 0);
    else
        reactionIDs{i} = 'Gene not found in the model';
    end
end

% Collect all reaction IDs into a single cell array
allReactionIDs = cat(1, reactionIDs{:});

% Remove duplicate and empty reaction IDs
allReactionIDs = unique(allReactionIDs);
allReactionIDs = allReactionIDs(~cellfun('isempty', allReactionIDs));
allReactionIDs = allReactionIDs(~startsWith(allReactionIDs, 'usage_prot'));

[~, allReactionIDXs] = ismember(allReactionIDs, model.rxns);

if ~isempty(sol.x)
    pFluxes = sol.x(allReactionIDXs);
    %Loop through all the provided enzymes
    for i=1:length(allReactionIDs)
        if ~isempty(allReactionIDs{i})
            rxnIndx = find(contains(model.rxns,allReactionIDs{i}));
            AllIndx = strcmpi(model.rxns,allReactionIDs{i});
            model = setParam(model, 'obj', rxnIndx, -1);
            sol   = solveLP(model);
            if ~isempty(sol.f)
                minFlux = sol.x(rxnIndx); 
                model   = setParam(model, 'obj', rxnIndx, +1);
                sol     = solveLP(model);
                if ~isempty(sol.f)
                   %disp(['Ready with enzyme #' num2str(i)])
                   maxFlux = sol.x(rxnIndx); 
                else
                   maxFlux = nan; 
                end
            else
                minFlux = nan;
                maxFlux = nan; 
            end
            ranges    = [ranges; (maxFlux-minFlux)];
            minUsgs   = [minUsgs; minFlux]; 
            maxUsgs   = [maxUsgs; maxFlux];
            % pFlux = [pFlux;pFluxes(AllIndx)];
        else
            ranges    = [ranges; 0];
            minUsgs   = [minUsgs; 0]; 
            maxUsgs   = [maxUsgs; 0];
            % pFlux = [pFlux;0];
        end
    end
else
    disp('Model is not feasible')
end
varNamesT = {'enzymes' 'ranges' 'minFlux' 'maxFlux' 'pFlux'};
FVAtable  = table(allReactionIDs,ranges,minUsgs,maxUsgs,pFlux,'VariableNames', varNamesT);
end