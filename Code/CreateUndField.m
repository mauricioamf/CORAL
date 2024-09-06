function ecModel = CreateUndField(ecModel, modelAdapter)

uniprotDB = loadDatabases('both', modelAdapter);
uniprotDB = uniprotDB.uniprot;

% Create .und field and populate und.rxns, und.genes and und.proteins
rxnWithGene  = find(sum(ecModel.rxnGeneMat,2));

und.rxns      = ecModel.rxns(rxnWithGene);
emptyCell    = cell(numel(rxnWithGene),1);
emptyCell(:) = {''};
emptyVect    = zeros(numel(rxnWithGene),1);

und.kcat      = emptyVect;
und.source    = emptyCell; % Strings, like 'dlkcat', 'manual', 'brenda', etc.
und.notes     = emptyCell; % Additional comments
und.eccodes   = emptyCell;
und.concs     = emptyVect;

uniprotCompatibleGenes = modelAdapter.getUniprotCompatibleGenes(ecModel.genes);
[Lia,Locb] = ismember(uniprotCompatibleGenes,uniprotDB.genes);

uniprot = modelAdapter.getUniprotIDsFromTable(uniprotCompatibleGenes);
if ~isequal(uniprot,uniprotCompatibleGenes)
    uniprot(cellfun(@isempty,uniprot)) = {''};
    [Lia,Locb] = ismember(uniprot,uniprotDB.ID);
end
noUniprot  = uniprotCompatibleGenes(~Lia);
und.genes        = ecModel.genes(Lia); %Will often be duplicate of model.genes, but is done here to prevent issues when it is not.
und.enzymes      = uniprotDB.ID(Locb(Lia));
und.mw           = uniprotDB.MW(Locb(Lia));
und.sequence     = uniprotDB.seq(Locb(Lia));
%Additional info
und.concs        = nan(numel(und.genes),1); % To be filled with proteomics data when available

und.rxnEnzMat = zeros(numel(rxnWithGene),numel(und.genes)); % Non-zeros will indicate the number of subunits
for r=1:numel(rxnWithGene)
    rxnGenes   = ecModel.genes(find(ecModel.rxnGeneMat(rxnWithGene(r),:)));
    [~,locEnz] = ismember(rxnGenes,und.genes); % Could also parse directly from rxnGeneMat, but some genes might be missing from Uniprot DB
    if locEnz ~= 0
        und.rxnEnzMat(r,locEnz) = 1; %Assume 1 copy per subunit or enzyme, can be modified later
    end
end

% Populate und.eccodes

rxnNames = und.rxns;
[~,rxnIdxs] = ismember(rxnNames,ecModel.rxns);

% Check if eccodes are valid
eccodes = ecModel.eccodes;
invalidEC = regexprep(eccodes,'(\d\.(\w|-)+\.(\w|-)+\.(\w|-)+)(;\w+\.(\w|-)+\.(\w|-)+\.(\w|-)+)*(.*)','$3');
invalidEC = ~cellfun(@isempty,invalidEC);
invalidECpos = find(invalidEC);
if any(invalidECpos)
    invalidEC = ecModel.eccodes(invalidEC);
    if nargout<2
        fprintf('Skipped incorrectly formatted EC numbers, rerun getECfromGEM with all outputs to get a list.\n')
    else
        fprintf('Skipped incorrectly formatted EC numbers.\n')
    end
    eccodes(invalidECpos)={''};
else
    invalidEC = [];
end

und.eccodes = eccodes(rxnIdxs);

% Populate und.kcats

enzymeTable = getEnzymeTable(ecModel);
enzymeTable = sortrows(enzymeTable,"enzRxns","ascend");

und.kcat = enzymeTable.kcats;

ecModel.und = und;
end