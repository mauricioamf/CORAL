function ecModelNoComplex = expandEnzComplex(ecModel)
%
%   Expands a model which uses 'AND' gene associations for one reaction.
%   Each such reaction is split into several reactions, each under the 
%   control of only one gene.
%
%   Here we take a different approach than for 'OR' gene associations by
%   breaking a reaction catalysed by an enzyme complex into separate 
%   partial reactions catalysed by one enzyme. The partial reactions 
%   include pseudo-metabolites to represent an intermediate state similar 
%   to arm reactions in the GECKO Toolbox 2.
%
%   Usage
%       ecModelNoComplex = expandEnzComplex(ecModel)
%
%   Parameters
%       ecModel                       (struct) a GECKO3 ecModel structure
%
%   Outputs
%       ecModelNoComplex              (struct) a GECKO3 ecModel structure 
%                                     without complex GPR rules containing
%                                     "AND" rules
%
% .. Author:
%       - Mauricio Ferreira           2023.09.29
%

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

%% Retrieve IDs and coefficients of mets and rxns associated with kcats
rxnsInd = num2cell(1:1:length(ecModel.rxns))';
[substrateMets, productMets, substrateCoeffs, productCoeffs] = getRxnData(ecModel, rxnsInd);

%% Loop over reactions with 'AND' GPR rules and slip into partial reactions
hasAND = {};

if isfield(ecModel,'rules')
    ecModel = rmfield(ecModel, "rules");
end

% Iterate over all reactions in the model
for i = 1:length(ecModel.rxns)
    currentRxnMets = substrateMets{i};

    % Check if the GPR contains an 'AND'
    if contains(ecModel.grRules{i}, ' and ') && any(startsWith(currentRxnMets, 'prot_'))

        % Print reaction
        disp(['Working on reaction ' ecModel.rxns{i}]);

        % Add reaction to hasAND cell
        hasAND = [hasAND ecModel.rxns{i}];

        % Split the GPR into individual genes
        genes = strsplit(ecModel.grRules{i}, ' and ')';

        % Create a pseudo-metabolite for each gene
        pMets = cellstr(strcat('pmet_', ecModel.rxns{i}, '_000', num2str((1:length(genes)-1)')));

        % Some pseudometabolites have 5 digits, this can cause problems later
        for psm = 1:length(pMets)
            % Check if the string ends with a 5-digit number
            if ~isempty(regexp(pMets{psm}, '\d{5}$', 'once'))
                fiveDigitNumber = regexp(pMets{psm}, '\d{5}$', 'match', 'once');
                trimmedNumber = fiveDigitNumber(2:end);
                modifiedString = regexprep(pMets{psm}, '\d{5}$', trimmedNumber);
                pMets{psm} = modifiedString;
            else
                pMets{psm} = pMets{psm};
            end
        end

        metsToAdd = struct();
        metsToAdd.mets = pMets;
        metsToAdd.compartments = cell(length(metsToAdd.mets), 1);
        for k = 1:length(metsToAdd.compartments)
            metsToAdd.compartments{k} = "c";
        end
        ecModel = addMets(ecModel, metsToAdd);

        % Create new reactions for each gene
        for j = 1:length(genes)

            if j == 1 % The first reaction. Uses the original substrates
                rxnsToAdd = struct();
                rxnsToAdd.rxns = [ecModel.rxns{i} '_Partial_000' num2str(j)];

                % Add substrates
                currentRxn = substrateMets{i};
                for m = 1:length(currentRxn)
                    if ~startsWith(currentRxn{m}, 'prot_')
                        rxnsToAdd.mets{m,1} = currentRxn{m};
                    end
                end
                for c = 1:length(rxnsToAdd.mets)
                    rxnsToAdd.stoichCoeffs(c) = substrateCoeffs{i}(c);
                end
                % Add protein pseudometabolite
                rxnProteins = cell(length(currentRxn), 1);
                for p = 1:length(currentRxn)
                    if startsWith(currentRxn{p}, 'prot_')
                        rxnProteins{p,1} = currentRxn{p};
                    end
                end
                rxnProteins = rxnProteins(~cellfun('isempty', rxnProteins));
                rxnsToAdd.mets{length(rxnsToAdd.mets)+1,1} = rxnProteins{1,1};
                % Add MW/kcat coefficient
                currentCoeffs = substrateCoeffs{i};
                coeffsProteins = cell(length(currentCoeffs), 1);
                for p = 1:length(currentCoeffs)
                    if startsWith(currentRxn(p), 'prot_')
                        coeffsProteins{p,1} = currentCoeffs(p);
                    end
                end
                coeffsProteins = coeffsProteins(~cellfun('isempty', coeffsProteins));
                rxnsToAdd.stoichCoeffs(length(rxnsToAdd.stoichCoeffs)+1) = coeffsProteins{1,1};
                % Add pseudo-product
                rxnsToAdd.mets{length(rxnsToAdd.mets)+1,1} = pMets{1,1};
                rxnsToAdd.stoichCoeffs(length(rxnsToAdd.stoichCoeffs)+1) = 1;
                % Copy annotations from original reaction
                rxnsToAdd = copyRxnInfo(ecModel, rxnsToAdd, i);
                % Add reactions to model
                ecModel = addRxns(ecModel, rxnsToAdd);

            elseif j == length(genes) % Last reaction. Produces the original products
                rxnsToAdd = struct();
                rxnsToAdd.rxns = [ecModel.rxns{i} '_Partial_000' num2str(j)];

                % Add pseudo-substrate
                rxnsToAdd.mets{1,1} = pMets{end,1};
                rxnsToAdd.stoichCoeffs(1,1) = -1;
                % Add protein pseudometabolite
                currentRxn = substrateMets{i};
                rxnProteins = cell(length(currentRxn), 1);
                for p = 1:length(currentRxn)
                    if startsWith(currentRxn{p}, 'prot_')
                        rxnProteins{p,1} = currentRxn{p};
                    end
                end
                rxnProteins = rxnProteins(~cellfun('isempty', rxnProteins));
                rxnsToAdd.mets{length(rxnsToAdd.mets)+1,1} = rxnProteins{end,1};
                % Add MW/kcat coefficient
                currentCoeffs = substrateCoeffs{i};
                coeffsProteins = cell(length(currentCoeffs), 1);
                for p = 1:length(currentCoeffs)
                    if startsWith(currentRxn(p), 'prot_')
                        coeffsProteins{p,1} = currentCoeffs(p);
                    end
                end
                coeffsProteins = coeffsProteins(~cellfun('isempty', coeffsProteins));
                rxnsToAdd.stoichCoeffs(length(rxnsToAdd.stoichCoeffs)+1) = coeffsProteins{end,1};
                % Add products
                currentRxn = productMets{i};
                for m = 1:length(currentRxn)
                    if ~startsWith(currentRxn{m}, 'prot_')
                        rxnsToAdd.mets{end+1,1} = currentRxn{m};
                    end
                end
                for c = 1:length(productCoeffs{i})
                    rxnsToAdd.stoichCoeffs(end+1) = productCoeffs{i}(c);
                end
                % Copy annotations from original reaction
                rxnsToAdd = copyRxnInfo(ecModel, rxnsToAdd, i);
                % Add reactions to model
                ecModel = addRxns(ecModel, rxnsToAdd);
            
            else % Intermediate reactions. Convert one pseudo-metabolite to the next
                rxnsToAdd = struct();
                rxnsToAdd.rxns = [ecModel.rxns{i} '_Partial_000' num2str(j)];
                % Add pseudo-substrate
                rxnsToAdd.mets{1,1} = pMets{j-1,1};
                rxnsToAdd.stoichCoeffs(1,1) = -1;
                % Add protein pseudometabolite
                currentRxn = substrateMets{i};
                rxnProteins = cell(length(currentRxn), 1);
                for p = 1:length(currentRxn)
                    if startsWith(currentRxn{p}, 'prot_')
                        rxnProteins{p,1} = currentRxn{p};
                    end
                end
                rxnProteins = rxnProteins(~cellfun('isempty', rxnProteins));
                rxnsToAdd.mets{length(rxnsToAdd.mets)+1,1} = rxnProteins{j,1};
                % Add MW/kcat coefficient
                currentCoeffs = substrateCoeffs{i};
                coeffsProteins = cell(length(currentCoeffs), 1);
                for p = 1:length(currentCoeffs)
                    if startsWith(currentRxn(p), 'prot_')
                        coeffsProteins{p,1} = currentCoeffs(p);
                    end
                end
                coeffsProteins = coeffsProteins(~cellfun('isempty', coeffsProteins));
                rxnsToAdd.stoichCoeffs(length(rxnsToAdd.stoichCoeffs)+1) = coeffsProteins{j,1};
                % Add pseudo-product
                rxnsToAdd.mets{length(rxnsToAdd.mets)+1,1} = pMets{j,1};
                rxnsToAdd.stoichCoeffs(length(rxnsToAdd.stoichCoeffs)+1) = 1;
                % Copy annotations from original reaction
                rxnsToAdd = copyRxnInfo(ecModel, rxnsToAdd, i);
                % Add reactions to model
                ecModel = addRxns(ecModel, rxnsToAdd);
            
            end
            
            % Assign the gene to the new reaction
            ecModel = changeGeneAssociation(ecModel, [ecModel.rxns{i} '_Partial_000' num2str(j)], genes{j});
        
        end
    end
end

% Remove original reactions from model
hasAND = hasAND';
ecModelNoComplex = removeReactions(ecModel, hasAND);

% Some reactions end in 5 digits, this can cause problems later
partialRxns = find(contains(ecModel.rxns, "_Partial_"));
partialRxnsIDs = ecModel.rxns(partialRxns);

for rxn5 = 1:length(partialRxnsIDs)
    % Check if the string ends with a 5-digit number
    if ~isempty(regexp(partialRxnsIDs{rxn5}, '\d{5}$', 'once'))
        fiveDigitNumber = regexp(partialRxnsIDs{rxn5}, '\d{5}$', 'match', 'once');
        trimmedNumber = fiveDigitNumber(2:end);
        modifiedString = regexprep(partialRxnsIDs{rxn5}, '\d{5}$', trimmedNumber);
        partialRxnsIDs{rxn5} = modifiedString;
    else
        partialRxnsIDs{rxn5} = partialRxnsIDs{rxn5};
    end
end

ecModel.rxns(partialRxns) = partialRxnsIDs;

if isfield(ecModel,'rules')
    ecModel = rmfield(ecModel, "rules");
end

%% All finished
numAND = size(hasAND,1);
disp(['Finished! A total of ' num2str(numAND) ' reactions were modified.'])
end