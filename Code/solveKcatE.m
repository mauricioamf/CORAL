function [SolutionEKcat] = solveKcatE(ecModel, biomassRxnID, vBio, blockUnd, verbose)
%
%   Optimises Z = v_bio - lambda * sum(theta_i)
%
%   Usage
%       [SolutionEKcat, SolutionTheta] = solveCoralLP(ecModel, biomassRxnID, vBio, lambda, verbose)
%
%   Parameters
%       ecModel             (struct) a GECKO3 ecModel structure
%       biomassRxnID        (char)   the rxn ID for biomass pseudoreaction
%       vBio                (double) specific growth rate
%       blockUnd            (logical) if underground reactions should be 
%                           blocked (opt, default = false)
%       verbose             (double) how much text should be reported (opt, default = 1)
%                               0 - No report
%                               1 - Report on results
%
%   Outputs
%       SolutionEKcat       (struct) a struct containing the solutions for
%                           the first condition
%
% .. Author:
%       - Mauricio Ferreira       2023.10.27

%% Defaut parameters and other checks
if nargin<5 || isempty(verbose)
    verbose=1;
end
if nargin<4 || isempty(blockUnd)
    blockUnd=false;
end
if vBio == 0
    warning("The LP might have unexpected behaviours with zero growth");
end

indexBiom=find(strcmp(biomassRxnID,ecModel.rxns),1);
if isempty(indexBiom)
    error(['Biomass pseudoreaction ' biomassRxnID ' does not exist in the model. Check for spelling mistakes and try again']);
end

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

%% Retrieve kcat values
enzymeTable = getEnzymeTable(ecModel);

%% Calculate enzyme usage distribution using a modified pFBA
modelpFBA = ecModel;
modelpFBA = buildRxnGeneMat(modelpFBA);

% gprAssociated = find(sum(modelpFBA.rxnGeneMat,2)>0);

% Fix biomass
modelpFBA = changeRxnBounds(modelpFBA, biomassRxnID, vBio, 'b');

enzUsageIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'usage_prot_')));
enzUsageIds(end) = [];

subpoolIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'subpool_')));

% Block underground reactions
if blockUnd == true
    undIndex = find(contains(modelpFBA.rxns, 'u0'));
    undIDs = modelpFBA.rxns(undIndex);
    modelpFBA = changeRxnBounds(modelpFBA, undIDs, 0, 'b');
end
% reacIndex = find(contains(modelpFBA.rxns, 'reac'));
% reacIDs = modelpFBA.rxns(reacIndex);
% 
% modelpFBA = changeRxnBounds(modelpFBA, reacIDs, 0, 'b');

pFBAmodel = buildLPproblemFromModel(modelpFBA);

[~,nRxns] = size(pFBAmodel.A);

pFBAmodel.lb(pFBAmodel.lb==-Inf) = -1000;
pFBAmodel.ub(pFBAmodel.ub==Inf) = 1000;

pFBAmodel.c = zeros(nRxns,1);
pFBAmodel.c(enzUsageIds) = 1;
pFBAmodel.c(enzUsageIds) = enzymeTable.kcats.*pFBAmodel.c(enzUsageIds);

pFBAmodel.osense = 1;

SolutionEKcat = solveCobraLP(pFBAmodel);

SolutionEKcat.x = SolutionEKcat.full;

% Get the solution(s)

if SolutionEKcat.stat == 1
    if verbose == 1
        disp('pFBA PREDICTIONS')

        printFluxes(ecModel, SolutionEKcat.x, true);

        % How many enzymes were predicted?
        fprintf('\n');
        [numE,~] = size(nonzeros(SolutionEKcat.x(enzUsageIds)));
        [numSub,~] = size(nonzeros(SolutionEKcat.x(subpoolIds)));
        formatSpecNUM = "Number of predicted enzymes: %f subenzymes and %f subpools";
        fprintf(formatSpecNUM, numE, numSub);
        fprintf('\n');
        
        % How many underground reactions were used?
        UndIds = find(~cellfun('isempty',strfind(modelpFBA.rxnNames,'u0'))); 
        [numUnd,~] = size(nonzeros(SolutionEKcat.x(UndIds)));
        formatSpecNUM = "Number of predicted underground reactions: %f";
        fprintf(formatSpecNUM, numUnd);
        fprintf('\n');
        fprintf('\n');
    end
else
    error("No feasible solution for pFBA")
end

end