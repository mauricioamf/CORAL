function [SolutionEKcat, SolutionTheta] = solveCoralLP(ecModel, biomassRxnID, vBio, lambda, verbose)
%
%   Optimises Z = v_bio - lambda * sum(theta_i). First, an modified
%   implementation of pFBA is solved to obtain the enzyme usage
%   distribution of native reactions. This modified implementation swaps
%   the objective function 'min Z = sum(v_g)' for all fluxes v with a gene
%   association g, to the objetive min Z = sum(Ei * Kcat_ij), for all
%   enzyme usage E of enzyme i catalysing a reaction j.
%
%   Using the predicted enzyme usage distribution, the second optimization
%   then multiplies it with the kcats of each enzyme enzyme to calculate 
%   theta_E1Kn.
%
%   Usage
%       [SolutionEKcat, SolutionTheta] = solveCoralLP(ecModel, biomassRxnID, vBio, lambda, verbose)
%
%   Parameters
%       ecModel             (struct) a GECKO3 ecModel structure
%       biomassRxnID        (char)   the rxn ID for biomass pseudoreaction
%       vBio                (double) specific growth rate
%       lambda              (double) weighting factor for theta (opt, default = 1)
%       verbose             (double) how much text should be reported (opt, default = 1)
%                               0 - No report
%                               1 - Report on results from CORAL LP
%                               2 - Report on all optimizations
%
%   Outputs
%       SolutionEKcat       (struct) a struct containing the solutions for
%                           the first condition
%       SolutionTheta       (struct) a struct containing the solutions for
%                           the second condition
%
% .. Author:
%       - Mauricio Ferreira       2023.10.27

%% Defaut parameters and other checks
if nargin<5 || isempty(verbose)
    verbose=1;
end
if nargin<4 || isempty(lambda)
    lambda=1;
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
undIndex = find(contains(modelpFBA.rxns, 'u0'));
undIDs = modelpFBA.rxns(undIndex);

modelpFBA = changeRxnBounds(modelpFBA, undIDs, 0, 'b');

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
    if verbose == 2
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

%% Calculate theta
E1Table = table();
E1IDs = enzymeTable.rxnUpIdx(find(endsWith(enzymeTable.enzNames, '_1')));
E1Table.enzUniprot = enzymeTable.enzUniprot(find(endsWith(enzymeTable.enzNames, '_1')));
E1Table.E1 = SolutionEKcat.x(E1IDs);

enzymeTable = outerjoin(enzymeTable, E1Table);

enzymeTable.Properties.VariableNames{5} = 'enzUniprot';
enzymeTable.enzUniprot_E1Table = [];

enzymeTable.fluxes = SolutionEKcat.x(enzymeTable.enzRxns);

enzymeTable.thetaKnE1 = enzymeTable.kcats .* enzymeTable.E1;

enzymeTable.En = SolutionEKcat.x(enzymeTable.rxnUpIdx);
enzymeTable.thetaKnEn = enzymeTable.kcats .* enzymeTable.En;

% Add "pseudocount" to enzymes with zero usage
% solutionArray = Solution.x(enzymeTable.rxnUpIdx);
% lowestEn = max(solutionArray(solutionArray<0));
% lowestEn = lowestEn / 1e6;
% 
% enzymeTable.En = enzymeTable.En + lowestEn;


%% Optimises Z = v_bio - lambda * sum(theta_i)
modelTheta = ecModel;

% Fix biomass pseudoreaction
modelTheta = changeRxnBounds(modelTheta, biomassRxnID, vBio, 'b');

% Construct optimization problem
modelTheta = buildLPproblemFromModel(modelTheta);

% [~,nRxns] = size(modelTheta.A);

% Define upper and/or lower bounds
modelTheta.lb(modelTheta.lb==-Inf) = -1000;
% modelTheta.lb(enzymeTable.enzIdx) = lambda.*enzymeTable.thetaKnE1;
% modelTheta.lb(enzymeTable.rxnUpIdx) = lambda.*enzymeTable.thetaKnE1;

for t = 1:length(enzymeTable.thetaKnE1)
    if enzymeTable.thetaKnE1(t) ~= 0
        modelTheta.lb(enzymeTable.enzIdx(t)) = lambda.*enzymeTable.thetaKnE1(t);
        modelTheta.ub(enzymeTable.enzIdx(t)) = lambda.*-enzymeTable.thetaKnE1(t);
        % modelTheta.lb(enzymeTable.enzIdx(t)) = lambda.*enzymeTable.thetaKnEn(t);
        % modelTheta.ub(enzymeTable.enzIdx(t)) = lambda.*-enzymeTable.thetaKnEn(t);
    end
end

modelTheta.ub(modelTheta.ub==Inf) = 1000;

% Define objective vector
modelTheta.c(:) = 0;
modelTheta.c(biomassRxnID) = vBio;
modelTheta.c(enzymeTable.enzIdx) = lambda.*enzymeTable.thetaKnE1;
% modelTheta.c(enzymeTable.enzIdx) = lambda.*enzymeTable.thetaKnEn;

% Other parameters
modelTheta.osense = -1;

% Solve the problem
SolutionTheta = solveCobraLP(modelTheta);
SolutionTheta.x = SolutionTheta.full;

% Get the solution(s)
if SolutionTheta.stat == 1
    if verbose >= 1
        disp('CORAL LP PREDICTIONS')

        printFluxes(ecModel, SolutionTheta.x, true);

        % How many enzymes were predicted?
        fprintf('\n');
        [numE,~] = size(nonzeros(SolutionTheta.x(enzUsageIds)));
        [numSub,~] = size(nonzeros(SolutionTheta.x(subpoolIds)));
        formatSpecNUM = "Number of predicted enzymes: %f subenzymes and %f subpools";
        fprintf(formatSpecNUM, numE, numSub);
        fprintf('\n');
        
        % How many underground reactions were used?
        UndIds = find(~cellfun('isempty',strfind(ecModel.rxnNames,'u0'))); 
        [numUnd,~] = size(nonzeros(SolutionTheta.x(UndIds)));
        formatSpecNUM = "Number of predicted underground reactions: %f";
        fprintf(formatSpecNUM, numUnd);
        fprintf('\n');
        fprintf('\n');
    end
else
    error("No feasible solution for modelTheta")
end

end