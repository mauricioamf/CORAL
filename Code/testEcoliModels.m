function [solution1, solution2, solution3] = testEcoliModels(ecModel, verbose)
%
%   Tests if E. coli ecModels are working properly and can predict growth
%   under standard simulation conditions. The first conditions is a FBA
%   without adding any constraints. The second condition is a FBA with
%   added constraints on exchange reactions to simulate the M9 medium
%   composition and nutrient availability. The third condition is a FBA
%   with a fixed growth rate to simulate chemostat conditions. This
%   function was built for and only tested with the iML1515 model.
%
%   Usage
%       [solution1, solution2, solution3] = testEcoliModels(model, verbose)
%
%   Parameters
%       ecModel             (struct) a GECKO3 ecModel structure
%       verbose             (double) how much text should be reported (opt, default = 0)
%                               0 - Minimal report, only mentions if models
%                               work or not
%                               1 - Expanded report
%                               2 - Expanded report and print exchange fluxes
%                               3 - Expanded report and print all fluxes
%
%   Outputs
%       solution1           (struct) a struct containing the solutions for 
%                           the first condition
%       solution2           (struct) a struct containing the solutions for 
%                           the second condition
%       solution3           (struct) a struct containing the solutions for 
%                           the third condition
%
% .. Author:
%       - Mauricio Ferreira       2023.10.04

if nargin<2
    verbose=0;
end

if ~isfield(ecModel,'ec')
    error('.ec field does not exist, ecModel must be a GECKO3-generated ecModel')
end

%% Condition 1: Checking fluxes without adjusting parameters 
if verbose >= 1
    disp('Solving FBA without adjusting parameters')
end

solution1 = optimizeCbModel(ecModel, 'max');

% enzymeIds = find(~cellfun('isempty',strfind(model.rxnNames,'subpool_')));

if verbose >= 1
    try
        if solution1.stat == 1
            disp('Optimal solution found!')
            if verbose == 2
                disp('First solution ');
                printFluxes(ecModel, solution1.x, true);
            elseif verbose == 3
                disp('First solution ');
                printFluxes(ecModel, solution1.x, false);
            end
        end
    catch
        warning("No feasible solution for FBA without adjusting parameters")
    end

    % fprintf('\n');
    % [numE,~] = size(nonzeros(solution1.x(enzymeIds)));
    % formatSpecNUM = "Number of predicted enzymes: %f";
    % fprintf(formatSpecNUM, numE);
    % fprintf('\n');
end

%% Condition 2: Setting model constraints and objetive funtion
if verbose >= 1
    fprintf('\n')
    disp('Solving FBA using M9 medium composition and maximizing the growth rate')
end

model_temp = ecModel;

% M9 media without carbon
exchangeIndex = find(contains(model_temp.rxnNames, "exchange"));
exchangeIDs = model_temp.rxns(exchangeIndex);
exchangeIDs(end) = [];

model_temp = changeRxnBounds(model_temp, exchangeIDs, 0, 'l');

M9_components = ["EX_pi_e", "EX_co2_e", "EX_fe3_e", "EX_h_e", ...
    "EX_mn2_e", "EX_fe2_e", "EX_zn2_e", "EX_mg2_e", ...
    "EX_ca2_e", "EX_ni2_e", "EX_cu2_e", "EX_sel_e", ...
    "EX_cobalt2_e", "EX_h2o_e", "EX_mobd_e", "EX_so4_e", ...
    "EX_nh4_e", "EX_k_e", "EX_na1_e", "EX_cl_e", ...
    "EX_o2_e", "EX_tungs_e", "EX_slnt_e"];

model_temp = changeRxnBounds(model_temp, M9_components, -1000, 'l');

% Carbon sources
model_temp = changeRxnBounds(model_temp, 'EX_glc__D_e', -10, 'l'); %D-glucose
% model = changeRxnBounds(model, 'EX_glyc_e', 0, 'l'); %glycerol
% model = changeRxnBounds(model, 'EX_xyl__D_e', 0, 'l'); %D-xylose
% model = changeRxnBounds(model, 'EX_fuc__L_e', 0, 'l'); %L-fucose
% model = changeRxnBounds(model, 'EX_arab__L_e', 0, 'l'); %L-arabinose
% model = changeRxnBounds(model, 'EX_ac_e', 0, 'l'); %L-arabinose

% Defining the objetive
model_temp.c(:) = 0;
model_temp = changeObjective(model_temp, 'BIOMASS_Ec_iML1515_core_75p37M', 1);

solution2 = optimizeCbModel(model_temp);

if verbose >= 1
    try
        if solution2.stat == 1
            disp('Optimal solution found!')
            if verbose == 2
                disp('Second solution ');
                printFluxes(ecModel, solution2.x, true);
            elseif verbose == 3
                disp('Second solution ');
                printFluxes(ecModel, solution2.x, false);
            end
        end
    catch
        warning("No feasible solution for FBA with M9 medium composition")
    end

    % fprintf('\n');
    % [numE,~] = size(nonzeros(solution2.x(enzymeIds)));
    % formatSpecNUM = "Number of predicted enzymes: %f";
    % fprintf(formatSpecNUM, numE);
    % fprintf('\n');
end

%% FBA under a fixed growth rate
if verbose >= 1
    fprintf('\n')
    disp('Solving FBA with a fixed growth rate')
end

model_temp = ecModel;

model_temp.c(:) = 0;
model_temp = changeObjective(model_temp, 'BIOMASS_Ec_iML1515_core_75p37M', 1);
model_temp = changeRxnBounds(model_temp, 'BIOMASS_Ec_iML1515_core_75p37M', 0.1, 'b'); 

solution3 = optimizeCbModel(model_temp);

if verbose >= 1
    try
        if solution3.stat == 1
            disp('Optimal solution found!')
            if verbose == 2
                disp('Second solution ');
                printFluxes(ecModel, solution3.x, true);
            elseif verbose == 3
                disp('Second solution ');
                printFluxes(ecModel, solution3.x, false);
            end
        end
    catch
        warning("No feasible solution for FBA with fixed growth rate")
    end

    % fprintf('\n');
    % [numE,~] = size(nonzeros(solution3.x(enzymeIds)));
    % formatSpecNUM = "Number of predicted enzymes: %f";
    % fprintf(formatSpecNUM, numE);
    % fprintf('\n');
end


%% Report the test results
if solution1.stat == 1 && solution2.stat == 1 && solution3.stat == 1
        fprintf('\n')
        disp("Model successfully predicts growth in all tested conditions")
    else
        error("Model fails to predict growth in one or all tested conditions")
end