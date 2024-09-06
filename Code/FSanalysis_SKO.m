%% Flux-sum analysis
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
% 

%%
clear;clc

%% Loading model and results
% load('SKO_redist_E1_glucose.mat')

% FSmatrix = readmatrix("FSmatrix_SKO_glucose.csv");
% FSmatrixWT = readmatrix("FSmatrix_SKO_glucose_WT.csv");

% load('../Models/eciML1515u_TurNuP_CORAL_Ratio.mat');
% 
% FSmatrix = readmatrix("FSmatrix_SKO_glucose_TurNuP.csv");
% FSmatrixWT = readmatrix("FSmatrix_SKO_glucose_TurNuP_WT.csv");

%% Calculate delta
FSmatrix_Delta = FSmatrix - FSmatrixWT;

% intervals = {[1834, 2847], [2892, 10153], [10296, 11821]};
intervals = {[1834, 2847], [2854, 4378], [4417, 11676]};

idxMet = true(size(model.mets, 1), 1);

for i = 1:length(intervals)
    idxMet(intervals{i}(1):intervals{i}(2)) = false;
end

modelMets = model.mets(idxMet);

max_values = max(FSmatrix_Delta);
columns_to_remove = max_values == 0;
FSmatrix_Delta = FSmatrix_Delta(:, ~columns_to_remove);
modelMetsDelta = modelMets(~columns_to_remove,:);

%% Export the table
filename = "FSmatrix_SKO_glucose_DLKcat_Delta.csv";
writematrix(FSmatrix_Delta, filename, 'Delimiter','\t')

%% Calculate ratio
max_values = max(FSmatrix);
columns_to_remove = max_values == 0;

FSmatrix_PreRatio = FSmatrix(:, ~columns_to_remove);
FSmatrixWT_Ratio = FSmatrixWT(:, ~columns_to_remove);

modelMetsRatio = modelMets(~columns_to_remove,:);

FSmatrix_PreRatio(FSmatrix_PreRatio==0) = eps;
FSmatrixWT_Ratio(FSmatrixWT_Ratio==0) = eps;

FSmatrix_Ratio = FSmatrixWT_Ratio ./ FSmatrix_PreRatio;

