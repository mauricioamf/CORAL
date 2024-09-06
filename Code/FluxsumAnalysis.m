%% Flux-sum analysis
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
initCobraToolbox(false)
cd("/home/mferreira/Doutorado/4.Underground_metabolism/underground_metabolism/Code/");

%%
clear;clc

%% Loading model and results
% load('../Models/eciML1515u_v2_stage2_DLKcat.mat');

% load('../Models/eciML1515u_CORAL_BRENDA.mat');
% load('../Models/eciML1515u_CORAL_DLKcat_Ratio.mat');
% load('SKO_redist_E1_glucose.mat');

load('SKO_redist_Delta_Ratio_glucose_TurNuP.mat', 'SolutionEKcat', 'KOresultsE1used');
load('../Models/eciML1515u_TurNuP_CORAL_Ratio.mat');

% load('sftp://mauricio@200.235.253.50/home/mauricio/mauricio_HDnovo/underground_metabolism/DKO_ratio_E1_glucose.mat', 'KOresultsE1used', 'DKOresultsE1used')

%% Calculate mean and median basal flux-sums for single KOs
tic

fluxsums = {};

for ko = 1:length(KOresultsE1used)
        fsKo1Ko2 = abs(model.S) * KOresultsE1used{ko,4};
        fsKo1Ko2 = fsKo1Ko2 * 0.5;
        fsKo1Ko2 = abs(fsKo1Ko2);

        % intervals = {[1834, 2847], [2892, 10153], [10296, 11821]};
        intervals = {[1834, 2847], [2854, 4378], [4417, 11676]};

        idx = true(size(fsKo1Ko2, 1), 1);

        for i = 1:length(intervals)
            idx(intervals{i}(1):intervals{i}(2)) = false;
        end
        
        fsKo1Ko2_filtered = fsKo1Ko2(idx);

        fluxsums = [fluxsums fsKo1Ko2_filtered];

        % index = (ko1 - 1) * length(DKOresultsE1used{ko1,1}) + ko2;

        % fluxsums(:,index) = num2cell(fsKo1Ko2(:,1));
 end

fluxsums = fluxsums';

FSmatrix = {};

for n = 1:length(fsKo1Ko2_filtered)
    valuesAtIndexN = zeros(length(fluxsums), 1);

    for i = 1:length(fluxsums)
        valuesAtIndexN(i) = fluxsums{i}(n);
    end

    % FSmatrix{:,n} = num2cell(valuesAtIndexN);
    FSmatrix = [FSmatrix num2cell(valuesAtIndexN)];

    % Calculate average and median for index n
    % averages(n,1) = mean(valuesAtIndexN);
    % medians(n,1) = median(valuesAtIndexN);
end

mets = model.mets(idx);

% save("FSmatrix.mat","FSmatrix","-v7.3")
toc

%% Calculate mean and median basal flux-sums for all double KO pairs
% tic
% 
% fluxsums = {};
% 
% for ko1 = 1:length(DKOresultsE1used)
%     for ko2 = 1:length(DKOresultsE1used{ko1,1})
%         fsKo1Ko2 = abs(model.S) * DKOresultsE1used{ko1,1}{ko2,4};
%         fsKo1Ko2 = fsKo1Ko2 * 0.5;
%         fsKo1Ko2 = abs(fsKo1Ko2);
% 
%         intervals = {[1834, 2847], [2892, 10153], [10296, 11821]};
% 
%         idx = true(size(fsKo1Ko2, 1), 1);
% 
%         for i = 1:length(intervals)
%             idx(intervals{i}(1):intervals{i}(2)) = false;
%         end
% 
%         fsKo1Ko2_filtered = fsKo1Ko2(idx);
% 
%         fluxsums = [fluxsums fsKo1Ko2_filtered];
% 
%         % index = (ko1 - 1) * length(DKOresultsE1used{ko1,1}) + ko2;
% 
%         % fluxsums(:,index) = num2cell(fsKo1Ko2(:,1));
%     end
% end
% 
% fluxsums = fluxsums';
% 
% FSmatrix = {};
% 
% for n = 1:length(fsKo1Ko2_filtered)
%     valuesAtIndexN = zeros(length(fluxsums), 1);
% 
%     for i = 1:length(fluxsums)
%         valuesAtIndexN(i) = fluxsums{i}(n);
%     end
% 
%     % FSmatrix{:,n} = num2cell(valuesAtIndexN);
%     FSmatrix = [FSmatrix num2cell(valuesAtIndexN)];
% 
%     % Calculate average and median for index n
%     % averages(n,1) = mean(valuesAtIndexN);
%     % medians(n,1) = median(valuesAtIndexN);
% end
% 
% mets = model.mets(idx);
% 
% % save("FSmatrix.mat","FSmatrix","-v7.3")
% toc

%%
% intervals = {[1834, 2847], [2892, 10153], [10296, 11821]};
intervals = {[1834, 2847], [2854, 4378], [4417, 11676]};

idxMet = true(size(model.mets, 1), 1);

for i = 1:length(intervals)
    idxMet(intervals{i}(1):intervals{i}(2)) = false;
end

modelMetsFiltered = model.mets(idxMet);

%% Export the table
filename = "FSmatrix_SKO_glucose_TurNuP.csv";
writecell(FSmatrix, filename, 'Delimiter','\t')

%%
% FSmatrix = readmatrix("FSmatrix.csv");

FSmatrix_filtered = FSmatrix;
FSmatrix_filtered = cell2mat(FSmatrix_filtered);
% FSmatrix_filtered = FSmatrix(:,any(FSmatrix_filtered));

% Remove values lower than 100
max_values = max(FSmatrix_filtered);
columns_to_remove = max_values == 0;
FSmatrix_filtered = FSmatrix_filtered(:, ~columns_to_remove);
modelMetsFiltered = modelMetsFiltered(~columns_to_remove,:);

% column_medians = median(FSmatrix_filtered);
% columns_0 = column_medians ~= 0;
% FSmatrix_filtered = FSmatrix_filtered(:, ~columns_0);
% 
% column_medians = median(FSmatrix_filtered);
% columns_1000 = column_medians == 1000;
% FSmatrix_filtered = FSmatrix_filtered(:, ~columns_1000);

%%
boxplot(FSmatrix_Ratio, 'PlotStyle','compact', 'Labels', modelMetsRatio);
% boxplot(FSmatrix_filtered, 'Labels', modelMetsFiltered);
ylabel('Basal flux-sums');
hLabels = findobj(gca, 'Type', 'text');
set(hLabels, 'HorizontalAlignment', 'right');

%% Calculate basal flux-sums for selected pairs

fluxSum_EKcat = abs(model.S) * SolutionEKcat.full;
fluxSum_EKcat = fluxSum_EKcat * 0.5;
fluxSum_EKcat = abs(fluxSum_EKcat);

fluxSum_EKcat = fluxSum_EKcat';
fluxSum_EKcat = fluxSum_EKcat(idx);
fluxSum_EKcat = fluxSum_EKcat';


% fluxSum_P0AFG8_P0ADG7 = abs(model.S) * DKOresultsE1used{99,1}{79,4};
% fluxSum_P0AFG8_P0ADG7 = fluxSum_P0AFG8_P0ADG7 * 0.5;
% fluxSum_P0AFG8_P0ADG7 = abs(fluxSum_P0AFG8_P0ADG7);
% 
% fluxSum_P0ABE9_P0A9U8 = abs(model.S) * DKOresultsE1used{61,1}{55,4};
% fluxSum_P0ABE9_P0A9U8 = fluxSum_P0ABE9_P0A9U8 * 0.5;
% fluxSum_P0ABE9_P0A9U8 = abs(fluxSum_P0ABE9_P0A9U8);
% 
% fluxSum_P0AGB0_P0AD61 = abs(model.S) * DKOresultsE1used{105,1}{73,4};
% fluxSum_P0AGB0_P0AD61 = fluxSum_P0AGB0_P0AD61 * 0.5;
% fluxSum_P0AGB0_P0AD61 = abs(fluxSum_P0AGB0_P0AD61);
% 
% fluxSum_P23721_P62707 = abs(model.S) * DKOresultsE1used{116,1}{147,4};
% fluxSum_P23721_P62707 = fluxSum_P23721_P62707 * 0.5;
% fluxSum_P23721_P62707 = abs(fluxSum_P23721_P62707);
% 
% fluxSum_P09831_P0ABJ1 = abs(model.S) * DKOresultsE1used{21,1}{64,4};
% fluxSum_P09831_P0ABJ1 = fluxSum_P09831_P0ABJ1 * 0.5;
% fluxSum_P09831_P0ABJ1 = abs(fluxSum_P09831_P0ABJ1);
% 
% %% Construct the output table
% FStable = table();
% FStable.mets = model.mets;
% FStable.metNames = model.metNames;
% 
% FStable.average = averages;
% FStable.median = medians;
% 
% FStable.P0AFG8_P0ADG7 = fluxSum_P0AFG8_P0ADG7;
% FStable.P0ABE9_P0A9U8 = fluxSum_P0ABE9_P0A9U8;
% FStable.P0AGB0_P0AD61 = fluxSum_P0AGB0_P0AD61;
% FStable.P23721_P62707 = fluxSum_P23721_P62707;
% FStable.P09831_P0ABJ1 = fluxSum_P09831_P0ABJ1;
% 
% patternToRemove1 = 'prot_';
% rowsToRemove = contains(FStable.mets, patternToRemove1);
% FStable(rowsToRemove, :) = [];
% patternToRemove2 = 'pmet_';
% rowsToRemove = contains(FStable.mets, patternToRemove2);
% FStable(rowsToRemove, :) = [];
% patternToRemove3 = 'subpool_';
% rowsToRemove = contains(FStable.mets, patternToRemove3);
% FStable(rowsToRemove, :) = [];
% patternToRemove3 = 'subpool_';
% rowsToRemove = contains(FStable.mets, patternToRemove3);
% FStable(rowsToRemove, :) = [];
% 
% %% Plot the first figure
% 
% tiledlayout(1,2)
% 
% nexttile
% hold on
% stem(FStable.average);
% title('Mean basal flux-sums')
% xlabel('Metabolite')
% ylabel('Flux-sum (mmol/gCDW·h')
% hold off
% 
% nexttile
% hold on
% stem(FStable.median);
% title('Median basal flux-sums')
% xlabel('Metabolite')
% hold off
% 
% %% Plot the second figure
% 
% tiledlayout(1,5)
% 
% nexttile
% hold on
% stem(FStable.P0AFG8_P0ADG7);
% title('P0AFG8 and P0ADG7')
% xlabel('Metabolite')
% ylabel('Flux-sum (mmol/gCDW·h')
% hold off
% 
% nexttile
% hold on
% stem(FStable.P0ABE9_P0A9U8);
% title('P0ABE9 and P0A9U8')
% xlabel('Metabolite')
% hold off
% 
% nexttile
% hold on
% stem(FStable.P0AGB0_P0AD61);
% title('P0AGB0 and P0AD61')
% xlabel('Metabolite')
% hold off
% 
% nexttile
% hold on
% stem(FStable.P23721_P62707);
% title('P23721 and P62707')
% xlabel('Metabolite')
% hold off
% 
% nexttile
% hold on
% stem(FStable.P09831_P0ABJ1);
% title('P09831 and P0ABJ1')
% xlabel('Metabolite')
% hold off
% 
% %% Export the table
% filename = "DKO_ratio_E1_glucose_basalfluxsums.csv";
% writetable(FStable, filename, 'Delimiter','\t')

