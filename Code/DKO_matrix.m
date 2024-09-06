%%
clc;clear

%%
DKOresultsE1 = readtable("../Results/DKO_ratio_glucose_E1_TurNuP.csv");
DKOresultsE2 = readtable("../Results/DKO_ratio_xylose_E2.csv");
DKOresultsE3 = readtable("../Results/DKO_ratio_xylose_E3.csv");
DKOresultsE4 = readtable("../Results/DKO_ratio_xylose_E4.csv");
DKOresultsE5 = readtable("../Results/DKO_ratio_xylose_E5.csv");

%%
% E1
DKOresultsE1.ratio = DKOresultsE1.DKOMu ./ 0.09;
proteinNamesE1 = unique([DKOresultsE1.FirstKO(:, 1); DKOresultsE1.DKO(:, 1)]);
numProteinsE1 = numel(proteinNamesE1);
matrixE1 = zeros(numProteinsE1);

for i = 1:size(DKOresultsE1, 1)
    rowIdx = find(strcmp(proteinNamesE1, DKOresultsE1.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNamesE1, DKOresultsE1.DKO(i, 1)));
    matrixE1(rowIdx, colIdx) = DKOresultsE1.ratio(i, 1);
end

% E2
DKOresultsE2.ratio = DKOresultsE2.DKOMu ./ 0.09;
proteinNamesE2 = unique([DKOresultsE2.FirstKO(:, 1); DKOresultsE2.DKO(:, 1)]);
numProteinsE2 = numel(proteinNamesE2);
matrixE2 = zeros(numProteinsE2);

for i = 1:size(DKOresultsE2, 1)
    rowIdx = find(strcmp(proteinNamesE2, DKOresultsE2.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNamesE2, DKOresultsE2.DKO(i, 1)));
    matrixE2(rowIdx, colIdx) = DKOresultsE2.ratio(i, 1);
end

% E3
DKOresultsE3.ratio = DKOresultsE3.DKOMu ./ 0.09;
proteinNamesE3 = unique([DKOresultsE3.FirstKO(:, 1); DKOresultsE3.DKO(:, 1)]);
numProteinsE3 = numel(proteinNamesE3);
matrixE3 = zeros(numProteinsE3);

for i = 1:size(DKOresultsE3, 1)
    rowIdx = find(strcmp(proteinNamesE3, DKOresultsE3.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNamesE3, DKOresultsE3.DKO(i, 1)));
    matrixE3(rowIdx, colIdx) = DKOresultsE3.ratio(i, 1);
end

% E4
DKOresultsE4.ratio = DKOresultsE4.DKOMu ./ 0.09;
proteinNamesE4 = unique([DKOresultsE4.FirstKO(:, 1); DKOresultsE4.DKO(:, 1)]);
numProteinsE4 = numel(proteinNamesE4);
matrixE4 = zeros(numProteinsE4);

for i = 1:size(DKOresultsE4, 1)
    rowIdx = find(strcmp(proteinNamesE4, DKOresultsE4.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNamesE4, DKOresultsE4.DKO(i, 1)));
    matrixE4(rowIdx, colIdx) = DKOresultsE4.ratio(i, 1);
end

% E5
DKOresultsE5.ratio = DKOresultsE5.DKOMu ./ 0.09;
proteinNamesE5 = unique([DKOresultsE5.FirstKO(:, 1); DKOresultsE5.DKO(:, 1)]);
numProteinsE5 = numel(proteinNamesE5);
matrixE5 = zeros(numProteinsE5);

for i = 1:size(DKOresultsE5, 1)
    rowIdx = find(strcmp(proteinNamesE5, DKOresultsE5.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNamesE5, DKOresultsE5.DKO(i, 1)));
    matrixE5(rowIdx, colIdx) = DKOresultsE5.ratio(i, 1);
end

%%
tiledlayout(1,5)
clims = [0 1];

nexttile
hold on
imagesc(matrixE1)
xlabel('E_s_,_1 knockout')
ylabel('E_j_,_1 knockout')
clim(clims);
hold off

nexttile
hold on
imagesc(matrixE2)
xlabel('E_s_,_2 knockout')
ylabel('E_j_,_2 knockout')
clim(clims);
hold off

nexttile
hold on
imagesc(matrixE3)
xlabel('E_s_,_3 knockout')
ylabel('E_j_,_3 knockout')
clim(clims);
hold off

nexttile
hold on
imagesc(matrixE4)
xlabel('E_s_,_4 knockout')
ylabel('E_j_,_4 knockout')
clim(clims);
hold off

nexttile
hold on
imagesc(matrixE5)
xlabel('E_s_,_5 knockout')
ylabel('E_j_,_5 knockout')
clim(clims);
hold off