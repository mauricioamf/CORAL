%%
clc;clear

%%
DKOresults1 = readtable("../Results/DKO_GEM_quarter1.csv");
DKOresults2 = readtable("../Results/DKO_GEM_quarter2.csv");
DKOresults3 = readtable("../Results/DKO_GEM_quarter3.csv");
DKOresults4 = readtable("../Results/DKO_GEM_quarter4.csv");

DKOresults1 = removevars(DKOresults1, "Ratio");
DKOresults2 = removevars(DKOresults2, "Ratio");
DKOresults3 = removevars(DKOresults3, "Ratio");
DKOresults4 = removevars(DKOresults4, "Ratio");

DKOresults = [DKOresults1; DKOresults2; DKOresults3; DKOresults4];

%%
DKOresults.ratio = DKOresults.DKOMu ./ 0.815;
proteinNames = unique([DKOresults.FirstKO(:, 1); DKOresults.DKO(:, 1)]);
numProteins = numel(proteinNames);
matrix = zeros(numProteins);

for i = 1:size(DKOresults, 1)
    rowIdx = find(strcmp(proteinNames, DKOresults.FirstKO(i, 1)));
    colIdx = find(strcmp(proteinNames, DKOresults.DKO(i, 1)));
    matrix(rowIdx, colIdx) = DKOresults.ratio(i, 1);
end

%%
hold on
imagesc(matrix)
xlabel('Gene 1 knockout')
ylabel('Gene 2 knockout')
hold off