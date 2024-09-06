%
% Plot enzyme usage FVA
%
%%
clc;clear

%%
noBio_noUnd = readtable("../Results/FVA_enzymes_eciML1515u_CORAL_TurNuP_noBio_noUnd_Ratio.csv");
noBio_Und = readtable("../Results/FVA_enzymes_eciML1515u_CORAL_TurNuP_noBio_Und_Ratio.csv");
Bio_noUnd = readtable("../Results/FVA_enzymes_eciML1515u_CORAL_TurNuP_Bio_noUnd_Ratio.csv");
Bio_Und = readtable("../Results/FVA_enzymes_eciML1515u_CORAL_TurNuP_Bio_Und_Ratio.csv");

%%
distributions = {noBio_noUnd.ranges, noBio_Und.ranges, Bio_noUnd.ranges, Bio_Und.ranges};
legends       = {'No underground reactions, no fixed biomass', 'Underground reactions, no fixed biomass','No underground reactions, fixed biomass', 'Underground reactions, fixed biomass'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);
