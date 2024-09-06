%%
clear;clc

%%
iJO1366_combined = readCbModel('combined.xml');
writeCbModel(iJO1366_combined, 'iJO1366_combined.mat');
iJO1366_combined = ravenCobraWrapper(iJO1366_combined);
exportToExcelFormat(iJO1366_combined, 'iJO1366_combined.xlsx');
clear;clc

%%
iJO1366_heterologous = readCbModel('heterologous.xml');
writeCbModel(iJO1366_heterologous, 'iJO1366_heterologous.mat');
writeCbModel(iJO1366_heterologous, 'format', 'xls', 'fileName', 'iJO1366_heterologous.xlsx');
% iJO1366_heterologous = ravenCobraWrapper(iJO1366_heterologous);
% exportToExcelFormat(iJO1366_heterologous, 'iJO1366_heterologous.xlsx');
clear;clc

%%
iJO1366_underground = readCbModel('underground.xml');
writeCbModel(iJO1366_underground, 'iJO1366_underground.mat');
iJO1366_underground = ravenCobraWrapper(iJO1366_underground);
exportToExcelFormat(iJO1366_underground, 'iJO1366_underground.xlsx');
clear;clc

%%
iML1515 = readCbModel('iML1515.xml');
writeCbModel(iML1515, 'iML1515.mat');
iML1515 = ravenCobraWrapper(iML1515);
exportToExcelFormat(iML1515, 'iML1515.xlsx');
clear;clc

%%
iML1515_adjusted = readCbModel('iML1515_model_adjusted_all_corrections.xml');
writeCbModel(iML1515_adjusted, 'iML1515_adjusted.mat');
iML1515_adjusted = ravenCobraWrapper(iML1515_adjusted);
exportToExcelFormat(iML1515_adjusted, 'iML1515_adjusted.xlsx');
clear;clc

%%
disp('All done')