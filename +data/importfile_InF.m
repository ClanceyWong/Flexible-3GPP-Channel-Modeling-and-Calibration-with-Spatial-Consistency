function R11909704IndoorindustrialchannelmodelcalibrationresultsS1 = importfile(workbookFile, sheetName, dataLines)
%IMPORTFILE 导入电子表格中的数据
%  R11909704INDOORINDUSTRIALCHANNELMODELCALIBRATIONRESULTSS1 =
%  IMPORTFILE(FILE) 读取名为 FILE 的 Microsoft Excel 电子表格文件的第一张工作表中的数据。
%  以表形式返回数据。
%
%  R11909704INDOORINDUSTRIALCHANNELMODELCALIBRATIONRESULTSS1 =
%  IMPORTFILE(FILE, SHEET) 从指定的工作表中读取。
%
%  R11909704INDOORINDUSTRIALCHANNELMODELCALIBRATIONRESULTSS1 =
%  IMPORTFILE(FILE, SHEET, DATALINES)按指定的行间隔读取指定工作表中的数据。对于不连续的行间隔，请将
%  DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  R11909704IndoorindustrialchannelmodelcalibrationresultsS1 = importfile("E:\channel_3GPP\Docs\R1-1909704_InF_calibration\R1-1909704 Indoor industrial channel model calibration results.xlsx", "Sub-scenario 1 - 3.5 GHz", [3, 103]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2020-09-17 20:00:14 自动生成

%% 输入处理

% 如果未指定工作表，则将读取第一张工作表
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% 如果未指定行的起点和终点，则会定义默认值。
if nargin <= 2
    dataLines = [3, 103];
end

%% 设置导入选项并导入数据
opts = spreadsheetImportOptions("NumVariables", 62);

% 指定工作表和范围
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":BJ" + dataLines(1, 2);

% 指定列名称和类型
opts.VariableNames = ["Percentile", "Ericsson", "Nokia", "ZTE", "Huawei", "Mean", "VarName7", "Percentile1", "Ericsson1", "Nokia1", "ZTE1", "Huawei1", "Mean1", "VarName14", "Percentile2", "Ericsson2", "Nokia2", "ZTE2", "Huawei2", "Mean2", "VarName21", "Percentile3", "Ericsson3", "Nokia3", "ZTE3", "Huawei3", "Mean3", "VarName28", "Percentile4", "Ericsson4", "Nokia4", "ZTE4", "Huawei4", "Mean4", "VarName35", "Percentile5", "Ericsson5", "Nokia5", "ZTE5", "Huawei5", "Mean5", "VarName42", "Percentile6", "Ericsson6", "Nokia6", "ZTE6", "Huawei6", "Mean6", "VarName49", "Percentile7", "Ericsson7", "Nokia7", "ZTE7", "Huawei7", "Mean7", "VarName56", "Percentile8", "Company1", "Company2", "VarName60", "VarName61", "Mean8"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "string", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "string", "double", "string", "string", "string", "string", "string"];

% 指定变量属性
opts = setvaropts(opts, ["VarName7", "VarName14", "Huawei2", "VarName21", "VarName28", "VarName35", "VarName42", "VarName49", "VarName56", "Company1", "Company2", "VarName60", "VarName61", "Mean8"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName7", "VarName14", "Huawei2", "VarName21", "VarName28", "VarName35", "VarName42", "VarName49", "VarName56", "Company1", "Company2", "VarName60", "VarName61", "Mean8"], "EmptyFieldRule", "auto");

% 导入数据
R11909704IndoorindustrialchannelmodelcalibrationresultsS1 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":BJ" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    R11909704IndoorindustrialchannelmodelcalibrationresultsS1 = [R11909704IndoorindustrialchannelmodelcalibrationresultsS1; tb]; %#ok<AGROW>
end
R11909704IndoorindustrialchannelmodelcalibrationresultsS1 = str2double(table2array(R11909704IndoorindustrialchannelmodelcalibrationresultsS1));
end