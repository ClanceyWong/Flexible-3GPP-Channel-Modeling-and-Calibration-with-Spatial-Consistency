function Phase1Calibrationv42CMCCS1 = importfile_p1_all(workbookFile, sheetName, dataLines)
%IMPORTFILE 导入电子表格中的数据
%  PHASE1CALIBRATIONV42CMCCS1 = IMPORTFILE(FILE) 读取名为 FILE 的 Microsoft
%  Excel 电子表格文件的第一张工作表中的数据。  以表形式返回数据。
%
%  PHASE1CALIBRATIONV42CMCCS1 = IMPORTFILE(FILE, SHEET) 从指定的工作表中读取。
%
%  PHASE1CALIBRATIONV42CMCCS1 = IMPORTFILE(FILE, SHEET,
%  DATALINES)按指定的行间隔读取指定工作表中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2
%  正整数标量数组。
%
%  示例:
%  Phase1Calibrationv42CMCCS1 = importfile("E:\channel_3GPP\Docs\R1-165974_large_scale_calibration\Phase1Calibration_v42_CMCC.xlsx", "UMi-6GHz", [25, 256]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2020-11-06 10:24:09 自动生成

%% 输入处理

% 如果未指定工作表，则将读取第一张工作表
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% 如果未指定行的起点和终点，则会定义默认值。
if nargin <= 2
    dataLines = [29, 129];
end

%% 设置导入选项并导入数据
opts = spreadsheetImportOptions("NumVariables", 93);

% 指定工作表和范围
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":CO" + dataLines(1, 2);

% 指定列名称和类型
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72", "VarName73", "VarName74", "VarName75", "VarName76", "VarName77", "VarName78", "VarName79", "VarName80", "VarName81", "VarName82", "VarName83", "VarName84", "VarName85", "VarName86", "VarName87", "VarName88", "VarName89", "VarName90", "VarName91", "VarName92", "VarName93"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 指定变量属性
opts = setvaropts(opts, ["VarName32", "VarName63"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName32", "VarName63"], "EmptyFieldRule", "auto");

% 导入数据
Phase1Calibrationv42CMCCS1 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":CO" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    Phase1Calibrationv42CMCCS1 = [Phase1Calibrationv42CMCCS1; tb]; %#ok<AGROW>
end
Phase1Calibrationv42CMCCS1 = table2array(Phase1Calibrationv42CMCCS1);
end