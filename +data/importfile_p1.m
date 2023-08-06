function DPhase1CalibrationGeogrDistancev20S1 = importfile_p1(workbookFile, sheetName, dataLines)
%IMPORTFILE
%  DPHASE1CALIBRATIONGEOGRDISTANCEV20S1 = IMPORTFILE(FILE)
%  Microsoft Excel
%
%  DPHASE1CALIBRATIONGEOGRDISTANCEV20S1 = IMPORTFILE(FILE, SHEET)
%  
%
%  DPHASE1CALIBRATIONGEOGRDISTANCEV20S1 = IMPORTFILE(FILE, SHEET, DATALINES)
%
%  DPhase1CalibrationGeogrDistancev20S1 = data.importfile(".\data_3GPP\3DPhase1CalibrationGeogrDistance_v20.xlsx", "3D-UMi (K=M=1)", [29, 129]);
%
%  refer to READTABLE¡£

    %% deal with input
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    if nargin <= 2
        dataLines = [29, 129];
    end
    opts = spreadsheetImportOptions("NumVariables", 93);

    % range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + dataLines(1, 1) + ":CO" + dataLines(1, 2);

    % name and type
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "VarName31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "VarName62", "Var63", "Var64", "Var65", "Var66", "Var67", "Var68", "Var69", "Var70", "Var71", "Var72", "Var73", "Var74", "Var75", "Var76", "Var77", "Var78", "Var79", "Var80", "Var81", "Var82", "Var83", "Var84", "Var85", "Var86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "VarName93"];
    opts.SelectedVariableNames = ["VarName31", "VarName62", "VarName93"];
    opts.VariableTypes = ["char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double"];
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92], "EmptyFieldRule", "auto");

    % import data
    DPhase1CalibrationGeogrDistancev20S1 = readtable(workbookFile, opts, "UseExcel", false);

    for idx = 2:size(dataLines, 1)
        opts.DataRange = "A" + dataLines(idx, 1) + ":CO" + dataLines(idx, 2);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        DPhase1CalibrationGeogrDistancev20S1 = [DPhase1CalibrationGeogrDistancev20S1; tb]; %#ok<AGROW>
    end
        DPhase1CalibrationGeogrDistancev20S1 = table2array(DPhase1CalibrationGeogrDistancev20S1);
end