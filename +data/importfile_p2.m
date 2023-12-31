function NOKNET3DPhase2CalibrationGeogrDistancev59S1 = importfile_p2(workbookFile, sheetName, dataLines)
%IMPORTFILE 
%  NOKNET3DPHASE2CALIBRATIONGEOGRDISTANCEV59S1 = IMPORTFILE(FILE) 
%  FILE is Microsoft Excel 
%
%  NOKNET3DPHASE2CALIBRATIONGEOGRDISTANCEV59S1 = IMPORTFILE(FILE, SHEET)
%
%  NOKNET3DPHASE2CALIBRATIONGEOGRDISTANCEV59S1 = IMPORTFILE(FILE, SHEET, DATALINES)
%
%  NOKNET3DPhase2CalibrationGeogrDistancev59S1 = data.importfile(".\data_3GPP\NOKNET_3DPhase2CalibrationGeogrDistance_v59.xls", "3D-UMi (K=1,M=2)", [29, 129]);
%
%  refer to READTABLEĦ£

    %% deal with input
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    if nargin <= 2
        dataLines = [29, 129];
    end
    opts = spreadsheetImportOptions("NumVariables", 221);

    % range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + dataLines(1, 1) + ":HM" + dataLines(1, 2);

    % name and type
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "VarName31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "VarName62", "Var63", "Var64", "Var65", "Var66", "Var67", "Var68", "Var69", "Var70", "Var71", "Var72", "Var73", "Var74", "Var75", "Var76", "Var77", "Var78", "Var79", "Var80", "Var81", "Var82", "Var83", "Var84", "Var85", "Var86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "VarName93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102", "Var103", "Var104", "Var105", "Var106", "Var107", "Var108", "Var109", "Var110", "Var111", "Var112", "Var113", "Var114", "Var115", "Var116", "Var117", "Var118", "Var119", "Var120", "Var121", "Var122", "Var123", "Var124", "VarName125", "Var126", "Var127", "Var128", "Var129", "Var130", "Var131", "Var132", "Var133", "Var134", "Var135", "Var136", "Var137", "Var138", "Var139", "Var140", "Var141", "Var142", "Var143", "Var144", "Var145", "Var146", "Var147", "Var148", "Var149", "Var150", "Var151", "Var152", "Var153", "Var154", "Var155", "Var156", "VarName157", "Var158", "Var159", "Var160", "Var161", "Var162", "Var163", "Var164", "Var165", "Var166", "Var167", "Var168", "Var169", "Var170", "Var171", "Var172", "Var173", "Var174", "Var175", "Var176", "Var177", "Var178", "Var179", "Var180", "Var181", "Var182", "Var183", "Var184", "Var185", "Var186", "Var187", "Var188", "VarName189", "Var190", "Var191", "Var192", "Var193", "Var194", "Var195", "Var196", "Var197", "Var198", "Var199", "Var200", "Var201", "Var202", "Var203", "Var204", "Var205", "Var206", "Var207", "Var208", "Var209", "Var210", "Var211", "Var212", "Var213", "Var214", "Var215", "Var216", "Var217", "Var218", "Var219", "Var220", "VarName221"];
    opts.SelectedVariableNames = ["VarName31", "VarName62", "VarName93", "VarName125", "VarName157", "VarName189", "VarName221"];
    opts.VariableTypes = ["char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double"];
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220], "EmptyFieldRule", "auto");

    % import data
    NOKNET3DPhase2CalibrationGeogrDistancev59S1 = readtable(workbookFile, opts, "UseExcel", false);

    for idx = 2:size(dataLines, 1)
        opts.DataRange = "A" + dataLines(idx, 1) + ":HM" + dataLines(idx, 2);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        NOKNET3DPhase2CalibrationGeogrDistancev59S1 = [NOKNET3DPhase2CalibrationGeogrDistancev59S1; tb]; %#ok<AGROW>
    end
        NOKNET3DPhase2CalibrationGeogrDistancev59S1 = table2array(NOKNET3DPhase2CalibrationGeogrDistancev59S1);
end