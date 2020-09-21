function Table = load_burn_excel()

[fileName, pathName] = uigetfile('*.xlsx');
Table = readtable(fullfile(pathName,fileName),'Sheet',2);
end

