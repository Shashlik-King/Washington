function MoveFiles(ID)
outputFolder = '\AppendixGenerationFiles\ProjectLocation';
outputFullFileName = fullfile(pwd,outputFolder, 'NFA_1.png');
inputFolder = '\output';
inputFiles = ['natural_frequency_',ID,'.png'];
% inputFiles = ['NFA_2_',ID{i},'.png']
% inputFiles = ['soil_properties_',ID{i},'.png']
inputFullFileName = fullfile(pwd,inputFolder,inputFiles);
copyfile(inputFullFileName,outputFullFileName);

outputFolder = '\AppendixGenerationFiles\ProjectLocation';
outputFullFileName = fullfile(pwd,outputFolder, 'NFA_2.png');
inputFolder = '\output';
% inputFiles = ['NFA_1_',ID{i},'.png']
inputFiles = ['relative_alpha_',ID,'.png'];
% inputFiles = ['soil_properties_',ID{i},'.png']
inputFullFileName = fullfile(pwd,inputFolder,inputFiles);
copyfile(inputFullFileName,outputFullFileName);

outputFolder = '\AppendixGenerationFiles\ProjectLocation';
outputFullFileName = fullfile(pwd,outputFolder, 'NFA_3.png');
inputFolder = '\output';
% inputFiles = ['NFA_1_',ID{i},'.png']
inputFiles = ['delta_alpha_',ID,'.png'];
% inputFiles = ['soil_properties_',ID{i},'.png']
inputFullFileName = fullfile(pwd,inputFolder,inputFiles);
copyfile(inputFullFileName,outputFullFileName);

outputFolder = '\AppendixGenerationFiles\ProjectLocation';
outputFullFileName = fullfile(pwd,outputFolder, 'soil_properties.png');
inputFolder = '\output';
inputFiles = ['soil_properties_',ID,'.png'];
inputFullFileName = fullfile(pwd,inputFolder,inputFiles);
copyfile(inputFullFileName,outputFullFileName);