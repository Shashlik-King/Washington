function AppendixGeneration(Settings,AppendixData)


%loc=Settings.Locations(any(cellfun(@(x)any(~isnan(x)),Settings.Locations(:,5)),2),:);

%    Data=GE_Data.(A.SimulationLable);
disp([' Creating Appendix for Location: ' , AppendixData.locationID]) 
% Appendix function for generating an appendix file for each location
path_LatexFolder = [pwd '\AppendixGenerationFiles\ProjectLocation'];

%% Generate files for the site specific information
% Generate location specific information
AppendixGeneration_SiteInput(Settings)

%% Generating appendix pdf
% Generate appendix for each location considered in defined run
% for c=1:size(ID,1) FKMV
    % Delete the files current located in the Location specific folder (ProjectLocation)
    dinfo = dir(path_LatexFolder);
    dinfo([dinfo.isdir]) = [];          % Skip directories
%     if ~isempty(dinfo)
%         filenames = fullfile(path_LatexFolder, {dinfo.name});
%         delete( filenames{:} )
%     end
    
    % Generate all .tex files for specific location considered
    %AppendixGeneration_LocationInput(results.ID{c},Settings, results.ID{c}.def.data.location, results.ID{c}.def.pile.length,PileGeometry)
    AppendixGeneration_LocationInput(Settings,AppendixData)
    % Copy and rename files to latex folder
%     dinfo = dir(fullfile(pwd,'Plots'));
%     counter = 0;
%     for i = 1:length(dinfo)
%         if length(dinfo(i).name) > length(char(Data.loc{c,1})) && strcmp(dinfo(i).name(1:length(char(Data.loc{c,1}))+1),{[char(Data.loc{c,1}) '_']})
%             counter = counter + 1;
%             names{counter} =  dinfo(i).name;
%             names_latex{counter} = dinfo(i).name(length(char(Data.loc{c,1}))+2:end);
%             %copyfile(fullfile(pwd,A.Folder,['Plots\' names{counter}]), fullfile(pwd,['AppendixGenerationFiles\ProjectLocation\' names_latex{counter}]));
%             copyfile(fullfile(pwd,['Plots\' names{counter}]), fullfile(path_LatexFolder,names_latex{counter}));
%         else
%         end
%     end
%     
    % Run LaTeX to create a pdf
    cmd = {'pdflatex AppendixTemplate.tex'};   % Command for running latex code and create a pdf
    [unix_1] = unix(char(cmd));
    pause(1)
    [unix_2] = unix(char(cmd));
    
    % Move the pdf file to storage of pdf documentation
    movefile(char(fullfile(pwd, 'AppendixTemplate.pdf')), char(fullfile(pwd,{['Appendix\AppendixGEO_' char(AppendixData.locationID) '.pdf']})));
    
    % Clear generation files from main folder
    dinfo_path = dir(fullfile(pwd));
    dinfo_path([dinfo_path.isdir]) = [];
    for i = 1:length(dinfo_path)
        if contains(dinfo_path(i).name,{'AppendixTemplate'}) && ~strcmp(dinfo_path(i).name,{'AppendixTemplate.tex'})
            delete(dinfo_path(i).name)
        end
    end
% end
clc;    % Clear command window after files are generated