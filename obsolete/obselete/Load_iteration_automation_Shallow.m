%% Load_iteration_automation
clc; close all; clear; 
addpath (genpath('library')); 
%% Settings

manual_file_name            = 'manual_data_input';
location                    = 'Deep';
load_iteration              = 'L11';
analysis                    = 'Load5'; % not really needed
load_excel_name             = 'L0SG25_11_MW';
load_excel_link             = '\\COWI.net\projects\A210000\A214109\20-Data\GEO_STR\';


%% Analyses

lateral_UR                  = 0;
lateral_deflection          = 0;
permanent_rot_100           = 0;
permanent_rot_10000         = 0;
axial_UR                    = 1;

%% Read excels

[loads, names, DLC] = extractloads_fkmv(load_iteration, load_excel_name, load_excel_link);

%% Run COSPIN functions
for i = 2:size(names,1)
    % uncomment these overwrites if you want to play around with the pile L and D
%     DLC.pile_L(i,:) = 28; % [m]
%     DLC.pile_D(i,:) = 9;  % [m]
    
    ID = strcat(names{i,1},names{i,3});
    save_name = strcat(names{i,2},names{i,3});

    if  permanent_rot_100 == 1 && strcmp(names{i,3},'enviromental ')
        [per_rot.output] = run_COSPIN_perm_rot_100(ID, location, manual_file_name, load_iteration, DLC.load_case{i}, analysis, DLC.pile_D(i), DLC.pile_L(i), DLC, names, save_name, i);
    end
    if  permanent_rot_10000 == 1 && strcmp(names{i,3},'enviromental ')
        [per_rot.output] = run_COSPIN_perm_rot_10000(ID, location, manual_file_name, load_iteration, DLC.load_case{i}, analysis, DLC.pile_D(i), DLC.pile_L(i), DLC, names, save_name, i);
    end
    if lateral_UR == 1
        [load_def.output] = run_COSPIN_load_deflection(ID, location, manual_file_name, load_iteration, DLC.load_case{i}, analysis, DLC.pile_D(i), DLC.pile_L(i), DLC, names, save_name, i);
    end
    if lateral_deflection == 1
       [lateral.output] = run_COSPIN_lateral_deflection(ID, location, manual_file_name, load_iteration, DLC.load_case{i}, analysis, DLC.pile_D(i), DLC.pile_L(i), DLC, names, save_name, i); 
    end
    if axial_UR == 1
        [per_rot.output] = run_COSPIN_axial_capacity(ID, location, manual_file_name, load_iteration, DLC.load_case{i}, analysis, DLC.pile_D(i), DLC.pile_L(i), DLC, names, save_name, i);
    end

close all;
end

%% Store results in excel