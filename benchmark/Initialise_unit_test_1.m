function [ID,settings,soil,pile,loads,data,plots] = Initialise()
ID = {'benchmark'};
%--------------------------------------------------------------------------
%% Soil and pile input
%--------------------------------------------------------------------------
settings.database               = 'Manual';       % if 'Database' the database module is activated. If 'Manual' input is taken from manual_data_input.xlsx
settings.PISA_database          = 1;              % Switch for PISA database usage (1 = on, 0 = off)
settings.psf_switch = 0;
% interim settings
settings.interimloads           = 0;
settings.interimgeometry        = 0;
soil.type_su                    = 0;              % Soil degradation switch (1 = to use the factors in cyclicDegradation.m , 0 = to not use)

%below only if reading from database is set
loads.type                      = 'ULS';          % Loads to be used (ULS or FLS)maybe also GEO
soil.psf                        = 0;              % Partial Safety Factors 0 reads the _geo columns from the load table in the database 
soil.type                       = 'BE';           % 'BE', 'LB', 'UB'

data.table_springs              = 'char' ;        % 'char', 'LB', 'UB', 'test' Table to update in the database
% WARNING: ONLY FOR SPECIAL USE (only applied if data.revision.global = -1):
data.revision.soil              = 0;              % revision no. of soil parameters to be used (1000 = latest revision)
data.revision.structure         = 0;              % revision no. of structure to be used (1000 = latest revision)
data.revision.loads             = 0;              % revision no. of ULS loads to be used (1000 = latest revision)
data.revision.output            = 0;             % revision no. for storing results into the database  

settings.interface              = 'GEO';          % FLS: FLS loads are applied, ULS: factored ULS loads are applied, 
                                                  % GEO: factored/unfactored loads are applied depending on the check
                                                  % that are carried out
settings.db_server              ='DKLYCOPILOD1';  % Databse server
settings.db_user                ='ewdb_user';        % Database user
settings.db_pass                ='ituotdewdb';  % Database pass
settings.db_name                ='ewdb';     % Database name                                      
%--------------------------------------------------------------------------
%% General calculation settings
%--------------------------------------------------------------------------
settings.nelem_factor           = 1;        % [m-1] minimum number of pile beam elements per meter 2.2
settings.j_max                  = 100;      % max. no. of iterations in a load step (lateral only)
settings.TOL                    = 1e-4;     % relative tolerance (lateral and axial)
settings.n_max                  = 1;        % number of load steps (lateral and axial)
%--------------------------------------------------------------------------
%% Calculation settings (axial)
%--------------------------------------------------------------------------
settings.axial_loading          = 1;        % calculate capacity for axial loading? 1 = yes, 0 = no
settings.clay_type              = 'alpha';  % 'alpha' or 'beta' for way to calculate skin friction in API clay
settings.analysis_type          = 0;        % 0 = normal axial assessment, 1 = settlements
settings.analysis_loading       = {'Comp'}; % {'Tens'} = tension only, {'Comp'} = compression only, {'Tens', 'Comp'} = both 
settings.Delta_lambda_bar       = 0.05;     % inital load increment for AL (Arc length) only used in axial calculation. Controls the number of load steps - higher value = fewer load steps, lower value = more load steps 
settings.axial_NR_AL            = 0;        % 0 = Newton-Raphson, 1 = Arc length IMPORTANT: A-L solver not fully QA-ed. Use not recomended

plots.res_vs_pilelength         = 0;        % 1 = yes, 0 = no
plots.normal_force_axial        = 0;        % Normal force plot. 1 = yes, 0 = no
plots.UR_axial                  = 0;        % utilization ratio for axial capacity
plots.load_deflection_axial		= 0;		% 1 = yes (should not be combined with other plots; minimum 300 load steps are recommended), 0 = no
plots.deflection_plot_axial     = 0;        % Axial displacement along the pile, 0 = no, 1 = yes
plots.py_curve_UR = 0;
pile.plug_unplug.tens           = 0;        % plug/unplug control in tension - 0 = auto, 1 = plugged, 2 = unplugged
pile.plug_unplug.comp           = 0;        % plug/unplug control in compression - 0 = auto, 1 = plugged, 2 = unplugged
pile.extra_L                    = 20;       % [m] extra pile length to be used in axial capacity plot

soil.zpeak                      = 0.094;    % usually = pile.diameter*0.01
soil.zres                       = 2;        % usually = 2
soil.zpeak_Qz                   = 0.94;     % usually = pile.diameter*0.1
soil.tres                       = 0.7;      % Usually = 0.7 - 0.9
%--------------------------------------------------------------------------
%% Calculation settings (torsional)
%--------------------------------------------------------------------------
settings.torsion                = 0;        % calculate torsional stiffness? 1 = yes, 0 = no. settings.axial_loading must also be set to 1
%--------------------------------------------------------------------------
%% Calculation settings (lateral)
%--------------------------------------------------------------------------
loads.static_cyclic             = 'cyclic'; % 'cyclic' or 'static' in accordance with DNV OS-J101
loads.n_cycles                  = 100;      % Number of cycles, relevant for Stiff clay w/o water only!
loads.A                         = 'API';    % use API or TUHH (Dührkop) approach for determination of A, relevant for API/Kirsch sand



if strcmp(soil.type,'UB')
loads.static_cyclic= 'static';
end 




plots.pilehead_vs_length        = 0;        % 1 = yes, 0 = no
plots.deflection_plot           = 1;        % 1 = yes, 0 = no
plots.utilization_ratio         = 0;        % 1 = yes, 0 = no, settings.toe_shear = 0 if only p-y UR is of interest
plots.deflection_bundle         = 0;        % 1 = yes, 0 = no
plots.toe_shear_graph           = 0;        % 1 = yes, 0 = no
plots.permanent_rot_def         = 0;        % 1 = yes (minimum 10 load steps is recommended), 0 = no -- cannot be used together with other analyses
plots.load_deflection			= 0;		% 1 = yes (should not be combined with other plots; minimum 50 load steps us recommended), 0 = no
plots.load_deflection_type      = 0;		% 0 = DNV-Gl / BSH , 1 = regular load-displacement
plots.moment_distribution		= 0;		% 1 = yes, 0 = no

settings.lateral_loading        = 1;        % calculate capacity for lateral loading? 1 = yes, 0 = no
settings.beam_theory            = 1;        % 1 = Timoshenko, 0 = Euler-Bernoulli
settings.toe_shear              = 0;        % include base shear and moment? 1 = yes, 0 = no  
settings.mteta	                = 0;        % include uniformly distributed moment? 1 = yes, 0 = no 
settings.Georgiadis				= 0;		% apply Georgiadis approach? 1 = yes, 0 = no
settings.lateralmultipliers 	= 0; 		% account for p and y multipliers: 0-> not accounted 1-> accounted
settings.rotationalmultipliers  = 0;        % account for m and theta multipliers: 0-> not accounted 1-> accounted
settings.ULS                    = 0;        % 1 = Integration of mobilisable soil resistance in accordance with EA-Phähle approach is plotted, 0 = no plot
settings.lat_cap_10_crit        = 0;        % 

% Critical pile length settings
pile.fixed_lenght_switch        = 0;        % if = 0 the pile length is not fixed (usually manual input), if = 1 pile length is fixed and shows plot with marker for this length (usually database)
pile.criterion                  = 20;       % [%] Maximum critical pile length criterion to be analysed. Only integers! All integer percentages below given number will be analysed

% Multi toe spring upload settings
pile.toe_max                    = 0;        % [m] upper range  
pile.toe_min                    = 10;       % [m] lower range

% Permanent rotation settings
if  plots.permanent_rot_def == 1
	settings.elasticstiffness 	= 'static'; % acceptable input: 'static' or 'cyclic' -- determination of initial stiffness based on static of cyclic curves
	settings.elasticmultipliers = 0;        % input details: 0-> lateral multipliers of 1 accounted in elastic stiffness, 1-> multipliers accounted in elastic stiffness.
end

% Load-displacement settings
if plots.load_deflection == 1 || plots.load_deflection_axial == 1
	settings.max_load_ratio     = 10;       % attempts to apply a load equal to the one imported from Excel/database multiplied by this factor
end
%--------------------------------------------------------------------------
%% Output files
%--------------------------------------------------------------------------
settings.PSI                    = 0;        % Create PSI file? 1 = yes, 0 = no
settings.ANSYS                  = 0;        % Create ANSYS ASAS file? 1 = yes, 0 = no
settings.appendix               = 0;        % Create appendix for report (calculation log)? 1 = yes, 0 = no
data.save_path                  = 'output\';% Saves files in defined working folder, for current 'pwd'
settings.SSI2db                 = 0;        % Save SSI-curves in database? 1 = yes, 0 = no
settings.update_db              = 0;        % update the results in MySQL database (1 = yes, 0 = no)
settings.save_plots             = 0;        % save plots in output folder (1 = yes, 0 =no)
settings.multi_toe_levels       = 0;        % upload to DB toe springs for more than one level (1 = yes, 0 = no)
settings.damping 				= 0;        % uploads needed input for soil damping to excel (1 = yes, 0 = no)
settings.SESAM                  = 0;        % Export Linearised springs for SESAM
settings.PISA_cal_save          = 0;        % Saving moment, shear, etc. for PISA SSI spring calibration
end