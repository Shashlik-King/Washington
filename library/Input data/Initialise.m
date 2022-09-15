function [settings,soil,pile,loads,data,plots, project_info] = Initialise(analysis_setting, revision_setting, soil_pile_setting, ...
    numerical_setting, axial_setting, lateral_setting, output_setting, project_info, database_setting)





data.analysis_id = num2str(analysis_setting.id);

%--------------------------------------------------------------------------
%% Soil and pile input
%--------------------------------------------------------------------------
settings.database               = analysis_setting.input_type{:};       % if 'Database' the database module is activated. If 'Manual' input is taken from manual_data_input.xlsx
settings.PISA_database          = soil_pile_setting.PISA_database;      % Switch for PISA database usage (1 = on, 0 = off)
settings.psf_switch             = soil_pile_setting.psf_switch;
settings.psf_undrained          = soil_pile_setting.psf_undrained;
settings.psf_drained            = soil_pile_setting.psf_drained;
% interim settings
settings.interimloads           = soil_pile_setting.interimloads;
settings.interimgeometry        = soil_pile_setting.interimgeometry;
soil.type_su                    = soil_pile_setting.type_su;              % Soil degradation switch (1 = to use the factors in cyclicDegradation.m , 0 = to not use)

%below only if reading from database is set
loads.type                      = soil_pile_setting.loads_type;          % Loads to be used (ULS or FLS)maybe also GEO
soil.psf                        = soil_pile_setting.soil_psf;            % Partial Safety Factors 0 reads the _geo columns from the load table in the database
soil.type                       = soil_pile_setting.soil_type{:};           % 'BE', 'LB', 'UB'

data.table_springs              = revision_setting.table_springs{:} ;      % 'char', 'LB', 'UB', 'test' Table to update in the database
% WARNING: ONLY FOR SPECIAL USE (only applied if data.revision.global = -1):
data.revision.soil              = revision_setting.rev_soil;               % revision no. of soil parameters to be used (1000 = latest revision)
data.revision.structure         = revision_setting.rev_structure;          % revision no. of structure to be used (1000 = latest revision)
data.revision.loads             = revision_setting.rev_loads;              % revision no. of ULS loads to be used (1000 = latest revision)
data.revision.output            = revision_setting.rev_output;             % revision no. for storing results into the database

settings.interface              = revision_setting.interface;              % FLS: FLS loads are applied, ULS: factored ULS loads are applied,
% GEO: factored/unfactored loads are applied depending on the check
% that are carried out
settings.db_server              = database_setting.db_server{:};   % Databse server % 'DKLYCOPILOD1'
settings.db_user                = database_setting.db_user{:};   % Database user
settings.db_pass                = database_setting.db_pass{:};    %'ituotdewdb';    % Database pass
settings.db_name                = database_setting.db_name{:};  %'ewdb';          % Database name
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% General calculation settings
%--------------------------------------------------------------------------
settings.nelem_factor           = numerical_setting.nelem_factor;        % [m-1] minimum number of pile beam elements per meter 2.2
settings.j_max                  = numerical_setting.j_max;% 100;      % max. no. of iterations in a load step (lateral only)
settings.TOL                    = numerical_setting.TOL; %1e-4;     % relative tolerance (lateral and axial)
settings.n_max                  = numerical_setting.n_max; % 100;        % number of load steps (lateral and axial)
%--------------------------------------------------------------------------
%% Calculation settings (axial)
%--------------------------------------------------------------------------
settings.axial_loading          = 0;                                  % calculate capacity for axial loading? 1 = yes, 0 = no
settings.clay_type              = axial_setting.clay_type{:};         % 'alpha' or 'beta' for way to calculate skin friction in API clay
settings.analysis_type          = axial_setting.analysis_type;        % 0 = normal axial assessment, 1 = settlements
settings.analysis_loading       = axial_setting.analysis_loading;     % {'Tens'} = tension only, {'Comp'} = compression only, {'Tens', 'Comp'} = both
settings.Delta_lambda_bar       = axial_setting.Delta_lambda_bar;     % inital load increment for AL (Arc length) only used in axial calculation. Controls the number of load steps - higher value = fewer load steps, lower value = more load steps
settings.axial_NR_AL            = axial_setting.axial_NR_AL;          % 0 = Newton-Raphson, 1 = Arc length IMPORTANT: A-L solver not fully QA-ed. Use not recomended

plots.res_vs_pilelength         = 0;        % 1 = yes, 0 = no
plots.normal_force_axial        = 0;        % Normal force plot. 1 = yes, 0 = no
plots.UR_axial                  = 0;        % utilization ratio for axial capacity
plots.load_deflection_axial		= 0;		% 1 = yes (should not be combined with other plots; minimum 300 load steps are recommended), 0 = no
plots.deflection_plot_axial     = 0;        % Axial displacement along the pile, 0 = no, 1 = yes
plots.py_curve_UR = 0;
pile.plug_unplug.tens           = 0;        % plug/unplug control in tension - 0 = auto, 1 = plugged, 2 = unplugged
pile.plug_unplug.comp           = 0;        % plug/unplug control in compression - 0 = auto, 1 = plugged, 2 = unplugged
pile.extra_L                    = axial_setting.extra_L;       % [m] extra pile length to be used in axial capacity plot

soil.zpeak                      = axial_setting.soil_zpeak;     % usually = pile.diameter*0.01
soil.zres                       = axial_setting.soil_zres ;     % usually = 2
soil.zpeak_Qz                   = axial_setting.soil_zpeak_Qz;  % usually = pile.diameter*0.1
soil.tres                       = axial_setting.soil_tres;      % Usually = 0.7 - 0.9


%--------------------------------------------------------------------------
%% Calculation settings (torsional)
%--------------------------------------------------------------------------
settings.torsion                = 0;        % calculate torsional stiffness? 1 = yes, 0 = no. settings.axial_loading must also be set to 1
%--------------------------------------------------------------------------
%% Calculation settings (lateral)
%--------------------------------------------------------------------------
loads.static_cyclic             = lateral_setting.loads_type{:}; % 'cyclic' or 'static' in accordance with DNV OS-J101
loads.n_cycles                  = lateral_setting.loads_n_cycles; % 100;      % Number of cycles, relevant for Stiff clay w/o water only!
loads.A                         = lateral_setting.loads_A{:}; % 'API';    % use API or TUHH (D�hrkop) approach for determination of A, relevant for API/Kirsch sand



if strcmp(soil.type,'UB')
    loads.static_cyclic= 'static';
end


plots.pilehead_vs_length        = 0;        % 1 = yes, 0 = no
plots.deflection_plot           = 0;        % 1 = yes, 0 = no
plots.utilization_ratio         = 0;        % 1 = yes, 0 = no, settings.toe_shear = 0 if only p-y UR is of interest
plots.deflection_bundle         = 0;        % 1 = yes, 0 = no
plots.toe_shear_graph           = 0;        % 1 = yes, 0 = no
plots.permanent_rot_def         = 0;        % 1 = yes (minimum 10 load steps is recommended), 0 = no -- cannot be used together with other analyses
plots.load_deflection			= 0;		% 1 = yes (should not be combined with other plots; minimum 50 load steps us recommended), 0 = no
plots.load_deflection_type      = 0;		% 0 = DNV-Gl / BSH , 1 = regular load-displacement
plots.moment_distribution		= 0;		% 1 = yes, 0 = no

settings.lateral_loading        = 0;        % calculate capacity for lateral loading? 1 = yes, 0 = no
settings.beam_theory            = lateral_setting.beam_theory;        % 1 = Timoshenko, 0 = Euler-Bernoulli
settings.toe_shear              = lateral_setting.toe_shear;        % include base shear and moment? 1 = yes, 0 = no
settings.mteta	                = lateral_setting.mteta;        % include uniformly distributed moment? 1 = yes, 0 = no
settings.Georgiadis				= lateral_setting.Georgiadis;		% apply Georgiadis approach? 1 = yes, 0 = no
settings.lateralmultipliers 	= lateral_setting.lateralmultipliers; 		% account for p and y multipliers: 0-> not accounted 1-> accounted
settings.rotationalmultipliers  = lateral_setting.rotationalmultipliers;        % account for m and theta multipliers: 0-> not accounted 1-> accounted
settings.ULS                    = lateral_setting.ULS;        % 1 = Integration of mobilisable soil resistance in accordance with EA-Ph�hle approach is plotted, 0 = no plot
settings.lat_cap_10_crit        = lateral_setting.lat_cap_10_crit;        %

% Critical pile length settings
pile.fixed_lenght_switch        = lateral_setting.fixed_lenght_switch;  % if = 0 the pile length is not fixed (usually manual input), if = 1 pile length is fixed and shows plot with marker for this length (usually database)
pile.criterion                  = lateral_setting.pile_criterion;       % [%] Maximum critical pile length criterion to be analysed. Only integers! All integer percentages below given number will be analysed

% Multi toe spring upload settings
pile.toe_max                    = lateral_setting.pile_toe_max;        % [m] upper range
pile.toe_min                    = lateral_setting.pile_toe_min;        % [m] lower range

%--------------------------------------------------------------------------
%% Output files
%--------------------------------------------------------------------------
settings.PSI                    = output_setting.PSI;        % Create PSI file? 1 = yes, 0 = no
settings.ANSYS                  = output_setting.ANSYS;        % Create ANSYS ASAS file? 1 = yes, 0 = no
settings.appendix               = output_setting.appendix;        % Create appendix for report (calculation log)? 1 = yes, 0 = no
data.save_path                  = output_setting.save_path{:}; %'output/';% Saves files in defined working folder, for current 'pwd'
settings.SSI2db                 = 0;                                % Save SSI-curves in database? 1 = yes, 0 = no
settings.update_db              = output_setting.update_db;        % update the results in MySQL database (1 = yes, 0 = no)
settings.save_plots             = output_setting.save_plots;        % save plots in output folder (1 = yes, 0 =no)
settings.multi_toe_levels       = output_setting.multi_toe_levels;        % upload to DB toe springs for more than one level (1 = yes, 0 = no)
settings.damping 				= output_setting.damping;        % uploads needed input for soil damping to excel (1 = yes, 0 = no)
settings.SESAM                  = output_setting.SESAM;        % Export Linearised springs for SESAM
settings.PISA_cal_save          = output_setting.PISA_cal_save;        % Saving moment, shear, etc. for PISA SSI spring calibration


%--------------------------------------------------------------------------
%% Analysis Types
%--------------------------------------------------------------------------

if strcmpi(analysis_setting.type,'lateral_load_deflection')
    
    settings.lateral_loading = 1;
    plots.load_deflection = 1;
    disp('Computing lateral load aganist deflection')
    
elseif strcmpi(analysis_setting.type,'lateral_deflection')
    
    settings.lateral_loading = 1;
    plots.deflection_plot = 1;
    disp('Computing lateral deflection along pile')

    
elseif strcmpi(analysis_setting.type,'lateral_load_rotation')
    
    settings.lateral_loading = 1 ;
    plots.permanent_rot_def  = 1 ;        % 1 = yes (minimum 10 load steps is recommended), 0 = no -- cannot be used together with other analyses
    disp('Computing lateral load against pile head rotation')

elseif strcmpi(analysis_setting.type,'export_springs')
    
    settings.lateral_loading = 1;
    settings.axial_loading = 1;
    plots.deflection_plot = 1;
    settings.SSI2db = 1;
    
    disp('Computing springs')


else
    
    error([analysis_setting.type{:},' was not defined !!!'])
    
end

%%

% Permanent rotation settings
if  plots.permanent_rot_def == 1
    settings.elasticstiffness 	= 'static'; % acceptable input: 'static' or 'cyclic' -- determination of initial stiffness based on static of cyclic curves
    settings.elasticmultipliers = 0;        % input details: 0-> lateral multipliers of 1 accounted in elastic stiffness, 1-> multipliers accounted in elastic stiffness.
end

% Load-displacement settings
if plots.load_deflection == 1 || plots.load_deflection_axial == 1
    settings.max_load_ratio     = numerical_setting.max_load_ratio;       % attempts to apply a load equal to the one imported from Excel/database multiplied by this factor
end










