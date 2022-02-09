%% CALCULATION ROUTINE FOR AXIALLY AND LATERALLY LOADED PILES
% RUN FILE FOR CMAT CALCULATION TOOL
% CAN PERFORM AXIAL AND LATERAL BEARING CAPACITY CALCULATIONS
clear
close all
clc
%--------------------------------------------------------------------------
% Units: kN, m, s, kPa
%--------------------------------------------------------------------------
% CHANGE LOG
% YYYY.MM.DD    USER    Task
% 2010.11.04    MMOL    Programming
% 2011.0X.XX    JALY    Streamlining code
% 2013.07.XX    MUOE    Cleaning and streamlining
%                       Adding t-z curves and SACS PSI-input module
%                       Adding weak rock p-y module
% 2013.10.16    MUOE    Adding print-out of
%                       documentation/assumptions during calculation
% 2014.07.15    MMOL    Cleaning and streamlining
% 2018.04.04	DATY	On/Off switch for Georgiadis approach
% 2019.10.25    ASSV    Changed the folder structure, creating a library,
%                       cleaning, add database input/output modules from SNA
%                       and corresponding commands in the main, removed the
%                       p-y saving on excel
%2019.11.20		ASSV	Corrected the pile length used in the PISA base elements functions
%                       Winkler solver modified to enhance convergence with
%                       base moment springs
%2019.11.22     ASSV    Database interaction fully working, all the inputs 
%                       can be read from mysql and PISA springs saved
%2020.04.06     FKMV    Settlement calculations using N-R and Arc-length solver 
%2020.08.17     FKMV    Added feature to use PISA database for parameters, generalisation of the code from project specific to M.V.
%--------------------------------------------------------------------------
addpath (genpath('library'));                           %make sure all the functions are available

% ID ={'L0GW01_LB'}; % API

Water_Depth='DEEP_PISA' ;   %'DEEP' or   'Shallow'
ID ={'EW1-69'}; % API Loads received 2020-11-06
loadtype='ULS-DEEP';  % or ULS-DEEP
[loads_total, names, DLC] = extractloads_fkmv(Water_Depth);
idx=find(and(contains(names(:,2),ID{1,1}),contains(names(:,3),loadtype)));
ID_Soil_in_db=names{idx,1};
loads.type                      = 'FLS';          % Loads to be used (ULS or FLS)maybe also GEO
soil.psf                        = 0;              % Partial Safety Factors 0 reads the _geo columns from the load table in the database
settings.psf_switch              =0;
settings.psf_undrained           =1.5;

data.Strength_B_condition       =1 ;   % only a condition valid for empire wind. in this case the phi and C are LB value of data base
ID_Soil_in_db='Z1_L08_BE_PISA';

for loc = 1:length(ID)
    
    
    
%% Project data
%--------------------------------------------------------------------------
data.project                    = 'EW OWF';   % 'Project name'
data.A_number                   = 'A205414';   % 'Project number'
data.location                   = ID_Soil_in_db;   % 'Project location'
data.db_location                = ID;
data.db_location_geo            = ID_Soil_in_db; %'Database location to read data
data.prepared_by                = 'PNGI';   %  Initials of preparer
data.zoneID                     =ID_Soil_in_db;
data.id                         =ID; % Id to store into the database   
data.revision.global            = 1;   % global revision no. for selected location to be used,
                                        %1000 detects corresponding soil, structure and load revision no. automatically
                                        %-1 reads the settings below 
%--------------------------------------------------------------------------
%% Soil and pile input
%--------------------------------------------------------------------------
settings.database               = 'Manual';       % if 'Database' the database module is activated. If 'Manual' input is taken from manual_data_input.xlsx
settings.PISA_database          = 1;              % Switch for PISA database usage (1 = on, 0 = off)

% interim settings
settings.interimloads           = 0;
settings.interimgeometry        = 0;
soil.type_su                    = 0;              % Soil degradation switch (1 = to use the factors in cyclicDegradation.m , 0 = to not use)

%below only if reading from database is set    
soil.type                       = 'BE';           % 'BE', 'LB', 'UB'

data.table_springs              = 'char' ;        % 'char', 'LB', 'UB', 'test' Table to update in the database
% WARNING: ONLY FOR SPECIAL USE (only applied if data.revision.global = -1):
data.revision.soil              = 0;              % revision no. of soil parameters to be used (1000 = latest revision)
data.revision.structure         = 0;              % revision no. of structure to be used (1000 = latest revision)
data.revision.loads             = 0;              % revision no. of ULS loads to be used (1000 = latest revision)
data.revision.output            = 0;             % revision no. for storing results into the database  

settings.interface              = 'FAC';          % FLS: FLS loads are applied, ULS: factored ULS loads are applied, 
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
settings.n_max                  = 200;        % number of load steps (lateral and axial)
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
loads.A                         = 'P_NGI';    % use API or TUHH (Dührkop) approach for determination of A, relevant for API/Kirsch sand



if strcmp(soil.type,'UB')
loads.static_cyclic= 'static';
end 


plots.pilehead_vs_length        = 0;        % 1 = yes, 0 = no
plots.deflection_plot           = 0;        % 1 = yes, 0 = no
plots.utilization_ratio         = 0;        % 1 = yes, 0 = no, settings.toe_shear = 0 if only p-y UR is of interest
plots.deflection_bundle         = 0;        % 1 = yes, 0 = no
plots.toe_shear_graph           = 0;        % 1 = yes, 0 = no
plots.permanent_rot_def         = 1;        % 1 = yes (minimum 10 load steps is recommended), 0 = no -- cannot be used together with other analyses
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
	settings.max_load_ratio     = 5;       % attempts to apply a load equal to the one imported from Excel/database multiplied by this factor
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
settings.multi_toe_levels       = 1;        % upload to DB toe springs for more than one level (1 = yes, 0 = no)
settings.damping 				= 0;        % uploads needed input for soil damping to excel (1 = yes, 0 = no)
settings.SESAM                  = 0;        % Export Linearised springs for SESAM
settings.PISA_cal_save          = 0;        % Saving moment, shear, etc. for PISA SSI spring calibration
%--------------------------------------------------------------------------
%% Collecting input data from source
%--------------------------------------------------------------------------
disp(['Position ',data.location,' at ',data.project])
[SF reduction] = factors(loads); %#ok<*NCOMMA>
if strcmp(settings.database,'Database')
    [pile loads settings plots scour soil] = database_input(pile,soil,data,...
        loads,settings,plots); % load soil and pile information from database (SNA module)
    soil.q_ur                = soil.q_ur_py; % [kPa]
    soil.delta_q_ur = zeros(length(soil.q_ur_py));
    if settings.interimloads 
    [loads] = interimLoads (plots,ID{loc},soil,loads,settings);
    disp('loads inserted manually in run_COSPIN - correct when loads are available in database')
    end
    if settings.interimgeometry 
    [pile] = interimGeometry (plots,ID{loc},soil,loads,settings,pile,scour);
    disp('MP geometry inserted manually in run_COSPIN - correct when MP geometry is available in database')
    end
elseif strcmp(settings.database,'Manual')
    [scour soil pile loads settings] = manual_data_input_excel(pile,data,soil,...
        loads,settings); % load soil and pile information from Exce-file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%over writing the loads from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%excel

if strcmp(loads.type,'ULS')
    loads.H=DLC.FXY.ULS(idx,1);   
    loads.M = -DLC.MXY.ULS(idx,1); % [kNm] overturning moment at pile head
    loads.Vc = DLC.Fz.ULS(idx,1); % [kN] compressive axial force at pile head
    loads.Vt = 0; % [kN] tensile axial force at pile head
    loads.Mz = DLC.Mz.ULS(idx,1); % [kN] torsional moment
elseif strcmp(loads.type,'FLS')
    loads.H=DLC.FXY.FLS(idx,1);   
    loads.M = -DLC.MXY.FLS(idx,1); % [kNm] overturning moment at pile head
    loads.Vc = DLC.Fz.FLS(idx,1); % [kN] compressive axial force at pile head
    loads.Vt = 0; % [kN] tensile axial force at pile head
    loads.Mz = DLC.Mz.FLS(idx,1); % [kN] torsional moment
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%over writing the pile geometry ###############################


pile.length_start=loads_total(idx,11);

if plots.pilehead_vs_length        == 1
pile.length_end=loads_total(idx,11)+30;
else    
pile.length_end=loads_total(idx,11) ;   
end 
pile.length  = (pile.length_start:pile.length_inc:pile.length_end); % [m] creates a vector for rotation vs. embedment if needed

pile.diameter=loads_total(idx,10);
pile.cross_section.diameter(:)=loads_total(idx,10);
pile.cross_section.thickness(:)  = loads_total(idx,12);
pile.cross_section.endarea  = pi*((pile.cross_section.diameter(end)/2)^2-(pile.cross_section.diameter(end)/2-pile.cross_section.thickness(end))^2); % [m^2] end bearing area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read PISA database
if settings.PISA_database
    [database.num,database.txt,database.raw] = xlsread('Data_Base.xlsx','Clay','A1:BP1000'); % reading of Database table
    database.txt = database.txt(4:end,:); % removing top 3 header lines of database
    database.num(isnan(database.num))=0; % NaN in read table due to a `-` or a different characted than a number to be converted to 0
else
    database = 0; % initialisation of the variable to avoind errors when not using PISA
end

output = []; %initialisation of output

% allow for setup of different calculation method for loading/unloading path
if  plots.permanent_rot_def == 1
	cases = [1,1]; % if permanent rotation calculations, we loop twice over same pile length
    % add change in load steps as well??
else
	cases = 1:length(pile.length); % Looping over all pile lengths...
                                   %(no reloading path to be determined for this kind of calculation)
end

iter_global = 0;
%%
 for i = cases
 iter_global = iter_global+1;

if plots.permanent_rot_def == 1 && iter_global == 1
    if settings.elasticmultipliers > settings.lateralmultipliers 
        error(['Permanent rotation calculation: Multipliers below 1 only used',...
        'for unloading stiffness - check variable:',...
        ' settings.lateralmultipliers.']) % safeguard if multipliers are 'on' 
                                          %for unloading stiffness and 'off' for loading stiffness
    end 
    fprintf('\nCalculation to evaluate initial stiffness for unloading part of permanent rotation check.\n\n');
	store1 							= loads.static_cyclic;
	store2 							= settings.lateralmultipliers;
	loads.static_cyclic 			= settings.elasticstiffness;
	settings.lateralmultipliers 	= settings.elasticmultipliers;
elseif plots.permanent_rot_def == 1 && iter_global == 2
    fprintf('\nCalculation to evaluate loading path of permanent rotation at pile head check.\n\n'); 
	loads.static_cyclic 			= store1;
	settings.lateralmultipliers 	= store2;
end

if loads.n_cycles==10000 && plots.permanent_rot_def == 1 && iter_global == 1
    loads.H=loads.H*0.3;
    loads.M=loads.M*0.3;
end

[element,node,pile] = elem_node_data(pile,soil,scour,settings,data,plots,i,database); % assign properties to elements and nodes
if settings.rotationalmultipliers
    [element,soil] = interimRotationalMultipliers (element,soil); % apply multipliers on rotational prings (PISA  only) usually not used
end
if soil.type_su
    [element] = cyclicDegradation(element); % apply cyclic degradation for different soil layers only to cu
end
%--------------------------------------------------------------------------
%% Axial calculation / t-z curves / Q-z curves
%--------------------------------------------------------------------------
if settings.axial_loading    
    disp('Axial calculation initiated')
    [fs qp A G] = skin_tip(element,reduction,pile,settings,loads,plots); % Skin friction and tip resistance calculation
    [plug_unplug G output pile] = axial_uls(data,A,G,fs,qp,element,pile,SF,reduction,ID,output,i,plots,loads); % Calculation of total pile resistance
    if settings.analysis_type
        for ii=1:length(settings.analysis_loading)
            if settings.axial_NR_AL == 0 
                [element Ndof Coord Mmax Es u ustep Ex Ey output] = winkler_axial(element,node,pile,loads,reduction,plug_unplug,G,settings,output,plots,SF,A,fs,qp,ii);
            elseif settings.axial_NR_AL == 1
                [element Ndof Coord Mmax Es u ustep Ex Ey output] = winkler_axial_arc_Chrisfield(element,node,pile,loads,reduction,plug_unplug,G,settings,output,plots,SF,A,fs,qp,ii); % Winkler analysis for vertical loads 
            end
            if plots.UR_axial == 1
                [output] = UR_axial(element,settings,pile,loads,output,reduction,plug_unplug,A,fs,ii); % Calculates utilization ratio for axial load 
            end
        end
    [output] = plot_functions_axial(element,pile,output,Coord,settings,plots,loads,Ex,Ey,Es,i,node,soil,pile.L,G,ii,data);
    end
    [t z] = t_z(fs,pile,element,reduction,plug_unplug,A,G); % calculation of t-z curves for soil layers according to API
    [Q zQ] = Q_z(qp,pile,element); % calculation of Q-z curves for soil layers according to API
    disp('Axial calculation completed')
%--------------------------------------------------------------------------
%% Torsional stiffness
%--------------------------------------------------------------------------
if settings.torsion
    data.torsion = torsion_stiffness(pile,element,t);
end
end
%--------------------------------------------------------------------------
%% Lateral calculation / p-y curves / lateral UR
%--------------------------------------------------------------------------
if settings.lateral_loading
    disp('Lateral calculation initiated')
    [element Ndof Coord Mmax Es u ustep output] = winkler(element,node,pile,loads,settings,output,plots,data); % calcutation routine for lateral loading
    for j = 1:element.nelem+1
        output.hor_defl(j,i) = u(j*3-2,1); % save horisontal pile deflection at final loading increment for each pile length
        output.rot(j,i) = u(j*3,1)*180/pi;
    end

    output.pilehead_rotation(1,i) = -u(3,1)/pi*180; % save pile head rotation
    [p y y_tot output toe_plot_u] = p_y(node,pile,element,loads,output,settings); % calculation of p-y curves for soil layers
    if settings.mteta
        [m teta output toe_plot_teta] = m_teta(node,pile,element,loads,output,settings);
    end
    if settings.toe_shear
        [p_toe,y_toe] = p_y_toe(node,pile,element,loads,output,settings);
        [m_toe,teta_toe] = m_teta_toe(node,pile,element,loads,output,settings);
    end
    if plots.load_deflection == 0
		[output] = UR_v2(element,settings,pile,loads,data,plots,ustep,y_tot,output,node);
    end
    disp('Lateral calculation completed')
    output.Coord = Coord;
    output.toe_plot_u = toe_plot_u;
end
%--------------------------------------------------------------------------
%% Plots and other output
%--------------------------------------------------------------------------
if plots.permanent_rot_def == 1 && iter_global == 1
	plots.node_rot_def= 1;                              % number of node to plot permanent rotations for, 1 = node at pile head
	F = linspace(0,1,settings.n_max+1); % this is valid because the load is applied in equally sized steps - the magnitude of the load doesn't matter, only the fact that it is applied in equally sized steps
    output.elasticstiff = (F(2)-F(1))/(output.deflections(3*plots.node_rot_def,2)-output.deflections(3*plots.node_rot_def,1)); % the unloading/reloading stiffness is calculated as the initial stiffness
else
	[output] = plot_functions(element,pile,node,soil,plots,output,settings,i, loads, data,SF);
end


if settings.PSI
    disp('Printing to SACS PSI-file')
    PSI(pile,data,element,scour,settings,t,z,p,y,Q,zQ,plug_unplug,i)
    disp('Finished printing to SACS PSI-file')
end
if settings.ANSYS
    disp('Printing to ANSYS ASAS-file')
	[p_TDA t_TDA level_TDA]   = TDA(node,p,t,data,scour,pile,y,z,settings);
    %ANSYS_ASAS(p,t,node,y,z,scour,pile,data,i)
	ANSYS_ASAS_TDA(p_TDA,t_TDA,level_TDA,y,z,scour,pile,data,i,Q,zQ,plug_unplug,element);
    disp('Finished printing to ANSYS ASAS-file')
end    

if settings.SSI2db
    if settings.lateral_loading
        database_write_springs(settings,element,loads,data,p,y,...
        strcat('soil_py_curves_',data.table_springs),'py'); % write py-curves into database
    end
    if settings.axial_loading
        database_write_springs(settings,element,loads,data,t,z,...
        strcat('soil_tz_curves_',data.table_springs),'tz'); % write tz-curves into database
    end
    if settings.mteta    
        database_write_Mtsprings(settings,element,loads,data,m,teta,...
            strcat('soil_Mt_curves_',data.table_springs),'Mt'); % write Mt-curves into database
    end
    if settings.toe_shear
        if settings.multi_toe_levels
            database_write_toesprings_multi(settings,element,loads,data,pile,soil,plots,output,i,p_toe,y_toe,m_toe,teta_toe,database) % writes toe springs for same revision for different depths
        else
            database_write_toesprings(settings,element,loads,data,p_toe,y_toe,...
                strcat('soil_py_toe_curves_',data.table_springs),'py'); % write py-curves into database
            database_write_toesprings(settings,element,loads,data,m_toe,teta_toe,...
                strcat('soil_Mt_toe_curves_',data.table_springs),'Mt'); % write Mt-curves into database
        end
    end
end
			 
if settings.appendix
    disp('Writing appendix to Word-file')
    
    disp('Finished writing appendix to Word-file')
end
    disp('--------------------------------------')
end
% auto_documentation(data,settings, soil, plots, loads);
%F = linspace(0,loads.H,settings.n_max+1)/1000; % This is for Hv plots in PISA (temporal)
if settings.damping
	soil_damping_input_write(data,Coord,output,loads)
end

if settings.PISA_cal_save
    save_results_PISA_calibration(data,element,output,plots,Es) % FKMV addition
%     PISA_calibration_plots(data, plots)
end

if settings.SESAM
    [Stiff_node_py,Stiff_node_mtheta] = SESAM(node,scour,pile,data,element,loads,u,ustep,settings,output); % FKMV addition
	QA_springs(data,loads,soil,scour,settings,pile)
end

end