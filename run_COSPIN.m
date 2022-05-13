
%% CALCULATION ROUTINE FOR AXIALLY AND LATERALLY LOADED PILES
% RUN FILE FOR CMAT CALCULATION TOOL
% CAN PERFORM AXIAL AND LATERAL BEARING CAPACITY CALCULATIONS
clear
close all
clc

if (~isdeployed)
    addpath (genpath('library'));
    addpath (genpath('input'));   %make sure all the functions are available
end
%--------------------------------------------------------------------------
% Units: kN, m, s, kPa
%--------------------------------------------------------------------------


%% Import data blocks from the cospin_input

[~,~,raw0]=xlsread('COSPIN_input.xlsx','PROJ');

[analysis_settings]=read_data_block('Analysis settings',raw0);
is_run = analysis_settings.is_run ==1; % filter the desired runs
if sum(is_run)==0
    error('No analysis was chosen to run/ There are empty cells for the chosen analysis !!! Check PROJ SHEET');
end

analysis_settings = analysis_settings(is_run,:);
[soil_pile_settings]=read_data_block('Soil and pile settings',raw0);
[revision_settings]=read_data_block('Revisions settings',raw0);
[database_settings]=read_data_block('Database settings',raw0);
[numerical_settings]=read_data_block('Numerical settings',raw0);
[axial_settings]=read_data_block('Axial settings',raw0);
[lateral_settings]=read_data_block('Lateral settings',raw0);
[output_settings]=read_data_block('Output settings',raw0);
[project_infos]=read_data_block('Project infos',raw0);
[pile_materials]=read_data_block('Pile material settings',raw0);

% Read LOCATIONS SHEET
[~,~,Locations]=xlsread('COSPIN_input.xlsx','LOCATIONS');
idx1=8; idx2=idx1+10*4-1;
geometries = Locations(1,idx1:idx2);  cond1 = ~ismissing(string(geometries));
geometries = geometries(cond1);

loc_info=Locations(4:end,:); % remove the headers
is_run=[loc_info{:,1}]'==1;
if sum(is_run)==0
    error('No location was chosen to run !!! Check LOCATIONS SHEET');
end
loc_info=loc_info(is_run,:); % filter the run locations
geometries_table=loc_info(:,idx1:idx2);
ID          = loc_info(:,2); %ID = {'benchmark'};
database_id = loc_info(:,3);
water_depth = loc_info(:,4);
scour_local = loc_info(:,5) ;
scour_ORD   = loc_info(:,6);


% Read LOADS SHEET
[~,~,Loads]=xlsread('COSPIN_input.xlsx','LOADS');
idx1=4; idx2=idx1+10*5-1;

load_types = Loads(1,idx1:idx2);  cond1 = ~ismissing(string(load_types));
load_types = load_types(cond1);

Loads_info=Loads(4:end,:); % remove the header
is_run=[Loads_info{:,1}]'==1;
Loads_table=Loads_info(is_run,idx1:idx2); % filter the run locations



for loc = 1:length(ID) % Loop through the locations
    
    disp(['**Running location : ', ID{loc},'**'])
    
    for jj=1:size(analysis_settings,1) % Loop through the analyses
        
        analysis_setting  = analysis_settings(jj,:);
        disp(['**Conducting analysis_id = ' num2str(analysis_setting.id),' - ',analysis_setting.label{:},'**']);

        revision_setting  = revision_settings(analysis_setting.revision_setting_id,:);
        soil_pile_setting = soil_pile_settings(analysis_setting.soil_pile_setting_id,:);
        numerical_setting = numerical_settings(analysis_setting.numerical_setting_id,:);
        axial_setting     = axial_settings(analysis_setting.axial_setting_id,:);
        lateral_setting   = lateral_settings(analysis_setting.lateral_setting_id,:);
        output_setting    = output_settings(analysis_setting.lateral_setting_id,:);
        project_info      = project_infos(analysis_setting.project_info_id,:);
        database_setting  = database_settings(analysis_setting.database_setting_id,:);
        pile_material     = pile_materials(analysis_setting.pile_material_id,:);
        
        [settings, soil, pile, loads, data, plots] = Initialise(analysis_setting, revision_setting, soil_pile_setting, ...
            numerical_setting, axial_setting, lateral_setting, output_setting, project_info, database_setting);
        
        %% Project data
        %--------------------------------------------------------------------------
        data.project                    = project_info.name{:};   % 'Project name'
        data.A_number                   = project_info.ATRnumber{:};  % 'Project number'
        data.location                   = ID{loc};   % 'Project location'
        data.db_location                = [project_info.STR_db_Loc_prefix{:},'_', ID{loc}];
        data.db_location_geo            = [project_info.GEO_db_Loc_prefix{:},'_',ID{loc}]; %'Database location to read data
        data.prepared_by                = project_info.prepared_by;   %  Initials of preparer
        data.zoneID                     = ID{loc};
        data.id                         = [project_info.STR_db_Loc_prefix{:},'',database_id{loc}]; % Id to store into the database
        data.revision.global            = analysis_setting.db_export_revision;   % global revision no. for selected location to be used,
        %1000 detects corresponding soil, structure and load revision no. automatically
        %-1 reads the settings below
        %--------------------------------------------------------------------------
        
        if length(data.id)>7
            error('Total sum of charcters for *STR_db_Loc_prefix* and *database_id* should be less than 7');
        end
        
        if length(data.revision.global )>2
            error('Total charcters for *db_revision_global* should be less than 2*');
        end
        
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
        
        %% Overwrite loads, pile material and geometry settings from manual_data_input_excel
        
        % overwrite_manual_data_input_loads
        cond1 = strcmpi(load_types, analysis_setting.load_type); idx = 1:length(load_types); idx0 = idx(cond1);
        if isempty(idx0)
           error(['In analysis settings, load_type = ' analysis_setting.load_type{:}, ' does not exist in LOADS sheet']) 
        end
        
        disp(['Loads were overwrriten. ' load_types{idx0} ' was applied!'])
        Range0 = (1+(idx0-1)*5):(5+(idx0-1)*5);
        Loads_table_loc = Loads_table(loc,Range0);
        loads.H  = Loads_table_loc{1};
        loads.M  = Loads_table_loc{2};
        loads.Vc = Loads_table_loc{3};
        loads.Vt = Loads_table_loc{4};
        loads.Mz = Loads_table_loc{5};
        
        % overwrite_manual_data_input_pile_material
        pile.density                    = pile_material.pile_density; % [kN/m^3] density of pile material (steel)
        pile.sigma_y                    = pile_material.sigma_y;      % [kPa] characteristic steel yield strength
        pile.E                          = pile_material.E;            % [kPa] Young's modulus of steel
        pile.G                          = pile_material.G;            % [kPa] Shear modulus
        pile.ksf                        = pile_material.ksf;          % [-] Shear correction factor for Timoshenko beam theory
        
        % overwrite_manual_data_input_scour
        scour.local = scour_local{loc};
        scour.ORD   = scour_ORD{loc}; % [m] overburden reduction depth
        scour.water_depth = water_depth{loc}; % [m] water depth
        disp('Scour info were overwritten according to LOCATIONS sheet')

        
        % overwrite_manual_data_input_geometry
        if ~strcmpi(analysis_setting.pile_geometry,'variable')
            disp('Pile Length, diameter and thickness were overwritten according to LOCATIONS sheet. Pile has uniform cross section.')

            cond1 = strcmpi(geometries, analysis_setting.pile_geometry); idx = 1:length(geometries); idx0 = idx(cond1);
            if isempty(idx0)
            error(['In analysis settings, pile_geometry = ' analysis_setting.pile_geometry{:}, ' does not exist in LOCATIONS sheet']) 
            end
            Range0 = (1+(idx0-1)*4):(4+(idx0-1)*4);
            geometries_table_loc = geometries_table(loc,Range0);
            
            Length      =  geometries_table_loc{1};
            Diameter    =  geometries_table_loc{2} ;
            Thickness   =  geometries_table_loc{3} ;
            pile_type   =  geometries_table_loc{4} ;
           
            pile.head                       = 0;         % [m]
            pile.length_start               = Length;    % [m] embedment length, may be a vector
            pile.length_end                 = Length;    % [m]
            pile.length_inc                 = 0.5;       % [m]
            pile.diameter                   = Diameter;  % [m] outer diameter of the pile
            pile.cross_section.toplevel     = 0;         % [m VREF]
            pile.cross_section.thickness    = Thickness; % [m]
            pile.length                     = (pile.length_start:pile.length_inc:pile.length_end);  % [m] creates a vector for rotation vs. embedment if needed
            pile.stick_up                   = abs(pile.head-0)+scour.local;                         % [m] calculate pile stick-up
            
            settings.pile_type              = geometries_table_loc{4}; % [text]
            if strcmp(settings.pile_type,'open')
                pile.cross_section.endarea  = pi*((pile.diameter/2)^2-(pile.diameter/2-pile.cross_section.thickness(end))^2); % [m^2] end bearing area
            elseif strcmp(settings.pile_type,'closed')
                pile.cross_section.endarea  = pi*(pile.diameter/2)^2;
                msgbox('Pile is closed-ended. For correct axial capacity calculation please set reducion.skin_inner = 0.0 in the file reduction_factors.m as this will exclude internal shaft friction')
            end
        end
        
        %%

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
        
        if plots.deflection_plot
            %     fid= fopen('output/deflection.txt','w');
            %     fprintf(fid,'%.5f',output.hor_defl);
            %     fclose(fid);
            writematrix(output.hor_defl,'output/deflection.txt','Delimiter','tab')
        elseif plots.load_deflection
            %     fid= fopen('output/disp.txt','w');
            %     fprintf(fid,'%.5f',output.def_calibration);
            %     fclose(fid);
            %     fid= fopen('output/load.txt','w');
            %     fprintf(fid,'%.5f',output.force_calibration);
            %     fclose(fid);
            writematrix(output.def_calibration,'output/disp.txt','Delimiter','tab')
            writematrix(output.force_calibration,'output/load.txt','Delimiter','tab')
        end
        
    end
    
end