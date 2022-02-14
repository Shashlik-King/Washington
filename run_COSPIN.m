%% CALCULATION ROUTINE FOR AXIALLY AND LATERALLY LOADED PILES
% RUN FILE FOR CMAT CALCULATION TOOL
% CAN PERFORM AXIAL AND LATERAL BEARING CAPACITY CALCULATIONS
clear
close all
clc
%--------------------------------------------------------------------------
% Units: kN, m, s, kPa
%--------------------------------------------------------------------------
addpath (genpath('library'));   
addpath (genpath('input'));   %make sure all the functions are available
[ID,settings,soil,pile,loads,data,plots] = Initialise();
for loc = 1:length(ID)
%% Project data
%--------------------------------------------------------------------------
data.project                    = 'XXXX';   % 'Project name'
data.A_number                   = 'AXXXXXXX';   % 'Project number'
data.location                   = ID{loc};   % 'Project location'
data.db_location                = ID{loc};
data.db_location_geo            = ID{loc}; %'Database location to read data
data.prepared_by                = 'XXXX';   %  Initials of preparer
data.zoneID                     = ID{loc};
data.id                         = ID{loc}; % Id to store into the database   
data.revision.global            = 1;   % global revision no. for selected location to be used,
                                        %1000 detects corresponding soil, structure and load revision no. automatically
                                        %-1 reads the settings below 
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