function [scour soil pile loads data settings] = database_input(pile,soil,data,loads,settings,plots) %#ok<*NCOMMA>
%% WIKINGER DATABASE MODULE
%--------------------------------------------------------------------------
% CHANGE LOG
% 2010.11.04    MMOL    Programming
% 2011.0XXXX    JALY    Streamlining code
% 2011.06.23    MMOL    Re-programming to match Dantysk OWF
% 2013.09.30    MUOE    Re-programming to match Wikinger OWF
% 2019.10.25    ASSV    Re-programming to take settings from the main and
% removed project specific lines, i.e. rotation plots, idea i that in the
% main can be point to manual input, the assignment of model_py is trough a
% legend _code, it can be used also for Arcadis but it needs to be agreed

% Access MySQL-database
mysql('open',settings.db_server,settings.db_user,settings.db_pass); % ('open','server','username','password')
mysql(['use ',settings.db_name]); % name of database

%% Database unique id's
%--------------------------------------------------------------------------
location_geo = ['''',data.db_location_geo,'''']; % name of location geo
location = ['''',data.db_location,'''']; % name of location
rev_global = data.revision.global; % global revision no. for specified location to be used
%%
%--------------------------------------------------------------------------
% Revision no. for soil, load and structure (table: configurations)
%--------------------------------------------------------------------------

if rev_global == -1 % then manual definition for selection of soil, structure and load revision
    rev_soil        = data.revision.soil; % revision no. of soil parameters to be used
    rev_structure   = data.revision.structure; % revision no. of structure to be used
    rev_load        = data.revision.loads; % revision no. of loads to be used
    if rev_soil ~= 1000 % check, if specified revision numbers are available for specified location
        [rev] = mysql(['select soil_revision from soil_data_detail where id=',location_geo]);
        if ismember(rev_soil,str2double(rev)) == 0 % if specified soil revision is not available then output an error message
            error('Specified soil revision number doesn''t exist for this location');
        end
    else % look for latest soil revision
        [rev] = mysql(['select soil_revision from soil_data_detail where id=',location_geo]);
        rev_soil = max(str2double(rev));
    end
    if rev_structure ~= 1000 % check, if specified revision numbers are available for specified location
        [rev] = mysql(['select structural_revision from mp_cans where id=',location]);
         
        if ismember(rev_structure,str2double(rev)) == 0 % if specified structure revision is not available then output an error message
            error('Specified structural revision number doesn''t exist for this location');
        end
    else
        [rev] = mysql(['select structural_revision from mp_cans where id=',location]);
        rev_structure = max(str2double(rev));
    end
    if rev_load ~= 1000 % check, if specified revision numbers are available for specified location
        [rev] = mysql(['select load_rev from loads_',loads.type,' where id=',location]);
        load_condition = 1;
        if ismember(rev_load,str2double(rev)) == 0 % if specified load revision is not available then output an error message
            warning('Specified load revision number doesn''t exist for this location');
            load_condition = 0;
        end
    else
        [rev] = mysql(['select load_rev from loads_',loads.type,' where id=',location]);
        rev_load = max(str2double(rev));
    end
else % automatically selection of revision for soil, structure and loads
    if rev_global == 1000 % then detect the latest revision no. for the specified location
        [rev] = mysql(['select rev from configurations where id=',location]);
        rev_global = max(str2double(rev));
    else % check, if specified global revision is available for specified location
        [rev] = mysql(['select rev from configurations where id=',location]);
%         if ismember(rev_global,str2double(rev)) == 0 % if specified global revision is not available then output an error message
%             error('Specified global revision number doesn''t exist for this location');
%         end
    end
    [rev_structure, rev_load_ULS, rev_load_FLS, rev_soil] = mysql(['select rev_MP, rev_load_ULS, rev_load_FLS, soil_data_detail_rev from configurations where id=',...
        location,' and rev=',num2str(rev_global)]);
    rev_structure = str2double(rev_structure);
    rev_soil = str2double(rev_soil);
    % check which load.type shall be used (ULS or FLS)
    if strcmp(loads.type,'ULS')
        rev_load = str2double(rev_load_ULS);
    elseif strcmp(loads.type,'FLS')
        rev_load = str2double(rev_load_FLS);
    end
end

% update of rev no.s specified by user
data.revision.global = rev_global;
data.revision.soil = rev_soil;
data.revision.structure = rev_structure;
data.revision.loads = rev_load;

if rev_soil<10
    rev_soil        = ['''0',num2str(rev_soil),'''']; % revision no. of soil parameters to be used
else
    rev_soil        = ['''',num2str(rev_soil),'''']; % revision no. of soil parameters to be used
end
if rev_structure<10
    rev_structure   = ['''0',num2str(rev_structure),'''']; % revision no. of structure to be used
else
    rev_structure   = ['''',num2str(rev_structure),'''']; % revision no. of structure to be used
end
if rev_load<10
    rev_load        = ['''0',num2str(rev_load),'''']; % revision no. of loads to be used
else
    rev_load        = ['''',num2str(rev_load),'''']; % revision no. of loads to be used
end
%%

%--------------------------------------------------------------------------
% Water depths and scour (table: water_depths)
%--------------------------------------------------------------------------
table = 'water_depths';
[water_depth,scour_local,scour_ord] = mysql(['select water_depth,',...
    'scour_local,scour_ord from ',table,' where id= ''',location(4:end)]);

scour.local = scour_local;
 
scour.ORD = scour_ord; 

scour.water_depth = abs(water_depth);

if strcmp(soil.type,'UB')
	scour.local = 0;
end
%%
%--------------------------------------------------------------------------
% Structure and steel (table: mp_cans and mp_main_dimensions)
%--------------------------------------------------------------------------
% getting pile dimensions (table: mp_main_dimensions)
table = 'mp_main_dimensions';
[pile_top, pile_tip] = mysql(['select pile_top, pile_tip from ',table,...
    ' where id=',location,' and structural_revision=',rev_structure]);
% getting elevation of top of first soil layer [m VREF] from db
table = 'soil_data_detail';
[top_level] = mysql...
    (['select top_level from ',table,' where id=',location_geo,...
    ' and layer =',num2str(1),' and soil_revision=',rev_soil]);


% definiton of pile main dimensions:
pile.tip = pile_tip;
pile.top = pile_top;
pile.head                       = 0; % [m VREF] pile head is set to be at mudline (internal COSPIN definition)
% pile.length_above_mudline       = pile_top + scour.water_depth; % [m] pile length above mudline (from top to mudline) 
if plots.pilehead_vs_length     == 0 % for general calculations only one pile embedment length is considered
    pile.length_start_or           = -scour.water_depth - pile_tip; % [m] original embedment length (pile_tip is defined in mLAT)
    
    if  settings.SSI2db ==1 
        if settings.multi_toe_levels == 1 % allows writing of toe springs for multiple levels in batch mode
           pile.length_start_or = 60;
           disp ('Pile lenght of 60 m is considered to provide p-y and/or t-z curves')
        elseif settings.multi_toe_levels == 2
           pile.length_start_or = pile.multi_length_start_or;
        else
            pile.length_start_or = 60;
            disp ('Pile lenght of 60 m is considered to provide p-y and/or t-z curves')
        end
       pile.length_start = pile.length_start_or ;% %increasing pile length for p-y and t-z curves so the STR team dont need updated curves for updated  pile lenght
    else  
       pile.length_start = pile.length_start_or;
    end
    
    pile.length_end             = pile.length_start; % [m] for general calculations only one pile embedment length is considered
    pile.length_inc             = 1; % [m] is not used for general calculations
  
elseif plots.pilehead_vs_length == 1 % determination of critical pile embedment length
    pile.length_start_or           = -scour.water_depth - pile_tip;
    if strcmp(data.location,'A04')
        pile.length_start               = 45;
        pile.length_end                 = 65;
        pile.length_inc                 = 1;
    elseif strcmp(data.location,'G01')
        pile.length_start               = 43;
        pile.length_end                 = 63;
        pile.length_inc                 = 1;
     elseif strcmp(data.location,'G04')
        pile.length_start               = 40;
        pile.length_end                 = 60;
        pile.length_inc                 = 1;
    else
        pile.length_start               = 45;
        pile.length_end                 = 65;
        pile.length_inc                 = 1;        
    end
end
pile.length                     = (pile.length_start:pile.length_inc:...
                                    pile.length_end); % [m] creates a vector for rotation vs. embedment if needed
pile.stick_up                   = abs(pile.head-top_level)+scour.local; % [m] calculate pile stick-up

%%

% pile geometry (table: mp_cans)
% getting pile geometry from db
table = 'mp_cans';
[can_id, top_od,bottom_od,height,wall_thickness,steel_grade] = mysql...
    (['select can_id, top_od, bottom_od, height, wall_thickness,steel_grade from ',table,' where id=',...
    location,' and structural_revision=',rev_structure,' ORDER BY can_id']);
	
pile.geometryTable = num2cell([can_id, top_od, bottom_od, height, wall_thickness/1000]);	


% looking for can that enters the seabed:
% the ref level is not necessarily LAT depends on the data in database
can_top_elevation(1,1) = pile_top; % [mLAT] elevation of each can above sealevel (negative below water level)
can_top_elevation_mud(1,1) = pile_top + scour.water_depth; % [m VREF] elevation of each can above mudline (negative below mudline)
for i = 2:size(can_id,1)
    can_top_elevation(i,1) = can_top_elevation(i-1) - height(i-1); % [mLAT] elevation of each can (-top) above sealevel
    can_top_elevation_mud(i,1) = can_top_elevation(i)+ scour.water_depth; % [m VREF] elevation of each can (-top) above mudline
end
% look for last can, which begins above or at mudline
i = 1;
while can_top_elevation(i) >= -scour.water_depth % scour.water_level is defined as positive value
    i = i+1;
end
i = i-1;    % can that enters the seabed

% definition of pile geometry:
k = 1;

pile.cross_section.diameter(k) = bottom_od(i); % bottom pile diameter of that can, which enters the seabed is used

% definition of cross section:

pile.cross_section.toplevel(k)     = 0;
pile.cross_section.thickness(k)    = wall_thickness(i) / 1000; % [m] wallthickness in db is in [mm]
% exclude pile shoe for determination of critical pile length
if plots.pilehead_vs_length == 1
    last_can = size(can_top_elevation_mud,1)-1;
else
    last_can = size(can_top_elevation_mud,1);
end
for j = (i+1):last_can
    % wallthickness will only be saved, if there is a change in wall thickness
    if wall_thickness(j)/1000 ~= pile.cross_section.thickness(k,1) 
        pile.cross_section.toplevel(k+1,1)     = can_top_elevation_mud(j); % [m VREF]
        pile.cross_section.thickness(k+1,1)    = wall_thickness(j) / 1000; % % [m] wallthickness in db is in [mm]
        pile.cross_section.diameter(k+1,1) = top_od(j); % pile diameter of that can, which enters the seabed is used

        k = k+1;
    end
end

% definition of pile type:
settings.pile_type              = 'open'; % only open ended piles are used at Merkur
if strcmp(settings.pile_type,'open')
    pile.cross_section.endarea  = pi*((pile.cross_section.diameter(end)/2)^2-(pile.cross_section.diameter(end)/2-pile.cross_section.thickness(end))^2); % [m^2] end bearing area
%     pile.cross_section.K        = 0.8;
elseif strcmp(settings.pile_type,'closed')
    pile.cross_section.endarea  = pi*(pile.cross_section.diameter(end)/2)^2;
%     pile.cross_section.K        = 1.0;
    msgbox('Pile is closed-ended. For correct axial capacity calculation please set reducion.skin_inner = 0.0 in the file reduction_factors.m as this will exclude internal shaft friction')
end  

pile.diameter = pile.cross_section.diameter(1,1);

% extraction of steel properties:
% steel grade of can, which is 2 cans below seabed will be used for the whole
% pile (sometimes can which enters the seabed has different steel properties)
temp = regexp(steel_grade{i+2},'\d*','match'); % extracts the number of steel grade from DB (example: S355NL/ML-Z35 -> [355, 35])
steel = str2double(temp{1});    % converts from cell to double
switch steel    % definition of different steel types
    case 355    % S355NL
        pile.density         = 78; % [kN/m^3] density of pile material (steel)
        pile.sigma_y         = 355000; % [kPa] characteristic steel yield strength 
        pile.E               = 2.1e8; % [kPa] Young's modulus of steel 
        pile.G               = 8.1e7; % [kPa] Shear modulus 
        pile.ksf             = 0.5652; % [-] Shear correction factor for Timoshenko beam theory
    case 420    % S420NL
        disp('No steel properties are defined for S420');
end

%%
%--------------------------------------------------------------------------
% General soil properties (table: soil_data_detail)
%--------------------------------------------------------------------------
table = 'soil_data_detail';
variables = ['layer, top_level, bottom_level, description, legend_code, gamma_eff, gamma_eff_LB, gamma_eff_UB, c_u, c_eff, c_eff_LB, c_eff_UB,'...
    ' delta_c_u, c_u_LB, c_u_UB, delta_c_u_LB, delta_c_u_UB, phi, phi_LB, phi_UB, delta_eff, delta_eff_LB,'...
    ' delta_eff_UB, E_s, delta_Es, delta_Es_LB, delta_Es_UB, E_s_LB, E_s_UB, K_0, epsilon_50, epsilon_50_LB, epsilon_50_UB, delta_epsilon_50, delta_epsilon_50_LB, delta_epsilon_50_UB, J,'...
    ' limit_skin, limit_skin_LB, limit_skin_UB, limit_alpha, limit_tip, limit_tip_LB, limit_tip_UB, G_0, G_0_LB, G_0_UB, poisson,'...
    ' delta_G_0, delta_G_0_LB, delta_G_0_UB, q_ur_axial, q_ur_axial_LB, q_ur_axial_UB, q_ur_py, q_ur_py_LB, q_ur_py_UB, krm, RQD,'...
    ' Nq, Nq_LB, Nq_UB, t_multiplier, z_multiplier, p_multiplier, y_multiplier, y_multiplier_10000'];
strSQL = ['select ',variables,' from ',table,' where  id=',location_geo,...
    ' and soil_revision=',rev_soil];

[layer, top_level, bottom_level, description, legend_code, gamma_eff, gamma_eff_LB, gamma_eff_UB, c_u, c_eff, c_eff_LB, c_eff_UB,...
    delta_c_u, c_u_LB, c_u_UB, delta_c_u_LB, delta_c_u_UB, phi, phi_LB, phi_UB, delta_eff, delta_eff_LB,...
    delta_eff_UB, E_s, delta_Es, delta_Es_LB, delta_Es_UB, E_s_LB, E_s_UB, K_0, epsilon_50, epsilon_50_LB, epsilon_50_UB, delta_epsilon_50, delta_epsilon_50_LB, delta_epsilon_50_UB, J,...
    limit_skin, limit_skin_LB, limit_skin_UB, limit_alpha, limit_tip, limit_tip_LB, limit_tip_UB, G_0, G_0_LB, G_0_UB, poisson,...
    delta_G_0, delta_G_0_LB, delta_G_0_UB, q_ur_axial, q_ur_axial_LB, q_ur_axial_UB, q_ur_py, q_ur_py_LB, q_ur_py_UB, krm, RQD,...
    Nq, Nq_LB, Nq_UB, t_multiplier, z_multiplier, p_multiplier, y_multiplier, y_multiplier_10000] = mysql(strSQL); %#ok<ASGLU>

%%
% definition of soil properties:
soil.layer                  = layer;									
soil.legendcode             = legend_code;%[-]
soil.toplevel               = -abs(top_level); % [m VREF]
soil.bottomlevel            = -abs(bottom_level); % [m VREF]
soil.description            = description; % soil unit

soil.gamma_eff              = gamma_eff; % [kN/m^3]
soil.c_eff                  = c_eff;
soil.K0                     = K_0; % [kPa]
soil.delta_epsilon50        = delta_epsilon_50; % [1/m]
soil.J                      = J; % [-]
soil.limit_alpha            = limit_alpha; % [-]
soil.poisson                = poisson; % [-]
soil.k_rm                   = krm; % [-]
soil.RQD                    = RQD; % [-]
soil.degradation.value_tz_t = t_multiplier; % [-]soil.degradation.value_tz_t = t_multiplier; % [m]
soil.degradation.value_tz_z = z_multiplier; % [-]
soil.degradation.value_py_p = p_multiplier; % [-]

if loads.n_cycles>100
    soil.degradation.value_py_y = y_multiplier_10000; % [-]
else
    soil.degradation.value_py_y = y_multiplier; % [-]
end

if strcmp(loads.static_cyclic,'static')
    soil.degradation.value_py_p = ones(length(1),length(soil.degradation.value_py_p));
    soil.degradation.value_py_y = ones(length(1),length(soil.degradation.value_py_y));
end

    soil.degradation.value_py_p = ones(length(1),length(soil.degradation.value_py_p));
    soil.degradation.value_py_y = ones(length(1),length(soil.degradation.value_py_y));
    
if settings.PSI==0
    soil.q_ur_axial             = q_ur_axial_LB; % [kPa]
    soil.Nq                     = Nq_LB; % [-]
    soil.limit_tip              = limit_tip_LB;% [kPa]
end


soil.cu_LB                  = c_u_LB;
soil.delta_cu_LB            = delta_c_u_LB; % [kPa/m]

% definition of soil properties assumed differently for BE, LB and UB:

if strcmp(soil.type,'BE')
    
soil.cu                     = c_u; % [kPa]
soil.cu_axial               = c_u; % [kPa]
soil.delta_cu               = delta_c_u; % [kPa/m]
soil.delta_cu_axial         = delta_c_u; % [kPa/m]
soil.phi                    = phi; % [degrees]
soil.delta_eff              = delta_eff; % [degrees]
soil.Es                     = E_s*1000; % [kPa]
soil.delta_Es               = delta_Es*1000 ; %%to be changed
soil.epsilon50              = epsilon_50/100; % [-]
soil.limit_tip              = limit_tip; % [kPa]
soil.limit_skin             = limit_skin; % [kPa]
soil.G0                     = G_0*1000; % [kPa]
soil.delta_G0               = delta_G_0*1000; % [kPa/m] (not used)
soil.q_ur_py                = q_ur_py; % [kPa]
soil.q_ur_axial             = q_ur_axial; % [kPa]
soil.Nq                     = Nq; % [-]
elseif strcmp(soil.type,'LB')
soil.cu                     = c_u_LB; % [kPa]
soil.cu_axial               = c_u_LB; % [kPa]
soil.delta_cu               = delta_c_u_LB; % [kPa/m]
soil.delta_cu_axial         = delta_c_u_LB; % [kPa/m]
soil.phi                    = phi_LB; % [degrees]
soil.delta_eff              = delta_eff_LB; % [degrees]
soil.Es                     = E_s_LB*1000; % [kPa]
soil.delta_Es               = delta_Es*1000 ; %%to be changed
soil.epsilon50              = epsilon_50_LB/100; % [-]
soil.limit_tip              = limit_tip_LB; % [kPa]
soil.limit_skin             = limit_skin_LB; % [kPa]
soil.G0                     = G_0_LB*1000; % [kPa]
soil.delta_G0               = delta_G_0_LB*1000; % [kPa/m] (not used)
soil.q_ur_py                = q_ur_py_LB; % [kPa]
soil.q_ur_axial             = q_ur_axial_LB; % [kPa]
soil.Nq                     = Nq_LB; % [-]
elseif strcmp(soil.type,'UB')
soil.cu                     = c_u_UB; % [kPa]
soil.cu_axial               = c_u_UB; % [kPa]
soil.delta_cu               = delta_c_u_UB; % [kPa/m]
soil.delta_cu_axial         = delta_c_u_UB; % [kPa/m]
soil.phi                    = phi_UB; % [degrees]
soil.delta_eff              = delta_eff_UB; % [degrees]
soil.Es                     = E_s_UB*1000; % [kPa]
soil.delta_Es               = delta_Es*1000 ; %%to be changed
soil.epsilon50              = epsilon_50_UB/100; % [-]
soil.limit_tip              = limit_tip_UB; % [kPa]
soil.limit_skin             = limit_skin_UB; % [kPa]
soil.G0                     = G_0_UB*1000; % [kPa]
soil.delta_G0               = delta_G_0_UB*1000; % [kPa/m] (not used)
soil.q_ur_py                = q_ur_py_UB; % [kPa]
soil.q_ur_axial             = q_ur_axial_UB; % [kPa]
soil.Nq                     = Nq_UB; % [-]

end


% use of soil parameters including psf
if soil.psf == 1    % then use of soil parameters including psf
    [SF reduction]  = factors(loads); %#ok<ASGLU>
    soil.phi        = atand(tand(soil.phi) / SF.drained);
    soil.delta_eff  = atand(tand(soil.delta_eff) / SF.drained);
    soil.cu         = soil.cu / SF.undrained; % [kPa]
    soil.c_eff      = soil.c_eff / SF.undrained; % [kPa]
    soil.delta_cu   = soil.delta_cu / SF.undrained; % [kPa]
    soil.q_ur_py    = q_ur_py/ SF.undrained; % [kPa]
end
%% p-y and t-z curves
%--------------------------------------------------------------------------
for i = 1:size(soil.description,1)
    if contains(soil.description{i},{'1','2','3'})
        % Bothkennar
        %soil.model_py{i} = 'PISA Bothkennar clay'; % [text]
        soil.model_py{i} = 'PISA Upper clay'; % [text]
		%soil.model_py{i} = 'Kirsch soft clay'; % [text]
        soil.model_axial{i} = 'API clay'; % [text]
        
    elseif contains(soil.description{i},{'4'})
        % Cowden
        %soil.model_py{i} = 'PISA Cowden clay'; % [text]
        soil.model_py{i} = 'PISA Till'; % [text]
		%soil.model_py{i} = 'Kirsch stiff clay'; % [text]
        soil.model_axial{i} = 'API sand'; % [text]
    else
		soil.model_py{i} = 'PISA Chalk'; % [text]
		%soil.model_py{i} = 'Kirsch stiff clay'; % [text]
        soil.model_axial{i} = 'API sand'; % [text]
	end
end

%%
%--------------------------------------------------------------------------
% Loads
%--------------------------------------------------------------------------
if load_condition
    % getting loads from db
    if loads.type == 'ULS' %#ok<*BDSCA>
        table = 'loads_ULS';
    elseif loads.type == 'FLS'
        table = 'loads_FLS';
    end

    if strcmp(settings.interface,'UNFAC') 
    variables = 'elev,FXY_geo,MXY_geo,FZ_geo,MZ_geo';
    strSQL = ['select ',variables,' from ',table,' where  id=',location,...
        ' and load_rev=',rev_load,' ORDER BY elev DESC'];
    [elev,FXY,MXY,FZ,MZ] = mysql(strSQL);

    else 
    variables = 'elev,FXY,MXY,FZ,MZ';
    strSQL = ['select ',variables,' from ',table,' where  id=',location,...
        ' and load_rev=',rev_load,' ORDER BY elev DESC'];
    [elev,FXY,MXY,FZ,MZ] = mysql(strSQL);
    end



    if settings.PSI == 1
        disp('Loads are changed, since settings.PSI = 1')
        loads.H = 10;
        loads.M = -100;
        loads.Vc = 10;
        loads.Vt = 10;
        loads.Mz = 0;
    else
    % interpolation of loads to specified pile head level:
    % look for last/first cross section (for loads) above/below mudline
        i = 1;
        while elev(i) >= (-scour.water_depth + pile.head) % scour.water_level is defined as positive value, pile.head in [m VREF]
        i = i+1;
        end
        i = i - 1;  % last cross section (for loads) above mudline
        j = i + 1;  % first cross section (for loads) below mudline
    % interpolation:
        factor   = 1 / (elev(j) - elev(i)) * ((-scour.water_depth + pile.head) - elev(i)); % factor for interpolation
        loads.H  = abs(FXY(i) + (FXY(j) - FXY(i)) * factor); % [kN] FXY at pile head, characteristic load
        loads.M  = -abs(MXY(i) + (MXY(j) - MXY(i)) * factor); % [kNm] MXY at pile head, moment must have opposite sign as FXY, since it is MONOPILE design at Merkur, characteristic load
    % Axial loading, Attention: DESIGN Load incl psf
        if FZ(i) < 0 % case for tension
        loads.Vc = 0; % [kN] FZ (compression) at pile head
        loads.Vt = abs(FZ(i) + (FZ(j) - FZ(i)) * factor); % [kN] FZ (tension) at pile head
        else % case for compression
        loads.Vc = abs(FZ(i) + (FZ(j) - FZ(i)) * factor); % [kN] FZ (compression) at pile head
        loads.Vt = 0; % [kN] FZ (tension) at pile head
        end
        loads.Mz = MZ(i) + (MZ(j) - MZ(i)) * factor; % [kNm] MZ at pile head, Attention: DESIGN Load incl. psf
    end
end
%%
% Close MySQL-database

% if plots.utilization_ratio
% 
%    loads.H_load_10= load(['output\rev0.1\mat_files\H_load_10_pos_',data.location,'.mat']);
%    loads.H=loads.H_load_10.H_load_10;
%    disp(['Lateral load taken at pile deflection equal to 10% of diameter'])
% 
% end
    
mysql('close')

