function [pile] = interimGeometry (plots,id,soil,loads,settings,pile,scour)
% CHANGE: FKMV - 27/02/2020 - loop over different locations added 
if strcmp('A04',id)
    length_start_original = 48.5;  % specifications of different pile lengths for each location
elseif strcmp('G01',id)
    length_start_original = 47;  % specifications of different pile lengths for each location
elseif strcmp('G04',id)
    length_start_original = 44;  % specifications of different pile lengths for each location
else
    length_start_original = 48.5;
    disp ('error - pile length for given location not specified, but assigned as A04 as it does not matter for SSI springs')
end
pile.diameter                   = 9.2;

    

if settings.axial_loading == 0
W_T=80/1000; % % [m] wallthickness in db is in [mm]
elseif settings.axial_loading == 1
W_T1=80/1000; % % [m] wallthickness in db is in [mm]
W_T2=100/1000; % % [m] wallthickness in db is in [mm]
end


if plots.pilehead_vs_length     == 0 % for general calculations only one pile embedment length is considered
%     pile.length_start_or           = 51; % [m] original embedment length (pile_tip is defined in mLAT)
    pile.length_start_or           = length_start_original; % [m] original embedment length (pile_tip is defined in mLAT)
    
    if  settings.SSI2db ==1 
       pile.length_start_or           = 60; % [m] original embedment length (pile_tip is defined in mLAT) %%%%%%%%%%%%%%%%%%%%%%%%%
       pile.length_start = pile.length_start_or ; 
       disp ('An additional pile lenght is considered to provide p-y and/or t-z curves')
    elseif plots.res_vs_pilelength         == 1
       pile.extra_L=  10;
        pile.length_start_or           = pile.length_start_or+pile.extra_L; % [m] original embedment length (pile_tip is defined in mLAT) %%%%%%%%%%%%%%%%%%%%%%%%%
       pile.length_start = pile.length_start_or ;   
        
    else
       pile.length_start = pile.length_start_or;
    end
    
    pile.length_end             = pile.length_start; % [m] for general calculations only one pile embedment length is considered
    pile.length_inc             = 1; % [m] is not used for general calculations
  
pile_tip=-scour.water_depth-pile.length_start_or;
pile.length                     = (pile.length_start:pile.length_inc:pile.length_end);
pile.tip = pile_tip; 
end

if settings.axial_loading == 0
pile.cross_section.toplevel     = 0;
pile.cross_section.thickness    = W_T; % % [m] wallthickness in db is in [mm]
pile.cross_section.diameter= [pile.diameter];


elseif settings.axial_loading == 1 
pile.cross_section.diameter= pile.diameter;
 disp ('Pile toe can WT is 100 mm')
  
pile.cross_section.toplevel     = [0;scour.water_depth+pile_tip+2.5];
pile.cross_section.thickness    = [W_T1;W_T1]; % % [m] wallthickness in db is in [mm]
pile.cross_section.diameter= [pile.diameter;pile.diameter];
pile.cross_section.thickness_toe= W_T2;


end

if settings.axial_loading == 1
if strcmp(settings.pile_type,'open')
    pile.cross_section.endarea  = pi*((pile.cross_section.diameter(end)/2)^2-(pile.cross_section.diameter(end)/2-pile.cross_section.thickness_toe)^2); % [m^2] end bearing area
%     pile.cross_section.K        = 0.8;
elseif strcmp(settings.pile_type,'closed')
    pile.cross_section.endarea  = pi*(pile.cross_section.diameter(end)/2)^2;
%     pile.cross_section.K        = 1.0;
    msgbox('Pile is closed-ended. For correct axial capacity calculation please set reducion.skin_inner = 0.0 in the file reduction_factors.m as this will exclude internal shaft friction')
end  
end
        
end