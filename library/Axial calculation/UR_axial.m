function output = UR_axial(element,settings,pile,loads,output,reduction,plug_unplug,A,fs,ii)
%% UTILIZATION RATIO
% CALCULATES THE AXIAL UTILIZATION RATIO FOR THE SOIL
%--------------------------------------------------------------------------
% CHANGE LOG
% 22.11.2016    EBGX - Programming 
%--------------------------------------------------------------------------
%% Pre-allocation
%--------------------------------------------------------------------------
nelem           = element.nelem;
output.fsu_UR   = zeros(nelem-1,2); 
output.fs_UR    = zeros(nelem-1,2); 
D               = pile.diameter;
%--------------------------------------------------------------------------
%% Calculation routine 
%--------------------------------------------------------------------------
for c = 1:nelem-1
    z_topbottom = [output.deflections(3*c-1,min(output.n_possible,settings.n_max)) output.deflections(3*c+2,min(output.n_possible,settings.n_max))];
    % Calculation of mobilised shaft resistance 
    [kstop ksbot]      = secspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,c,ii);
    output.fs_UR(c,1)  = kstop*z_topbottom(1); 
    output.fs_UR(c,2)  = ksbot*z_topbottom(2); 
    % Calculation of ultimate shaft resistance 
    if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
        output.fsu_UR(c,1)   = fs.o(c,1)*D*pi;          % Ultimate outer shaft resistance [kN/m] in top
    elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
        output.fsu_UR(c,1)   = fs.o(c,1)*D*pi+fs.o(c,1)*(A.si(c,1)/A.so(c,1))*D*pi*reduction.skin_inner;
    elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
        output.fsu_UR(c,1)   = fs.o(c,1)*D*pi;
    elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
        output.fsu_UR(c,1)   = fs.o(c,1)*D*pi+fs.o(c,1)*(A.si(c,1)/A.so(c,1))*D*pi*reduction.skin_inner;
    end 
    if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
        output.fsu_UR(c,2)   = fs.o(c,2)*D*pi;          % and bottom of pile segment, determined in skin_tip.m
    elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
        output.fsu_UR(c,2)   = fs.o(c,2)*D*pi+fs.o(c,2)*(A.si(c,2)/A.so(c,2))*D*pi*reduction.skin_inner;
    elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
        output.fsu_UR(c,2)   = fs.o(c,2)*D*pi;
    elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
        output.fsu_UR(c,2)   = fs.o(c,2)*D*pi+fs.o(c,2)*(A.si(c,2)/A.so(c,2))*D*pi*reduction.skin_inner;
    end    
 
    if z_topbottom(1) > element.zpeak(c) 
        output.fsu_UR(c,1) = output.fs_UR(c,1); 
    end 
    if z_topbottom(2) > element.zpeak(c) 
        output.fsu_UR(c,2) = output.fs_UR(c,2); 
    end 

%% Calculation of utilization ratio for plot 
output.UR_axial(c,:) = output.fs_UR(c,:)/output.fsu_UR(c,:);
if output.UR_axial(c,:) == 0
    disp('Warning of rank deficiency due to scour')
end
end 