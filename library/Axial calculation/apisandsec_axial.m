function [kstop ksbot] = apisandsec_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying t-z curves according to API(2011).
% 
% INPUT:  fs.o   : Outer shaft resistance [kN/m^2]
%         element.tres  : Residual value on t-z curve [-]
%         pile.diameter : Outer pile diameter [m],
%         z_topbottom   : Vertical displacement in top and bottom of
%                         element
%         i             : Counter referring to element number
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : EBGX
% APPROVED        : 
%
% LAST MODIFIED   : EBGX    24.10.2016   Programming
%                   
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------

ztop    = z_topbottom(1);            % Vertical disp. at the top of the
                                    % pile segment [m]

zbot    = z_topbottom(2);            % Vertical disp. at the bottom of the
                                    % pile segment [m]

D       = pile.diameter;            % Outer diameter [m]
zpeak   = element.zpeak(i);                   % Displacement to maximum shaft resistance according to API (2011) 
zres    = element.zres(i);          % Displacement to residual value given as zres/zpeak [-]

if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
    tmaxtop   = fs.o(i,1)*D*pi;          % Ultimate outer shaft resistance [kN/m] in top
elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
    tmaxtop   = fs.o(i,1)*D*pi+fs.o(i,1)*(A.si(i,1)/A.so(i,1))*D*pi*reduction.skin_inner;
elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
    tmaxtop   = fs.o(i,1)*D*pi;
elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
    tmaxtop   = fs.o(i,1)*D*pi+fs.o(i,1)*(A.si(i,1)/A.so(i,1))*D*pi*reduction.skin_inner;
end 
if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
    tmaxbot   = fs.o(i,2)*D*pi;          % and bottom of pile segment, determined in skin_tip.m
elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') == 1
    tmaxbot   = fs.o(i,2)*D*pi+fs.o(i,2)*(A.si(i,2)/A.so(i,2))*D*pi*reduction.skin_inner;
elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
    tmaxbot   = fs.o(i,2)*D*pi;
elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') == 1
    tmaxbot   = fs.o(i,2)*D*pi+fs.o(i,2)*(A.si(i,2)/A.so(i,2))*D*pi*reduction.skin_inner;
end    

tres      =element.tres(i);              % Residual value (t=tres*tmax) of the t-z curve for each element [-]

z = [0 0.16 0.31 0.57 0.80 1 zres 10]*zpeak;                  % Definition of t-z curve according to API(2011)

t_top = [0.00 0.30 0.50 0.75 0.90 1.00 tres tres]*tmaxtop; % Definition of t-z curve according to API(2011)
t_bot = [0.00 0.30 0.50 0.75 0.90 1.00 tres tres]*tmaxbot; % Definition of t-z curve according to API(2011)                                                                  

% -------- Secant stiffness, kstop and ksbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness t/z [kN/m^2] for z >0. For z = 0
% kstop is the tangent stiffness between origo and the first point defined 
% on the t-z curve. 
% Similar for kbot

if ztop == 0;
    % Expression based on secant stiffness between origo and first point
    % defined on t-z curve by API(2011)
    kstop = (t_top(2)-t_top(1))/(z(2)-z(1)); 
elseif ztop>0 
    t=interp1(z,t_top,ztop);
    kstop=t/ztop;
elseif ztop<0
    t=interp1(-z,-t_top,ztop);
    kstop=t/ztop;
end


if zbot == 0;
    ksbot = (t_bot(2)-t_bot(1))/(z(2)-z(1));
elseif zbot>0 
    t=interp1(z,t_bot,zbot);
    ksbot=t/zbot;
elseif zbot<0
    t=interp1(-z,-t_bot,zbot);
    ksbot=t/zbot;
end

end 