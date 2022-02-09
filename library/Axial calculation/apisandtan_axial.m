function [kttop ktbot] = apisandtan_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii)
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
%--------------------------------------------------------------------------

% -------- Initializing pile and soil parameters --------------------------

ztop    = z_topbottom(1);            % Vertical disp. at the top of the
                                    % pile segment [m]

zbot    = z_topbottom(2);            % Vertical disp. at the bottom of the
                                    % pile segment [m]

D       = pile.diameter;            % Outer diameter [m]
zpeak   = element.zpeak(i);                   % Displacement to maximum shaft friction according to API (2011) 
zres    = element.zres(i);          % Displacement to residual value given as zres/zpeak [-]

if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') ==1
    tmaxtop   = fs.o(i,1)*D*pi;          % Ultimate outer shaft resistance [kN/m] in top
elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') ==1
    tmaxtop   = fs.o(i,1)*D*pi+fs.o(i,1)*(A.si(i,1)/A.so(i,1))*D*pi*reduction.skin_inner;
elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') ==1
    tmaxtop   = fs.o(i,1)*D*pi;
elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') ==1
    tmaxtop   = fs.o(i,1)*D*pi+fs.o(i,1)*(A.si(i,1)/A.so(i,1))*D*pi*reduction.skin_inner;
end 
if strcmp(plug_unplug.comp_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Comp') ==1
    tmaxbot   = fs.o(i,2)*D*pi;          % and bottom of pile segment, determined in skin_tip.m
elseif strcmp(plug_unplug.comp_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Comp') ==1
    tmaxbot   = fs.o(i,2)*D*pi+fs.o(i,2)*(A.si(i,2)/A.so(i,2))*D*pi*reduction.skin_inner;
elseif strcmp(plug_unplug.tens_index,'Plugged') && strcmp(settings.analysis_loading(ii),'Tens') ==1
    tmaxbot   = fs.o(i,2)*D*pi;
elseif strcmp(plug_unplug.tens_index,'Unplugged') && strcmp(settings.analysis_loading(ii),'Tens') ==1
    tmaxbot   = fs.o(i,2)*D*pi+fs.o(i,2)*(A.si(i,2)/A.so(i,2))*D*pi*reduction.skin_inner;
end  

tres      =element.tres(i);              % Residual value (t=tres*tmax) of the t-z curve for each element, defined in input [-]

z = [0 0.16 0.31 0.57 0.80 1 zres 10]*zpeak;

t_top = [0.00 0.30 0.50 0.75 0.90 1.00 tres tres]*tmaxtop;
t_bot = [0.00 0.30 0.50 0.75 0.90 1.00 tres tres]*tmaxbot;                

% -------- Tangent stiffness, kttop and ktbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dt/dz [kN/m^2]. 
% Similar for ktbot
% Expressions based on API(2011)

% kttop

if ztop>=z(1) && ztop<=z(2)
    kttop=(t_top(2)-t_top(1))/(z(2)-z(1)); 
elseif ztop>z(2) && ztop<=z(3)
    kttop=(t_top(3)-t_top(2))/(z(3)-z(2)); 
elseif ztop>z(3) && ztop<=z(4)
    kttop=(t_top(4)-t_top(3))/(z(4)-z(3));
elseif ztop>z(4) && ztop<=z(5)
    kttop=(t_top(5)-t_top(4))/(z(5)-z(4));
elseif ztop>z(5) && ztop<=z(6) 
    kttop=(t_top(6)-t_top(5))/(z(6)-z(5));
elseif ztop>z(6) && ztop<=z(7)
    kttop=(t_top(7)-t_top(6))/(z(7)-z(6));
elseif ztop>z(7) && ztop<=z(8)
    kttop=(t_top(8)-t_top(7))/(z(8)-z(7));
elseif ztop>z(8)
    kttop=0;
    %error('t-z curve not defined for vertical displacement')
elseif ztop<z(1) && ztop>=-z(2)
    kttop=-(t_top(2)-t_top(1))/(z(2)-z(1)); 
elseif ztop<-z(2) && ztop>=-z(3)
    kttop=-(t_top(3)-t_top(2))/(z(3)-z(2)); 
elseif ztop<-z(3) && ztop>=-z(4)
    kttop=-(t_top(4)-t_top(3))/(z(4)-z(3));
elseif ztop<-z(4) && ztop>=-z(5)
    kttop=-(t_top(5)-t_top(4))/(z(5)-z(4));
elseif ztop<-z(5) && ztop>=-z(6)
    kttop=-(t_top(6)-t_top(5))/(z(6)-z(5));
elseif ztop<-z(6) && ztop>=-z(7)
    kttop=-(t_top(7)-t_top(6))/(z(7)-z(6));
elseif ztop<-z(7) && ztop>=-z(8)
    kttop=-(t_top(8)-t_top(7))/(z(8)-z(7));
elseif ztop<-z(8)
    kttop=0;
%     %error('t-z curve not defined for vertical displacement')
end 

% ktbot

if zbot>=z(1) && zbot<=z(2)
    ktbot=(t_bot(2)-t_bot(1))/(z(2)-z(1)); 
elseif zbot>z(2) && zbot<=z(3)
    ktbot=(t_bot(3)-t_bot(2))/(z(3)-z(2)); 
elseif zbot>z(3) && zbot<=z(4)
    ktbot=(t_bot(4)-t_bot(3))/(z(4)-z(3));
elseif zbot>z(4) && zbot<=z(5)
    ktbot=(t_bot(5)-t_bot(4))/(z(5)-z(4));
elseif zbot>z(5) && zbot<=z(6)
    ktbot=(t_bot(6)-t_bot(5))/(z(6)-z(5));
elseif zbot>z(6) && zbot<=z(7)
    ktbot=(t_bot(7)-t_bot(6))/(z(7)-z(6));
elseif zbot>z(7) && zbot<=z(8)
    ktbot=(t_bot(8)-t_bot(7))/(z(8)-z(7));
elseif zbot>z(8)
    ktbot=0;
%     %error('t-z curve not defined for vertical displacement')
elseif zbot<z(1) && zbot>=-z(2)
    ktbot=-(t_bot(2)-t_bot(1))/(z(2)-z(1)); 
elseif zbot<-z(2) && zbot>=-z(3)
    ktbot=-(t_bot(3)-t_bot(2))/(z(3)-z(2)); 
elseif zbot<-z(3) && zbot>=-z(4)
    ktbot=-(t_bot(4)-t_bot(3))/(z(4)-z(3));
elseif zbot<-z(4) && zbot>=-z(5)
    ktbot=-(t_bot(5)-t_bot(4))/(z(5)-z(4));
elseif zbot<-z(5) && zbot>=-z(6)
    ktbot=-(t_bot(6)-t_bot(5))/(z(6)-z(5));
elseif zbot<-z(6) && zbot>=-z(7)
    ktbot=-(t_bot(7)-t_bot(6))/(z(7)-z(6));
elseif zbot<-z(7) && zbot>=-z(8)
    ktbot=-(t_bot(8)-t_bot(7))/(z(8)-z(7));
elseif zbot<-z(8)
    ktbot=0;
%     %error('t-z curve not defined for vertical displacement')    
end