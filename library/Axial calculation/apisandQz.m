function  Q = apisandQz(element,pile,plug_unplug,zbot,qp,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying t-z curves.
% 
% INPUT:  element       : Element information about soil 
%         element.degradation: Degradation of t-z curve to take cyclic into
%                              account 
%         pile          : Pile information for each element
%         z_topbottom   : Vertical displacement in top and bottom of
%                         element
%         i             : Counter referring to element number
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    24.10.2016   Programming
%--------------------------------------------------------------------------

zpeak_Qz = element.zpeak_Qz(i); % Displacement to maximum tip resitance [m]

if strcmp(plug_unplug.comp_index,'Plugged')
    Qmax = qp(i,2)*pi/4*pile.diameter^2;
elseif strcmp(plug_unplug.comp_index,'Unplugged')
    Qmax = qp(i,2)*pile.cross_section.endarea;
end 

z = [0 0.02 0.13 0.42 0.73 1]*zpeak_Qz; 

Q_API = [0 0.25 0.5 0.75 0.9 1]*Qmax; 
if zbot > zpeak_Qz 
    Q = Qmax;
else 
    Q=interp1(z,Q_API,zbot); 
end 
end 
