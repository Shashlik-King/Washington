function Q = apiclayQz(element,pile,plug_unplug,zbot,qp,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the mobilised end bearing resistance using the Q-z curve from API.
% 
% INPUT:  element       : Element information about soil 
%         pile          : Pile information for each element
%         z_bot         : Displacement in the bottom node [m]
%         plug_unplug   : Determines wether the pile is plugged or
%                         unplugged
%         i             : Counter referring to element number
%
% OUTPUT: Q             : Mobilised end bearing resistance [kN]
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