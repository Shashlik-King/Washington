function Q_end = Qz_spring(element,pile,plug_unplug,zbot,qp,i)
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

% Determination of the model to apply

if strcmp(element.model_axial(i),'API sand')
    Q = apisandQz(element,pile,plug_unplug,zbot,qp,i);
 
elseif strcmp(element.model_axial(i),'API clay')
    Q = apiclayQz(element,pile,plug_unplug,zbot,qp,i);
    
elseif strcmp(element.model_axial(i),'Zero soil')
    Q = 0;
    
else
    error('The specified soil model is not supported')
    
end

if isnan(Q)
    Q_end = 0; 
else
    Q_end=Q;
end 
end 

