function [kttop ktbot] = tanspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom 
% of each pile segment by applying t-z curves.
% 
% INPUT:  element       : Element information about soil 
%         element.degradation: Degradation of t-z curve to take cyclic into
%                              account 
%         pile          : Pile information for each element
%         z_topbottom   : Vertical displacement in top and bottom of
%                         element
%         i             : Counter referring to element number
%
% OUTPUT: kttop   : Soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    24.10.2016   Programming
%--------------------------------------------------------------------------
% Determination of the model to apply

if strcmp(element.model_axial(i),'API sand')
    [kttop ktbot]   = apisandtan_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom/element.degradation_tz_z(i),A,fs,i,ii);
    
elseif strcmp(element.model_axial(i),'API clay')
    [kttop ktbot]   = apiclaytan_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom/element.degradation_tz_z(i),A,fs,i,ii);
    
elseif strcmp(element.model_axial(i),'Zero soil')
    kttop           = 0;
    ktbot           = 0;
    
else
   error('The specified soil model is not supported')
   
end
kttop = kttop*element.degradation_tz_t(i)/element.degradation_tz_z(i);
ktbot = ktbot*element.degradation_tz_t(i)/element.degradation_tz_z(i);