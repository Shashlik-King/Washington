function [kstop ksbot] = secspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii)
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
    [kstop ksbot] = apisandsec_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom/element.degradation_tz_z(i),A,fs,i,ii);
    
elseif strcmp(element.model_axial(i),'API clay')
    [kstop ksbot] = apiclaysec_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom/element.degradation_tz_z(i),A,fs,i,ii);
    
elseif strcmp(element.model_axial(i),'Zero soil')
    kstop         = 0;
    ksbot         = 0;
    
else
    error('The specified soil model is not supported')
    
end

kstop = kstop*element.degradation_tz_t(i)/element.degradation_tz_z(i);
ksbot = ksbot*element.degradation_tz_t(i)/element.degradation_tz_z(i);