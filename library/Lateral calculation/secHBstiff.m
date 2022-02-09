function [kstop ksbot] = secHBstiff(element,pile,y_topbottom,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-v curve.
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kstop   : Soil secant stiffness from the base shear at the top of the element [kN/m/m]
%         ksbot   : Soil secant stiffness from the base shear at the bottom of the element [kN/m/m]
%
% Log: 
% EVVA    28.07.2017      Setting up the function
%--------------------------------------------------------------------------

% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[kstop ksbot] = PISAsandsecHb(element,pile,y_topbottom,i);

elseif strcmp(element.model_py(i),'API sand') || strcmp(element.model_py(i),'Kirsch sand') || strcmp(element.model_py(i),'Kallehave sand')
	[kstop ksbot] = sandsecHb(element,pile,y_topbottom,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [kstop ksbot] = PISAclaysecHb(element,pile,y_topbottom,i);
    
elseif strcmp(element.model_py(i),'API clay') || strcmp(element.model_py(i),'Stiff clay w/o free water') || strcmp(element.model_py(i),'Kirsch soft clay') || strcmp(element.model_py(i),'Kirsch stiff clay') || strcmp(element.model_py(i),'Reese stiff clay')
	[kstop ksbot] = claysecHb(element,pile,y_topbottom,i);
	
elseif strcmp(element.model_py(i),'Zero soil')
    kstop         = 0;
    ksbot         = 0;
    
else
   disp('The specified soil model is not supported in secHBstiff.m and stiffness of toe shear is therefore set to zero')
    kstop           = 0;
    ksbot           = 0;
end
kstop = kstop;
ksbot = ksbot;


 