function [ksmomtop ksmombot] = secMBstiff(element,pile,teta_topbottom,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] from the base moment in the top and the bottom of
% the last pile segment by applying Mb-v curve.
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ksmomtop   : Soil secant stiffness from the base shear at the top of the element [kN/m/m]
%         ksmombot   : Soil secant stiffness from the base shear at the bottom of the element [kN/m/m]
%
% Log: 
% EVVA    28.07.2017      Setting up the function
%--------------------------------------------------------------------------

% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[ksmomtop ksmombot] = PISAsandsecMb(element,pile,teta_topbottom,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [ksmomtop ksmombot] = PISAclaysecMb(element,pile,teta_topbottom,i);
  
elseif strcmp(element.model_py(i),'Zero soil')
    ksmomtop         = 0;
    ksmombot         = 0;
    
else
   disp('The specified soil model is not supported in secMBstiff.m and stiffness of base moment is therefore set to zero')
    ksmomtop           = 0;
    ksmombot           = 0;  
end
ksmomtop = ksmomtop;
ksmombot = ksmombot;


 