function [ktmomtop ktmombot] = tanMBstiff(element,pile,teta_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent base moment stiffness [kN/m/m] from the base moment in the top and the bottom of
% the last pile segment by applying Mb-teta curve.
% 
% INPUT:  Mbu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%						  load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ktmomtop    : Soil tangent stiffness at the top of the element [kN/m/m]
%         ktmombot   : Soil tangent stiffness at the bottom of the element [kN/m/m]
%
% Log  : 
% EVVA    11.09.2017   Setting up the code
%--------------------------------------------------------------------------
% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')  
	[ktmomtop ktmombot] = PISAsandtanMb(element,pile,teta_topbottom,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [ktmomtop ktmombot] = PISAclaytanMb(element,pile,teta_topbottom,i);

elseif strcmp(element.model_py(i),'Zero soil')
    ktmomtop           = 0;
    ktmombot           = 0;
    
else
   disp('The specified soil model is not supported in tanMBstiff.m and the stiffness of base moment is set to zero')
    ktmomtop           = 0;
    ktmombot           = 0;  
end
ktmomtop = ktmomtop;
ktmombot = ktmombot;
