function [kttop ktbot] = tanHBstiff(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-v curve.
% 
% INPUT:  u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kttop   : Soil tangent stiffness at the top of the element [kN/m/m]
%         ktbot   : Soil tangent stiffness at the bottom of the element [kN/m/m]
%
% Log  : 
% EVVA    11.09.2017   Setting up the code
%--------------------------------------------------------------------------
% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[kttop ktbot] = PISAsandtanHb(element,pile,y_topbottom,i);
	
elseif strcmp(element.model_py(i),'API sand') || strcmp(element.model_py(i),'Kirsch sand') || strcmp(element.model_py(i),'Kallehave sand')
	[kttop ktbot] = sandtanHb(element,pile,y_topbottom,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [kttop ktbot] = PISAclaytanHb(element,pile,y_topbottom,i);

elseif strcmp(element.model_py(i),'API clay') || strcmp(element.model_py(i),'Stiff clay w/o free water') || strcmp(element.model_py(i),'Kirsch soft clay') || strcmp(element.model_py(i),'Kirsch stiff clay') || strcmp(element.model_py(i),'Reese stiff clay')   
	[kttop ktbot] = claytanHb(element,pile,y_topbottom,i);

elseif strcmp(element.model_py(i),'Zero soil')
    kttop           = 0;
    ktbot           = 0;
    
else
   disp('The specified soil model is not supported in tanHBstiff.m and the stiffness of toe shear is set to zero')
    kttop           = 0;
    ktbot           = 0;
   
end
kttop = kttop;
ktbot = ktbot;
