function [ktmomtop ktmombot]  = tanmomstiff(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent uniformly distributed moment stiffness [kNm/m] in the top and the bottom of
% each pile segment by applying m-teta curves.
% 
% INPUT:  mu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%						  load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ktmomtop   : Tangent soil moment stiffness at the top of the element [kNm/m/m]
%         ktmombot   : Tanmgent soil moment stiffness at the bottom of the element [kNm/m/m]
%
% Log: 
% EVVA    20.09.2017   Setting up the function
%--------------------------------------------------------------------------
% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[ktmomtop ktmombot] = PISAsandtanm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [ktmomtop ktmombot] = PISAclaytanm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);

elseif strcmp(element.model_py(i),'Zero soil')
    ktmomtop           = 0;
    ktmombot           = 0;
    
else
   disp('The specified soil model is not supported in tanmomstiff.m and the stiffness of the distrubuted moment is set to zero')
    ktmomtop           = 0;
    ktmombot           = 0;  
end
if isfield(element,'degradation_mt_m')
    ktmomtop = ktmomtop*element.degradation_mt_m(i);
    ktmombot = ktmombot*element.degradation_mt_m(i);
else
    ktmomtop = ktmomtop;
    ktmombot = ktmombot;
end


 