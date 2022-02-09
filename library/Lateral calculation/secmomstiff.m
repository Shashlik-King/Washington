function [ksmomtop ksmombot]  = secmomstiff(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant uniformly distributed moment stiffness [kNm/m] in the top and the bottom of
% each pile segment by applying m-teta curves.
% 
% INPUT:  mu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%                         load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ksmomtop   : Secant soil moment stiffness at the top of the element [kNm/m/m]
%         ksmombot   : Secant soil moment stiffness at the bottom of the element [kNm/m/m]
%
% Log:
% EVVA    17.08.2016   Setting up the function
%--------------------------------------------------------------------------
% Determination of the model to apply

if element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[ksmomtop ksmombot] = PISAsandsecm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
    [ksmomtop ksmombot] = PISAclaysecm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);    

elseif strcmp(element.model_py(i),'Zero soil')
    ksmomtop         = 0;
    ksmombot         = 0;
    
else
   disp('The specified soil model is not supported in secmomstiff.m and the stiffness of the distributed moment is set to zero')
    ksmomtop           = 0;
    ksmombot           = 0;  
end
if isfield(element,'degradation_mt_m')
    ksmomtop = ksmomtop*element.degradation_mt_m(i);
    ksmombot = ksmombot*element.degradation_mt_m(i);
else
    ksmomtop = ksmomtop;
    ksmombot = ksmombot;
end

 