function [kttop ktbot] = tanspringstiff(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom 
% of each pile segment by applying p-y curves.
% 
% INPUT:  springinput   : cf. mainfile
%         putot         : capacity [kN/m] in top and bottom of element
%         heqv          : Equivalent length for top and bottom of element
%         zrtot         : Transition depth [m] - moderate to deep
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kttop   : Soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    08.11.2007   Programming
%                   AHA    15.11.2007   Introduction of sigma.
%                   AHA    21.05.2008   Including layered soil
%                   AHA    01.08.2008   Free standing length
%                   MMOL   04.11.2010   Integrating into calculation routine + stiff clay
%                   JALY   29.04.2011   Streamlining code
%--------------------------------------------------------------------------
% Determination of the model to apply

if strcmp(element.model_py(i),'API sand')
    [kttop ktbot]   = apisandtan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch sand')
    [kttop ktbot]   = kirschsandtan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
	
elseif strcmp(element.model_py(i),'Kirsch sand 2015')
    [kttop ktbot]   = kirsch2015sandtan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kallehave sand')
    [kttop ktbot]   = kallehavesandtan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'API clay')
    [kttop ktbot]   = apiclaytan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch soft clay')
    [kttop ktbot] = kirschsoftclaytan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Reese stiff clay')
    [kttop ktbot] = reesestiffclaytan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Stiff clay w/o free water')
    [kttop ktbot] = stiffclay_w_o_water_tan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch stiff clay')
    [kttop ktbot] = kirschstiffclaytan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Modified Weak rock')
    [kttop ktbot] = weakrocktan(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Sand') 
	[kttop ktbot] = PISAsandtan(element,pile,y_topbottom/element.degradation_py_y(i),i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
	[kttop ktbot] = PISAclaytan(element,pile,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Zero soil')
    kttop           = 0;
    ktbot           = 0;
    
else
   error('The specified soil model is not supported')
   
end
kttop = kttop*element.degradation_py_p(i)/element.degradation_py_y(i);
ktbot = ktbot*element.degradation_py_p(i)/element.degradation_py_y(i);