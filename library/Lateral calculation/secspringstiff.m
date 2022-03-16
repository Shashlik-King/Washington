function [kstop ksbot] = secspringstiff(element,pile,loads,y_topbottom,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves.
% 
% INPUT:  springinput   : cf. mainfile
%         putot         : capacity [kN/m] in top and bottom of element
%         heqv          : Equivalent length for top and bottom of element
%         zrtot         : Transition depth [m] - moderate to deep
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    05.11.2007   Programming
%                   AHA    15.11.2007   Including sigma
%                   AHA    21.05.2008   Including layered soil
%                   AHA    01.08.2008   Free standing length
%                   MMOL   04.11.2010   Integrating into calculation routine + stiff clay
%                   JALY   29.04.2011   Streamlining code
%--------------------------------------------------------------------------

% Determination of the model to apply

if strcmp(element.model_py(i),'API sand')
    [kstop ksbot] = apisandsec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch sand')
    [kstop ksbot] = kirschsandsec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
	
elseif strcmp(element.model_py(i),'Kirsch sand 2015')
    [kstop ksbot] = kirsch2015sandsec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);

    
elseif strcmp(element.model_py(i),'Kallehave sand')
    [kstop ksbot] = kallehavesandsec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'API clay')
    [kstop ksbot] = apiclaysec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch soft clay')
    [kstop ksbot] = kirschsoftclaysec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);

elseif strcmp(element.model_py(i),'Reese stiff clay')
    [kstop ksbot] = reesestiffclaysec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Stiff clay w/o free water')
    [kstop ksbot] = stiffclay_w_o_water_sec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Jeanjean clay')
    [kstop ksbot] = jeanjeanclaysec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Kirsch stiff clay')
    [kstop ksbot] = kirschstiffclaysec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i);
    
elseif strcmp(element.model_py(i),'Modified Weak rock')
    [kstop ksbot] = weakrocksec(element,pile,loads,y_topbottom/element.degradation_py_y(i),i); 
  
elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Sand')
	[kstop ksbot] = PISAsandsec(element,pile,y_topbottom/element.degradation_py_y(i),i);

elseif element.PISA_switch==1 && strcmp(element.type{i}, 'Clay')
	[kstop ksbot] = PISAclaysec(element,pile,y_topbottom/element.degradation_py_y(i),i); 
    
elseif strcmp(element.model_py(i),'Zero soil')
    kstop         = 0;
    ksbot         = 0;
    
else
    error('The specified soil model is not supported')
    
end
kstop = kstop*element.degradation_py_p(i)/element.degradation_py_y(i);
ksbot = ksbot*element.degradation_py_p(i)/element.degradation_py_y(i);