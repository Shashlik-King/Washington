function [kttop ktbot] = apisandtan(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom
% of each pile segment by applying p-y curves according to API(1993).
% 
% INPUT:  springinput   : cf. mainfile [D phi gamma 1]
%                         1 = static behaviour
%                         2 = cyclic behaviour
%         putot         : capacity [kN/m] in top and bottom of element
%         heqv          : Equivalent length for top and bottom of element
%         u             : Global displacment vector
%         i             : Counter referring to element number
%
% OUTPUT: kttop   : Soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA  

% LAST MODIFIED   : AHA    11.07.2007   Programming
%                   AHA    18.07.2007   Correcting - ytop,ybot
%                   AHA    22.08.2007   Lines 116 and 123 - the way the 
%                                       information is printed on screen is
%                                       changed.
%                   AHA    08.11.2007   Name of function is changed from
%                                       tanspringstiff to apisandtan
%                   AHA    15.11.2007   Introduction of sigma.
%                   AHA    20.05.2008   Introduction of layered soil.
%                   MMOL   04.11.2010   Integrating into calculation routine
%                   JALY   29.04.2011   Streamlining code
%--------------------------------------------------------------------------

% -------- Initializing pile and soil parameters --------------------------

ytop    = y_topbottom(1);           % Horizontal disp. at the top of the
                                    % pile segment [m]

ybot    = y_topbottom(2);           % Horizontal disp. at the bottom of the
                                    % pile segment [m]

xtop    = element.heqv(i,1);        % Distance from pile head to the top
                                    % of the considered pile segment [m].
                                    % Layered soil has been taken into
                                    % account.

xbot    = element.heqv(i,2);        % Distance from pile head to the bottom
                                    % of the considered pile segment [m].
                                    % Layered soil has been taken into
                                    % account.

D       = pile.diameter;            % Outer diameter [m]
phi     = element.phi(i);           % Triaxial design friction angle [degrees]

putop   = element.pu(i,1);          % Ultimate resistance [kN/m] in top
pubot   = element.pu(i,2);          % and bottom of pile segment, determined in layer.m                   

% -------- Determination of A ---------------------------------------------

element = A_factor(element,pile,loads,i);
Atop    = element.A(i,1);
Abot    = element.A(i,2);

% -------- Initial modulus of subgrade reaction ---------------------------

subgrade = 2;

if subgrade == 1
    % exponential fit
    k = 158.89*exp(0.1411*phi); % static subgrade modulus
elseif subgrade == 2
    % polynomial fit
    if phi <= 28
        k = 5400;
    else
        k = (0.1978*phi^2-10.232*phi+136.82)*1000;
    end
end

% -------- Tangent stiffness, kttop and ktbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dp/dy [kN/m^2]. 
% Similar for ktbot
% Expressions based on API(1993) p.65

if xtop == 0;
    % This statement is included because Kttop = dp/dy cannot be evaluated 
    % at xtop = 0 since putop in this case turns zero. However, the stiff-
    % ness at the soil surface is zero
    kttop = 0; 
else
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kttop = Atop/0.9*k*xtop/(cosh(k*xtop*ytop/(0.9*putop)))^2;
	else
		kttop = k*xtop/(cosh(k*xtop*ytop/(Atop*putop)))^2;
	end
end
if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	ktbot = Abot/0.9*k*xbot/(cosh(k*xbot*ybot/(0.9*pubot)))^2;
else
	ktbot = k*xbot/(cosh(k*xbot*ybot/(Abot*pubot)))^2;
end