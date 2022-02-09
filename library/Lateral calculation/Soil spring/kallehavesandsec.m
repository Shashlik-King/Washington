function [kstop ksbot] = kallehavesandsec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to Kallehave (2012).
% 
% INPUT:  springinput   : cf. mainfile [D phi gamma 1]
%                         1 = static behaviour
%                         2 = cyclic behaviour
%         putot         : capacity [kN/m] in top and bottom of element
%         heqv          : Equivalent length for top and bottom of element
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : MUOE
% APPROVED        : MMOL

% LAST MODIFIED   : MUOE    16.09.2015   Programming
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

% Parameters for modifying stiffness

D0 = 0.61; % [m]
x0 = 2.5; % [m]
n = 0.5; % site specific parameter in the range 0.4-0.7 acc. to Kallehave (2012)
m = 0.5; % parameter is a function of y. 0.5 is a good approximate value acc. to Kallehave (2012)

k_mod_top = k*(xtop/x0)^n*(D/D0)^m;
k_mod_bot = k*(xbot/x0)^n*(D/D0)^m;

% -------- Secant stiffness, kstop and ksbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness p/y [kN/m^2] for y >0. For y = 0
% kstop is the tangent stiffness, i.e. kstop = dp/dy evaluated at y = 0. 
% Similar for kbot

if ytop == 0;
    % Expression based on API(1993) p.65, kstop = dp/dy at y = 0
    kstop = k_mod_top*x0; 
elseif xtop == 0;
    kstop = 0;
else
    % Expression based on API(1993) p.65
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kstop = Atop*putop/ytop*tanh(k_mod_top*x0*ytop/(0.9*putop));
	else
		kstop = Atop*putop/ytop*tanh(k_mod_top*x0*ytop/(Atop*putop));
	end
end

if ybot == 0;
    ksbot = k_mod_bot*x0;
else
    if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		ksbot = Abot*pubot/ybot*tanh(k_mod_bot*x0*ybot/(0.9*pubot));
	else
		ksbot = Abot*pubot/ybot*tanh(k_mod_bot*x0*ybot/(Abot*pubot));
	end
end