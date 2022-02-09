function [kttop ktbot] = kirsch2015sandtan(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom of 
% each pile segment by applying p-y curves according to Kirsch (2014).
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
% CODE            : MUOE
% APPROVED        : MMOL  
%
% LAST MODIFIED   : MUOE   23.06.2015   Programming
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
f       = 0.5;                      % Reduction factor for friction angle
phi     = element.phi(i)-f*(D-2);   % Modified triaxial design friction angle [degrees]
Estop   = element.Es(i,1)/1E3;      % Stiffness [MPa] in top          
Esbot   = element.Es(i,1)/1E3;      % and bottom of pile segment - must be in [MPa] for Kirsch equations to work   

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

Estop=k*max(xtop,1)/1000;
Esbot=k*max(ybot,1)/1000;



Esdtop = 10^(-0.42*log10(0.0006*Estop))*Estop; % small strain stiffness for top 
Esdbot = 10^(-0.42*log10(0.0006*Esbot))*Esbot; % and bottom of pile segment

if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	ptop = Atop*putop*tanh(k*xtop*ytop/(0.9*putop)); % starting p-value for for top
	pbot = Abot*pubot*tanh(k*xbot*ybot/(0.9*pubot)); % and bottom of pile segment
else
	ptop = Atop*putop*tanh(k*xtop*ytop/(Atop*putop)); % starting p-value for for top
	pbot = Abot*pubot*tanh(k*xbot*ybot/(Abot*pubot)); % and bottom of pile segment
end

kmodtop = k*(1+(1-ptop/(Atop*putop))*(Esdtop/Estop-1)); % starting value for modified subgrade modulus for top
kmodbot = k*(1+(1-pbot/(Abot*pubot))*(Esdbot/Esbot-1)); % and bottom of pile segment

% -------- Tangent stiffness, kttop and ktbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dp/dy [kN/m^2]. 
% Similar for ktbot

tol = 0.01*k; % tolerance for iteration of k

% iteration to find ktop
if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	ktop = k*(1+(1-tanh(kmodtop*xtop*ytop/(0.9*putop)))*(Esdtop/Estop-1));
	while abs(kmodtop-ktop) > tol  
		kmodtop = (kmodtop+ktop)/2;
		ktop = k*(1+(1-tanh(kmodtop*xtop*ytop/(0.9*putop)))*(Esdtop/Estop-1));
	end
else
	ktop = k*(1+(1-tanh(kmodtop*xtop*ytop/(Atop*putop)))*(Esdtop/Estop-1));
	while abs(kmodtop-ktop) > tol  
		kmodtop = (kmodtop+ktop)/2;
		ktop = k*(1+(1-tanh(kmodtop*xtop*ytop/(Atop*putop)))*(Esdtop/Estop-1));
	end
end

% iteration to find kbot
if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(0.9*pubot)))*(Esdbot/Esbot-1));
	while abs(kmodbot-kbot) > tol
		kmodbot = (kmodbot+kbot)/2;
		kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(Abot*pubot)))*(Esdbot/Esbot-1));
	end
else
	kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(0.9*pubot)))*(Esdbot/Esbot-1));
	while abs(kmodbot-kbot) > tol
		kmodbot = (kmodbot+kbot)/2;
		kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(Abot*pubot)))*(Esdbot/Esbot-1));
	end
end

if xtop == 0;
	% This statement is included because Kttop = dp/dy cannot be evaluated 
	% at xtop = 0 since putop in this case turns zero. However, the stiff-
	% ness at the soil surface is zero
	kttop = 0; 
else
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kttop = Atop/0.9*ktop*xtop/(cosh(ktop*xtop*ytop/(0.9*putop)))^2;
	else
		kttop = ktop*xtop/(cosh(ktop*xtop*ytop/(Atop*putop)))^2;
	end
end
if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	ktbot = Abot/0.9*kbot*xbot/(cosh(kbot*xbot*ybot/(0.9*pubot)))^2;
else
	ktbot = kbot*xbot/(cosh(kbot*xbot*ybot/(Abot*pubot)))^2;
end