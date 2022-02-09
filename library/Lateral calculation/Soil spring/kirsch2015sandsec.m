function [kstop ksbot] = kirsch2015sandsec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to Kirsch (2014).
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
	

% -------- Secant stiffness, kstop and ksbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness p/y [kN/m^2] for y > 0. For y = 0
% kstop is the tangent stiffness, i.e. kstop = dp/dy evaluated at y = 0. 
% Similar for kbot

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

% determine kstop from iterated ktop
if ytop == 0;
	% Expression based on API(1993) p.65, kstop = dp/dy at y = 0
	kstop = ktop*xtop; 
elseif xtop == 0;
	kstop = 0;
else
	% Expression based on API(1993) p.65    
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kstop = Atop*putop/ytop*tanh(ktop*xtop*ytop/(0.9*putop));
	else
		kstop = Atop*putop/ytop*tanh(ktop*xtop*ytop/(Atop*putop));
	end
end
if isnan(kstop)
	kstop = 0;
% 	disp(['kstop = NaN is substituted by kstop = 0 for element no. ',num2str(i),' in kirschsandsec.m due to putop = 0'])
end

% iteration to find kbot
if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(0.9*pubot)))*(Esdbot/Esbot-1));
	while abs(kmodbot-kbot) > tol
		kmodbot = (kmodbot+kbot)/2;
		kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(0.9*pubot)))*(Esdbot/Esbot-1));
	end
else
	kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(Abot*pubot)))*(Esdbot/Esbot-1));
	while abs(kmodbot-kbot) > tol
		kmodbot = (kmodbot+kbot)/2;
		kbot = k*(1+(1-tanh(kmodbot*xbot*ybot/(Abot*pubot)))*(Esdbot/Esbot-1));
	end
end

% determine ksbot from iterated kbot
if ybot == 0;
	ksbot = kbot*xbot;
else
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		ksbot = Abot*pubot/ybot*tanh(kbot*xbot*ybot/(0.9*pubot));
	else
		ksbot = Abot*pubot/ybot*tanh(kbot*xbot*ybot/(Abot*pubot));
	end
end