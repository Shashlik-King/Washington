function [kstop ksbot] = apisandsec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to Jeanjean (2009).
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
% CODE            : MDGI
% APPROVED        : 

% LAST MODIFIED   : MDGI    15.03.2022   Programming
% ATTENTION!!!! factor of 0.8 for cyclic and fs =100 is hard coded here !!!
%--------------------------------------------------------------------------

% -------- Initializing pile and soil parameters --------------------------

ytop    = y_topbottom(1);            % Horizontal disp. at the top of the
                                    % pile segment [m]

ybot    = y_topbottom(2);            % Horizontal disp. at the bottom of the
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
% phi     = element.phi(i);           % Triaxial design friction angle [degrees]

putop   = element.pu(i,1);          % Ultimate resistance [kN/m] in top
pubot   = element.pu(i,2);          % and bottom of pile segment, determined in layer.m

G0top   = element.G0(i,1);          % Small Strain Shear Modulus [kPa] in top
G0bot   = element.G0(i,2);          % Small Strain Shear Modulus in bottom

cutop   = element.cu(i,1);          % Small Strain Shear Modulus [kPa] in top
cubot   = element.cu(i,2);          % Small Strain Shear Modulus in bottom

fs=100;                             % value fo 100 was taken from Jeanjean 2009
k_limit = 10;                    % without this limit there might be convergancy problem

% -------- Determination of A ---------------------------------------------

element = A_factor(element,pile,loads,i);
Atop    = element.A(i,1);
Abot    = element.A(i,2);

% -------- Initial modulus of subgrade reaction ---------------------------

% subgrade = 2;
% 
% if subgrade == 1
%     % exponential fit
%     k = 158.89*exp(0.1411*phi); % static subgrade modulus
% elseif subgrade == 2
%     % polynomial fit
%     if phi <= 28
%         k = 5400;
%     else
%         k = (0.1978*phi^2-10.232*phi+136.82)*1000;
%     end
% end

% -------- Secant stiffness, kstop and ksbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness p/y [kN/m^2] for y >0. For y = 0
% kstop is the tangent stiffness, i.e. kstop = dp/dy evaluated at y = 0. 
% Similar for kbot

if ytop == 0
    % to aviod problem with calculation at y = 0
    kstop = k_limit*G0top; 
elseif xtop == 0
    kstop = 0;
else
    % Expression based on 
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kstop = 0.8*putop/ytop*tanh( sqrt(ytop/D) * G0top/(fs*cutop));
	else
		kstop = putop/ytop*tanh( sqrt(ytop/D) * G0top/(fs*cutop));
    end
    
    %mkae sure siffness has a cap
    if kstop>(k_limit*G0top)
        kstop = k_limit*G0top;
    end    
end

if ybot ==0
    % to aviod problem with calculation at y = 0
    ksbot = k_limit*G0bot;
else
    if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		ksbot = 0.8*pubot/ybot*tanh( sqrt(ybot/D) * G0bot/(fs*cubot));
	else
		ksbot = pubot/ybot*tanh( sqrt(ybot/D) * G0bot/(fs*cubot));
    end
    %mkae sure siffness has a cap
    if ksbot>(k_limit*G0bot)
        ksbot = k_limit*G0bot;
    end
end