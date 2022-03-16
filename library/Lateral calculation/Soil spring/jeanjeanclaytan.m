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
% CODE            : MDGI
% APPROVED        : 

% LAST MODIFIED   : MDGI    15.03.2022   Programming

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

G0top   = element.G0(i,1);          % Small Strain Shear Modulus [kPa] in top
G0bot   = element.G0(i,2);          % Small Strain Shear Modulus in bottom

cutop   = element.cu(i,1);          % Small Strain Shear Modulus [kPa] in top
cubot   = element.cu(i,2);          % Small Strain Shear Modulus in bottom

fs = 100;% a factor to be calibrated from FEM, but here it is assumed.
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

% -------- Tangent stiffness, kttop and ktbot ------------------------------

% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dp/dy [kN/m^2]. 
% Similar for ktbot
% Expressions based on API(1993) p.65

if ytop == 0;
    % to aviod problem with calculation at y = 0
    ytop = 1e-5; 
end


if xtop == 0;
    % This statement is included because Kttop = dp/dy cannot be evaluated 
    % at xtop = 0 since putop in this case turns zero. However, the stiff-
    % ness at the soil surface is zero
    kttop = 0; 
else
	if strcmp(loads.static_cyclic,'cyclic') % cyclic case
		kttop = 0.8*( putop*G0top/(2*fs*cutop*sqrt(D*ytop))  )  /(cosh(sqrt(ytop/D) * G0top/(fs*cutop) ))^2;
	else
		kttop = 0.8*( putop*G0top/(2*fs*cutop*sqrt(D*ytop))  )  /(cosh(sqrt(ytop/D) * G0top/(fs*cutop) ))^2;
	end
end

if ybot == 0;
    % to aviod problem with calculation at y = 0
    ybot = 1e-5; 
end

if strcmp(loads.static_cyclic,'cyclic') % cyclic case
	ktbot = 0.8*( pubot*G0bot/(2*fs*cubot*sqrt(D*ybot))  )  /(cosh(sqrt(ybot/D) * G0bot/(fs*cubot) ))^2;
else
	ktbot = ( pubot*G0bot/(2*fs*cubot*sqrt(D*ybot))  )  /(cosh(sqrt(ybot/D) * G0bot/(fs*cubot) ))^2;
end