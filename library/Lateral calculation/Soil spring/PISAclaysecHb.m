function [kstop ksbot] = PISAclaysecHb(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-y curve according to PISA clay model. 
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% Log: 
% EVVA    12.09.2017  Programming
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------
D    = pile.diameter;
L   = abs(element.level(end,1));	   % pile.length_end;
%G_top   = element.G0(i,1);             % Shear modulus at the top of the element [kPa]
G0   = element.G0(i,2);             % Shear modulus at the bottom of the element [kPa]
Hbu_n   = element.HBu(i,2);            % Normalized base shear ultimate resistance [-]
su  = element.cu(i,2);   		% Undrained shear strength in the bottom node[kPa]
ybot = abs(y_topbottom(2));            % Horizontal disp. at the bottom of the pile segment [m]
xbot = element.heqv(i,2);        % Distance from the seabed to the bottom of the considered pile segment [m].

% -------- Initial stiffness moduli -------------------------------------------------
kH_n = element.PISA_param(i,17);   % Normalized initial stiffness for base shear

kH = kH_n*D*G0;          % Non-normalized initial stiffness for base shear

% -------- Curvature ------------------------------------------------------
nH = element.PISA_param(i,18);

% -------- Normalized ultimate nodal displacement and rotation -----------------------------
v_u_n  = element.PISA_param(i,15); % Same for Cowden clay and Bothkennar clay

% --------- Normalized nodal displacement and rotation -------------------------------------
v_n = ybot/D*(G0/su);


% -------- Calculation of base shear Hb ---------------------------
aH_n = 1-2*nH;
bH_n = 2*nH*(v_n/v_u_n)-(1-nH)*(1+(v_n*kH_n/Hbu_n));
cH_n = (v_n*kH_n/Hbu_n)*(1-nH)-nH*(v_n/v_u_n).^2;

if v_n <= v_u_n
Hb_n = Hbu_n.*((2*cH_n)./(-bH_n+sqrt(bH_n.^2-4*aH_n*cH_n))); %Normalized pile resistance 
Hb = Hb_n.*su*D^2; %Non-normalized pile resistance 
else
Hb_n = Hbu_n; %Normalized ultimate pile resistance
Hb= Hb_n.*su*D^2; %Non-normalized ultimate pile resistance 
end


% -------- Secant stiffness, kstop and ksbot -------------------------------
% Determination of base shear Hb spring stiffness at the top and bottom of the last pile
% segment, kstop is the secant stiffness Hb/y [kN/m^2] for y >0. For y = 0
% kstop is the tangent stiffness, i.e. kstop = dHb/dy evaluated at y = 0. 

if ybot == 0;
    ksbot = kH.*xbot;
else
	ksbot = Hb./ybot; %kN/m
end
ksbot=ksbot/(element.level(end,1)-element.level(end,2)); %kN/m^2
kstop = ksbot;