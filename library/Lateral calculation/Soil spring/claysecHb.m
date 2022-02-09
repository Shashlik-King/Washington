function [kstop ksbot] = claysecHb(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-y curve according to DNV-OS-J101. 
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% Log: 
% EVVA    25.09.2017  Programming
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------
D  = pile.diameter;
A  = pi*(D/2)^2; 						% Cross sectional area. Used for calculating the toe shear
R  = D/2; 								% Radius of foundation [m]
su = element.cu(i,2);
OCR = 2;
sigma_toe = element.sigma_v_eff (i,2);  % Effective vertical stress at the base of the pile [kPa]
%G0 = element.G0(i,2);                   % Shear modulus at the bottom of the element [kPa]
G0 = 600*su-170*su*sqrt(OCR-1);   
ybot = abs(y_topbottom(2));            % Horizontal disp. at the bottom of the pile segment [m]


pu_cyc      = 0.5*A*su;					% kN
p_ts        = 2.5*G0*R*ybot;			% kN


if p_ts<0
    Hb     = max(-pu_cyc,p_ts);			 % kN
else
    Hb     = min(pu_cyc,p_ts);			 % kN
end

% -------- Secant stiffness, kstop and ksbot -------------------------------
% Determination of base shear Hb spring stiffness at the top and bottom of the last pile
% segment, ks is the secant stiffness Hb/y [kN/m^2] for y >0. For y = 0 Hb=0, ks=0.

if ybot == 0;
    ksbot = 0;
else
	ksbot = Hb./ybot; %kN/m
end
ksbot=ksbot/(element.level(end,1)-element.level(end,2));  %kN/m/m
kstop = ksbot;