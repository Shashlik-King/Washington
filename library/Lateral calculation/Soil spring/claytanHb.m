function [kttop ktbot] = claytanHb(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-y curve according to DNV-OS-J101.
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kttop   : Tangent soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Tangent soil stiffness at the bottom of the element [kN/m/m]
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
ybot1 = abs(y_topbottom(2));            % Horizontal disp. at the bottom of the pile segment [m]
ybot2 = ybot1+0.00001;
xtop = element.heqv(i,1);        		% Distance from the seabed to the bottom of the considered pile segment [m].


pu_cyc      = 0.5*A*su;					% kN
p_ts1        = 2.5*G0*R*ybot1;			% kN

if p_ts1<0
    Hb1     = max(-pu_cyc,p_ts1);
else
    Hb1     = min(pu_cyc,p_ts1);
end

p_ts2        = 2.5*G0*R*ybot2;			% kN

if p_ts2<0
    Hb2     = max(-pu_cyc,p_ts2);
else
    Hb2     = min(pu_cyc,p_ts2);
end

% -------- Change in base shear Hb ---------------------------
Hb = Hb2-Hb1;

% -------- Tangent stiffness, kttop and ktbot -------------------------------

if xtop == 0; %Actually xbotp=0 is not possible, because we're calculating the last pile segment. 
    kttop = 0; 
else
	kttop = Hb./0.00001; %kN/m
end
	kttop=kttop/(element.level(end,1)-element.level(end,2)); %kN/m/m
	ktbot = kttop;