function [kttop ktbot] = PISAclaytanHb(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-y curve according to PISA clay model. 
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kttop   : Tangent soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Tangent soil stiffness at the bottom of the element [kN/m/m]
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
xtop = element.heqv(i,1);        % Distance from the seabed to the bottom of the considered pile segment [m].

% -------- Initial stiffness moduli -------------------------------------------------
kH_n = element.PISA_param(i,17);   % Normalized initial stiffness for base shear

kH = kH_n*D*G0;          % Non-normalized initial stiffness for base shear

% -------- Curvature ------------------------------------------------------
nH = element.PISA_param(i,18);

% -------- Normalized ultimate nodal displacement and rotation -----------------------------
v_u_n  = element.PISA_param(i,15); % Same for Cowden clay and Bothkennar clay

% --------- Normalized nodal displacement and rotation -------------------------------------
v_n1 = ybot/D*(G0/su);
v_n2 = (ybot+0.00001)/D*(G0/su);

% -------- Calculation of base shear Hb ---------------------------
aH_n = 1-2*nH;

% Parameter b1
bH_n1 = 2*nH*(v_n1/v_u_n)-(1-nH)*(1+(v_n1*kH_n/Hbu_n));
% Parameter b2
bH_n2 = 2*nH*(v_n2/v_u_n)-(1-nH)*(1+(v_n2*kH_n/Hbu_n));
% Parameter c1
cH_n1 = (v_n1*kH_n/Hbu_n)*(1-nH)-nH*(v_n1/v_u_n).^2;
% Parameter c2
cH_n2 = (v_n2*kH_n/Hbu_n)*(1-nH)-nH*(v_n2/v_u_n).^2;


if v_n1 <= v_u_n
Hb_n1 = Hbu_n.*((2*cH_n1)./(-bH_n1+sqrt(bH_n1.^2-4*aH_n*cH_n1))); %Normalized pile resistance 
Hb1 = Hb_n1.*su*D^2; %Non-normalized pile resistance 
else
Hb_n1 = Hbu_n; %Normalized ultimate pile resistance
Hb1= Hb_n1.*su*D^2; %Non-normalized ultimate pile resistance 
end

if v_n2 <= v_u_n
Hb_n2 = Hbu_n.*((2*cH_n2)./(-bH_n2+sqrt(bH_n2.^2-4*aH_n*cH_n2))); %Normalized pile resistance 
Hb2 = Hb_n2.*su*D^2; %Non-normalized pile resistance 
else
Hb_n2 = Hbu_n; %Normalized ultimate pile resistance
Hb2= Hb_n2.*su*D^2; %Non-normalized ultimate pile resistance 
end

Hb = Hb2-Hb1;

% -------- Tangent stiffness, kttop and ktbot -------------------------------

if xtop == 0; %Actually xtop=0 is not possible, because we're calculating the last pile segment. 
    kttop = 0; 
else
	kttop = Hb./0.00001; %kN/m
end
	kttop=kttop/(element.level(end,1)-element.level(end,2)); %kN/m^2
	ktbot = kttop;


% if xtop == 0;
%     % This statement is included because Kttop = dp/dy cannot be evaluated 
%     % at xtop = 0 since putop_n in this case turns zero. However, the stiff-
%     % ness at the soil surface is zero/low
%     kttop = 0; 
% else
% 	kttop = Hb./0.00001;
% end
% 	kttop=kttop/(element.level(end,1)-element.level(end,2));%SPSO added
% 	ktbot = Hb./0.00001;
% 	ktbot=ktbot/(element.level(end,1)-element.level(end,2));%SPSO added