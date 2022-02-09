function [kstop ksbot] = PISAsandsecHb(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] from the base shear in the top and the bottom of
% the last pile segment by applying Hb-y curve according to PISA sand model. 
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% Log: 
% EVVA    11.09.2017  Programming
% FKMV    17.08.2020  Update to use database
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------
D    = pile.diameter;
L   = abs(element.level(end,1));	   % pile.length_end;
%G_top   = element.G0(i,1);             % Shear modulus at the top of the element [kPa]
G0   = element.G0(i,2);                 % Shear modulus at the bottom of the element [kPa]
Hbu_n   = element.HBu(i,2);            	% Normalized base shear ultimate resistance [-]
sigma_v = element.sigma_v_eff (i,2);   	% Effective vertical stress at the base of the pile [kPa]
ybot = abs(y_topbottom(2));            	% Horizontal disp. at the bottom of the pile segment [m]
xbot = element.heqv(i,2);        		% Distance from the seabed to the bottom of the considered pile segment [m].

% % This is only to replicate PISA results
% B = 875; % Fitted parameter in PISA report [-]
% e = 0.628; % Void ratio  [-]
% K0 = 0.4; % Lateral eath pressure coefficient [-]
% p_ref = 101.3; % Reference pressure in PISA report [kPa]
% p = (sigma_v+2*K0*sigma_v)/3; % Effective mean stress at the element top [kPa]

% G0 = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p/p_ref); % PISA report pg. 83/269

% -------- Initial stiffness moduli -------------------------------------------------
kH_n = element.PISA_param(i,17);   % Normalized initial stiffness for base shear
% kH_n = -0.38*L/D+3.02;   % Normalized initial stiffness for base shear - According to PISA 2 valid for Dr=75%
% kH_n = -0.37*L/D+3.07;   % Normalized initial stiffness for base shear - According to PISA 2 valid for Dr=45%
% kH_n = -0.37*L/D+3.05;   % Normalized initial stiffness for base shear - According to PISA 2 valid for Dr=60%
% kH_n = -0.38*L/D+3.02;   % Normalized initial stiffness for base shear - According to PISA 2 valid for Dr=90%
kH = kH_n*D*G0;          % Non-normalized initial stiffness for base shear

% -------- Curvature ------------------------------------------------------
nH = element.PISA_param(i,18);
% nH = -0.05*L/D+0.94; % - According to PISA 2 valid for Dr=75%
% nH = -0.04*L/D+0.90; % - According to PISA 2 valid for Dr=45%
% nH = -0.06*L/D+0.94; % - According to PISA 2 valid for Dr=60%
% nH = -0.06*L/D+0.95; % - According to PISA 2 valid for Dr=90%

% -------- Normalized ultimate nodal displacement and rotation -----------------------------
v_u_n  = element.PISA_param(i,15); % Same for Cowden clay and Bothkennar clay
% v_u_n  = -0.29*L/D+2.31; % - According to PISA 2 valid for Dr=75%
% v_u_n  = -0.20*L/D+2.17; % - According to PISA 2 valid for Dr=45%
% v_u_n  = -0.20*L/D+2.07; % - According to PISA 2 valid for Dr=60%
% v_u_n  = -0.48*L/D+3.33; % - According to PISA 2 valid for Dr=90%

% --------- Normalized nodal displacement and rotation -------------------------------------
v_n = ybot/D*(G0/sigma_v);


% -------- Calculation of base shear Hb ---------------------------
aH_n = 1-2*nH;
bH_n = 2*nH*(v_n/v_u_n)-(1-nH)*(1+(v_n*kH_n/Hbu_n));
cH_n = (v_n*kH_n/Hbu_n)*(1-nH)-nH*(v_n/v_u_n).^2;

if v_n <= v_u_n
Hb_n = Hbu_n.*((2*cH_n)./(-bH_n+sqrt(bH_n.^2-4*aH_n*cH_n))); %Normalized pile resistance 
Hb = Hb_n.*sigma_v*D^2; %Non-normalized pile resistance 
else
Hb_n = Hbu_n; %Normalized ultimate pile resistance
Hb= Hb_n.*sigma_v*D^2; %Non-normalized ultimate pile resistance 
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
ksbot=ksbot/(element.level(end,1)-element.level(end,2));  %kN/m^2
kstop = ksbot;