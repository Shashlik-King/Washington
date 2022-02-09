function [ksmomtop ksmombot] = PISAclaysecMb(element,pile,teta_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant stiffness [kN/m/rad] from the base moment in the top and the bottom of
% the last pile segment by applying Mb-teta curve according to PISA sand model. 
% 
% INPUT:  
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: ksmomtop   : Soil stiffness at the top of the (last) element [kN/m/rad]
%         ksmombot   : Soil stiffness at the bottom of the (last) element [kN/m/rad]
%
% Log: 
% EVVA    18.09.2017  Programming
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------
D    = pile.diameter;
L   = abs(element.level(end,1));	   % pile.length_end;
%G_top   = element.G0(i,1);             % Shear modulus at the top of the element [kPa]
G0   = element.G0(i,2);             % Shear modulus at the bottom of the element [kPa]
Mbu_n   = element.MBu(i,2);            % Normalized base shear ultimate resistance [-]
su  = element.cu(i,2);   		% Undrained shear strength in the bottom node[kPa]
tetabot = abs(teta_topbottom(2));                  % Horizontal disp. at the bottom of the pile segment [m]
xbot = element.heqv(i,2);        % Distance from the seabed to the bottom of the considered pile segment [m].


% -------- Initial stiffness moduli -------------------------------------------------
kM_n = element.PISA_param(i,21);  % Normalized initial stiffness for base moment

kMb = kM_n*D^3*G0;       % Non-normalized initial stiffness for base moment

% -------- Curvature ------------------------------------------------------
nM = element.PISA_param(i,22);

% -------- Normalized ultimate nodal displacement and rotation -----------------------------
teta_u_n  = element.PISA_param(i,19); % Same for Cowden clay and bothkennar clay

% --------- Normalized nodal displacement and rotation -------------------------------------
teta_n = tetabot*(G0/su);

% -------- Calculation of base moment Mb ---------------------------
aM_n = 1-2*nM;
bM_n = 2*nM*(teta_n/teta_u_n)-(1-nM)*(1+(teta_n*kM_n/Mbu_n));
cM_n = (teta_n*kM_n/Mbu_n)*(1-nM)-nM*(teta_n/teta_u_n).^2;

if teta_n <= teta_u_n
Mb_n = Mbu_n.*((2*cM_n)./(-bM_n+sqrt(bM_n.^2-4*aM_n*cM_n))); %Normalized pile resistance 
Mb = Mb_n.*su*D^3; %Non-normalized pile resistance 
else
Mb_n= Mbu_n; %Normalized ultimate pile resistance
Mb = Mb_n.*su*D^3; %Non-normalized ultimate pile resistance 
end

% -------- Secant stiffness, ksmomtop and ksmombot -------------------------------
% Determination of base moment Mb spring stiffness at the top and bottom of the last pile
% segment, ksmomtop is the secant stiffness Mb/teta [kNm/m/rad] for teta >0. For teta = 0
% ksmomtop is the tangent stiffness, i.e. ksmomtop = dMb/dteta evaluated at teta = 0. 

if tetabot == 0;
    ksmombot = kMb.*xbot; 
else
	ksmombot = Mb./tetabot; %kNm/rad
end
ksmombot=ksmombot/(element.level(end,1)-element.level(end,2)); %kNm/rad/m
ksmomtop = ksmombot;

