function  [ktmomtop ktmombot] = PISAclaytanMb(element,pile,teta_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent stiffness [kN/m/rad] from the base moment in the top and the bottom of
% the last pile segment by applying Mb-teta curve according to PISA clay model. 
%
% INPUT:  Mbu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%						  load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ktmomtop   : Soil stiffness at the top of the (last) element [kN/m/rad]
%         ktmombot   : Soil stiffness at the bottom of the (last) element [kN/m/rad]
%
% Log: 
% EVVA   12.09.2017   Writting the code
%--------------------------------------------------------------------------
D   = pile.diameter;
L   = abs(element.level(end,1));	   % pile.length_end;
G0       = element.G0(i,2);            % Shear modulus at the base of the pile [kPa]
Mbu_n   = element.MBu(i,2);            % Normalized base shear ultimate resistance [-]
su   = element.cu(i,2);              % Undrained shear strength in the bottom node[kPa]
tetabot = abs(teta_topbottom(2));                  % Horizontal disp. at the bottom of the pile segment [m]
xtop = element.heqv(i,1);        % Distance from the seabed to the bottom of the considered pile segment [m].

% -------- Initial stiffness moduli -------------------------------------------------
kM_n = element.PISA_param(i,21);  % Normalized initial stiffness for base moment

kMb = kM_n*D^3*G0;       % Non-normalized initial stiffness for base moment

% -------- Curvature ------------------------------------------------------
nM = element.PISA_param(i,22);

% -------- Normalized ultimate nodal displacement and rotation -----------------------------
teta_u_n  = element.PISA_param(i,19); % Same for Cowden clay and bothkennar clay

% --------- Normalized nodal displacement and rotation -------------------------------------
teta_n1 = tetabot*(G0/su);
teta_n2 = (tetabot+0.00001)*(G0/su);

% -------- Calculation of base moment Mb ---------------------------
aM_n = 1-2*nM;
% Parameter b1
bM_n1 = 2*nM*(teta_n1/teta_u_n)-(1-nM)*(1+(teta_n1*kM_n/Mbu_n));
% Parameter b2
bM_n2 = 2*nM*(teta_n2/teta_u_n)-(1-nM)*(1+(teta_n2*kM_n/Mbu_n));
% Parameter c1
cM_n1 = (teta_n1*kM_n/Mbu_n)*(1-nM)-nM*(teta_n1/teta_u_n).^2;
% Parameter c2
cM_n2 = (teta_n2*kM_n/Mbu_n)*(1-nM)-nM*(teta_n2/teta_u_n).^2;


if teta_n1 <= teta_u_n
Mb_n1 = Mbu_n.*((2*cM_n1)./(-bM_n1+sqrt(bM_n1.^2-4*aM_n*cM_n1))); %Normalized pile resistance 
Mb1 = Mb_n1.*su*D^3; %Non-normalized pile resistance 
else
Mb_n1= Mbu_n; %Normalized ultimate pile resistance
Mb1 = Mb_n1.*su*D^3; %Non-normalized ultimate pile resistance 
end

if teta_n2 <= teta_u_n
Mb_n2 = Mbu_n.*((2*cM_n2)./(-bM_n2+sqrt(bM_n2.^2-4*aM_n*cM_n2))); %Normalized pile resistance 
Mb2 = Mb_n2.*su*D^3; %Non-normalized pile resistance 
else
Mb_n2= Mbu_n; %Normalized ultimate pile resistance
Mb2 = Mb_n2.*su*D^3; %Non-normalized ultimate pile resistance 
end

Mb = Mb2-Mb1;

% -------- Tangent stiffness, ktmomtop and ktmombot -------------------------------

if xtop == 0; %Actually xtop=0 is not possible, because we're calculating the last pile segment. 
    ktmomtop = 0; 
else
	ktmomtop = Mb./0.00001; %kNm/rad
end
	ktmomtop=ktmomtop/(element.level(end,1)-element.level(end,2)); %kN/rad/m
	ktmombot = ktmomtop;

% if xtop == 0;
%     % This statement is included because Kttop = dp/dy cannot be evaluated 
%     % at xtop = 0 since putop_n in this case turns zero. However, the stiff-
%     % ness at the soil surface is zero/low
%     ktmomtop = 0; 
% else
% 	ktmomtop = Mb./0.00001;
% end
% 	ktmomtop=ktmomtop/(element.level(end,1)-element.level(end,2));%SPSO added
% 	ktmombot = Mb./0.00001;
% 	ktmombot=ktmombot/(element.level(end,1)-element.level(end,2));%SPSO added