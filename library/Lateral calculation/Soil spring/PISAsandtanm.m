function [ktmomtop ktmombot]  = PISAsandtanm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent uniformly distributed moment stiffness [kNm/m] in the top and the bottom of
% each pile segment by applying m-teta curves.
% 
% INPUT:  mu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%						  load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ktmomtop   : Tangent soil moment stiffness at the top of the element [kNm/m/m]
%         ktmombot   : Tanmgent soil moment stiffness at the bottom of the element [kNm/m/m]
%
% Log: 
% EVVA    17.08.2016   Setting up the function
% FKMV    17.08.2020  Update to use database
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------

tetatop    = abs(teta_topbottom(1));     % Nodal rotation at the top of the
                                         % pile segment [rad]

tetabot    = abs(teta_topbottom(2));     % Nodal rotation at the bottom of the
                                         % pile segment [rad]
                                    
ytop    = abs(y_topbottom(1));           % Horizontal disp. at the top of the
                                         % pile segment [m]

ybot    = abs(y_topbottom(2));           % Horizontal disp. at the bottom of the
                                         % pile segment [m]

xtop    = element.heqv(i,1);             % Distance from the seabed to the top
                                         % of the considered pile segment [m].

xbot    = element.heqv(i,2);             % Distance from the seabed to the bottom
                                         % of the considered pile segment [m].

D       = pile.diameter;                 % Outer diameter [m]
G_top   = element.G0(i,1);                 % Shear modulus at the top of the element [kPa]
G_bot   = element.G0(i,2);                 % Shear modulus at the bottom of the element [kPa]
mutop_n   = element.mu(i,1);             % Normalized ultimate resistance [-] in top
mubot_n   = element.mu(i,2);             % and bottom of pile element, determined in momlayer.m
sigma_v_top= element.sigma_v_eff(i,1);   % Effective vertical stress at the element top [kPa]
sigma_v_bot= element.sigma_v_eff(i,2);   % Effective vertical stress at the element bottom [kPa]
su_top   = element.cu(i,1);              % Undrained shear strength in the top node[kPa]
su_bot   = element.cu(i,2);              % Undrained shear strength in the bottom node[kPa]

% % ---------Calculating G0 ----------------------------------------------------------
% % This is only to replicate PISA results
% B = 875; % Fitted parameter in PISA report [-]
% e = 0.628; % Void ratio  [-]
% K0 = 0.4; % Lateral eath pressure coefficient [-]
% p_ref = 101.3; % Reference pressure in PISA report [kPa]
% p_top = (sigma_v_top+2*K0*sigma_v_top)/3; % Effective mean stress at the element top [kPa]
% p_bot = (sigma_v_bot+2*K0*sigma_v_bot)/3; % Effective mean stress at the element bottom [kPa]

% G_top = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_top/p_ref);% PISA report pg. 83/269
% G_bot = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_bot/p_ref);

% -------- Non-normalized distributed load---------------------------------
ptop = kstop.*ytop;
pbot = ksbot.*ybot;

% -------- Initial modulus -------------------------------------------------
% k_n = 20; % Normalized initial stiffness - According to PISA 2 valid for Dr=45, 60, 75 and 90%
k_top_n = element.PISA_param(i,11); % Normalized initial stiffness
k_bot_n = element.PISA_param(i,12); % Normalized initial stiffness
k_top = k_top_n*D*(ptop.*G_top)./sigma_v_top; % Non-normalized initial stiffness
k_bot = k_bot_n*D*(pbot.*G_bot)./sigma_v_bot; % Non-normalized initial stiffness

% -------- Curvature ------------------------------------------------------
% n = 0;% - According to PISA 2 valid for Dr=45, 60, 75 and 90%
ntop = element.PISA_param(i,13);
nbot = element.PISA_param(i,14);

% -------- Normalized ultimate nodal rotation -----------------------------
% teta_u_n  = 20;% - According to PISA 2 valid for Dr=45, 60, 75 and 90%
teta_u_n  = element.PISA_param(i,8); % Sane for Cowden and Bothkennar clay

% --------- Normalized nodal rotation  ------------------------------------
% Determination of moment stiffness at the top and bottom of each pile
% segment, ktmomtop is the tangent stiffness dm/dteta [kNm/m/rad] for teta >0. Real expression 
% for dm/dteta is very complex; thus, dteta=0.0001rad, m2=m(teta+0.0001) m1=m(teta) and kt=(m(teta+0.0001)-m(teta))/0.0001.

% Parameter teta1
tetatop_n1 = tetatop.*(G_top./sigma_v_top);
tetabot_n1 = tetabot.*(G_bot./sigma_v_bot);

% Parameter teta2
tetatop_n2 = (tetatop+0.0001).*(G_top./sigma_v_top);
tetabot_n2 = (tetabot+0.0001).*(G_bot./sigma_v_bot);


% -------- Calculation of parameters a, b and c ---------------------------
% a = 1-2*n;
a_top = 1-2*ntop;
a_bot = 1-2*nbot;

% Parameter b1
b_top_n1 = 2*ntop*(tetatop_n1/teta_u_n)-(1-ntop)*(1+(tetatop_n1.*k_top_n./mutop_n));
b_bot_n1 = 2*nbot*(tetabot_n1/teta_u_n)-(1-nbot)*(1+(tetabot_n1.*k_bot_n./mubot_n));
% Parameter b2
b_top_n2 = 2*ntop*(tetatop_n2/teta_u_n)-(1-ntop)*(1+(tetatop_n2.*k_top_n./mutop_n));
b_bot_n2 = 2*nbot*(tetabot_n2/teta_u_n)-(1-nbot)*(1+(tetabot_n2.*k_bot_n./mubot_n));

% Parameter c1
c_top_n1 = (tetatop_n1.*k_top_n./mutop_n)*(1-ntop)-ntop*(tetatop_n1./teta_u_n).^2;
c_bot_n1 = (tetabot_n1.*k_bot_n./mubot_n)*(1-nbot)-nbot*(tetabot_n1./teta_u_n).^2;
% Parameter c2
c_top_n2 = (tetatop_n2.*k_top_n./mutop_n)*(1-ntop)-ntop*(tetatop_n2./teta_u_n).^2;
c_bot_n2 = (tetabot_n2.*k_bot_n./mubot_n)*(1-nbot)-nbot*(tetabot_n2./teta_u_n).^2;

% -------- distr. moment -------------------------------------------------
% top, calculating m1 at the top of the element, normalized and non-normalized
if tetatop_n1 <= teta_u_n
mtop_n1 = mutop_n.*((2*c_top_n1)./(-b_top_n1+sqrt(b_top_n1.^2-4*a_top*c_top_n1))); %Normalized distr. moment resistance 
mtop1 = mtop_n1.*ptop*D; %Non-normalized distr. moment resistance 
else
mtop_n1 = mutop_n; %Normalized ultimate distr. moment resistance
mtop1 = mtop_n1.*ptop*D; %Non-normalized ultimate distr. moment resistance 
end
% top, calculating m2 at the top of the element, normalized and non-normalized
if tetatop_n2 <= teta_u_n
mtop_n2 = mutop_n.*((2*c_top_n2)./(-b_top_n2+sqrt(b_top_n2.^2-4*a_top*c_top_n2))); %Normalized distr. moment resistance 
mtop2 = mtop_n2.*ptop*D; %Non-normalized distr. moment resistance 
else
mtop_n2 = mutop_n; %Normalized ultimate distr. moment resistance
mtop2 = mtop_n2.*ptop*D; %Non-normalized ultimate distr. moment resistance 
end

mtop=mtop2-mtop1;

% bottom, calculating m1 at the bottom of the element, normalized and non-normalized
if tetabot_n1 <= teta_u_n
mbot_n1 = mubot_n.*((2*c_bot_n1)./(-b_bot_n1+sqrt(b_bot_n1.^2-4*a_bot*c_bot_n1))); %Normalized distr. moment resistance 
mbot1 = mbot_n1.*pbot*D; %Non-normalized distr. moment resistance 
else
mbot_n1 = mubot_n; %Normalized ultimate distr. moment resistance 
mbot1 = mbot_n1.*pbot*D; %Non-normalized ultimate distr. moment resistance 
end
% bottom, calculating m2 at the bottom of the element, normalized and non-normalized
if tetabot_n2 <= teta_u_n
mbot_n2 = mubot_n.*((2*c_bot_n2)./(-b_bot_n2+sqrt(b_bot_n2.^2-4*a_bot*c_bot_n2))); %Normalized distr. moment resistance 
mbot2 = mbot_n2.*pbot*D; %Non-normalized distr. moment resistance 
else
mbot_n2 = mubot_n; %Normalized ultimate distr. moment resistance 
mbot2 = mbot_n2.*pbot*D; %Non-normalized ultimate distr. moment resistance 
end

mbot=mbot2-mbot1;

% -------- Tangent stiffness, kttop and ktbot -------------------------------

if xtop == 0;
    % This statement is included because Kttop = dm/dteta cannot be evaluated 
    % at xtop = 0 since mutop_n in this case turns zero. However, the stiff-
    % ness at the soil surface is zero
    ktmomtop = 0; 
else
	ktmomtop = mtop/0.00001;
end

	ktmombot = mbot/0.00001;

ktmomtop;



 