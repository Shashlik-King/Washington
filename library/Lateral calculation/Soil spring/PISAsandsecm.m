function [ksmomtop ksmombot]  = PISAsandsecm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i)
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant uniformly distributed moment stiffness [kNm/m] in the top and the bottom of
% each pile segment by applying m-teta curves.
% 
% INPUT:  mu            : Normalized ultimate distr. moment capacity [-] in top and bottom of element
%         kstop ksbot   : Secant stiffness of the uniformly distributed
%                         load p [kN/m/m]
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: ksmomtop   : Secant soil moment stiffness at the top of the element [kNm/m/m]
%         ksmombot   : Secant soil moment stiffness at the bottom of the element [kNm/m/m]
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
G_top   = element.G0(i,1);               % Shear modulus at the top of the element [kPa]
G_bot   = element.G0(i,2);               % Shear modulus at the bottom of the element [kPa]
mutop_n   = element.mu(i,1);             % Normalized ultimate resistance [-] in top
mubot_n   = element.mu(i,2);             % and bottom of pile element, determined in PISAlayer.m
sigma_v_top= element.sigma_v_eff(i,1);   % Effective vertical stress at the element top [kPa]
sigma_v_bot= element.sigma_v_eff(i,2);   % Effective vertical stress at the element bottom [kPa]
su_top   = element.cu(i,1);     	     % Undrained shear strength in the top node[kPa]
su_bot   = element.cu(i,2);   		     % Undrained shear strength in the bottom node[kPa]

% % ---------Calculating G0 -------------------------------------------------
% % This is only to replicate PISA results
% B = 875; % Fitted parameter in PISA report [-]
% e = 0.628; % Void ratio  [-]
% K0 = 0.4; % Lateral eath pressure coefficient [-]
% p_ref = 101.3; % Reference pressure in PISA report [kPa]
% p_top = (sigma_v_top+2*K0*sigma_v_top)/3; % Effective mean stress at the element top [kPa]
% p_bot = (sigma_v_bot+2*K0*sigma_v_bot)/3; % Effective mean stress at the element bottom [kPa]

% G_top = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_top/p_ref); % PISA report pg. 83/269
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
% teta_u_n  = 20; % - According to PISA 2 valid for Dr=45, 60, 75 and 90%
teta_u_n  = element.PISA_param(i,8); % Sane for Cowden and Bothkennar clay

% --------- Normalized nodal rotation -------------------------------------
tetatop_n = tetatop.*(G_top./sigma_v_top);
tetabot_n = tetabot.*(G_bot./sigma_v_bot);


% -------- Calculation of parameters a, b and c ---------------------------
% a = 1-2*n;
a_top = 1-2*ntop;
a_bot = 1-2*nbot;
b_top_n = 2*ntop*(tetatop_n/teta_u_n)-(1-ntop)*(1+(tetatop_n.*k_top_n./mutop_n));
b_bot_n = 2*nbot*(tetabot_n/teta_u_n)-(1-nbot)*(1+(tetabot_n.*k_bot_n./mubot_n));
c_top_n = (tetatop_n.*k_top_n./mutop_n)*(1-ntop)-ntop*(tetatop_n./teta_u_n).^2;
c_bot_n = (tetabot_n.*k_bot_n./mubot_n)*(1-nbot)-nbot*(tetabot_n./teta_u_n).^2;

% -------- distr. moment -------------------------------------------------
%top
if tetatop_n <= teta_u_n
mtop_n = mutop_n.*((2*c_top_n)./(-b_top_n+sqrt(b_top_n.^2-4*a_top*c_top_n))); %Normalized distr. moment 
mtop = mtop_n.*ptop*D; %Non-normalized distr. moment 
else
mtop_n = mutop_n; %Normalized ultimate distr. moment
mtop = mtop_n.*ptop*D; %Non-normalized ultimate distr. moment 
end
%bottom
if tetabot_n <= teta_u_n
mbot_n = mubot_n.*((2*c_bot_n)./(-b_bot_n+sqrt(b_bot_n.^2-4*a_bot*c_bot_n))); %Normalized distr. moment 
mbot = mbot_n.*pbot*D; %Non-normalized distr. moment 
else
mbot_n = mubot_n; %Normalized ultimate distr. moment 
mbot = mbot_n.*pbot*D; %Non-normalized ultimate distr. moment 
end


% -------- Secant stiffness, ksmomtop and ksmombot -------------------------------
% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness m/teta [kNm/m/rad] for teta >0. For teta = 0
% ksmomtop is the tangent stiffness, i.e. ksmomtop = dm/dteta evaluated at teta = 0. 
% Similar for ksmombot

if tetatop == 0;
    %ksmomtop = k_bot*xtop;%gives infinity for 0 displacement and 0 depth
    ksmomtop = 0;
elseif xtop == 0;
    ksmomtop = 0;
else
	ksmomtop = mtop/tetatop;
end

if tetabot == 0;
    ksmombot = k_bot*xbot;
else
	ksmombot = mbot/tetabot;
end
 
ksmomtop;
ksmombot;


 