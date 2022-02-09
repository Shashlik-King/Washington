function [ktmomtop ktmombot]  = PISAclaytanm(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i)
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
L   = abs(element.level(end,1));	   % pile.length_end;
G_top   = element.G0(i,1);            % Shear modulus at the top of the element [kPa]
G_bot   = element.G0(i,2);            % Shear modulus at the bottom of the element [kPa]
mutop_n   = element.mu(i,1);             % Normalized ultimate resistance [-] in top
mubot_n   = element.mu(i,2);             % and bottom of pile element, determined in momlayer.m
sigma_v_top= element.sigma_v_eff(i,1);   % Effective vertical stress at the element top [kPa]
sigma_v_bot= element.sigma_v_eff(i,2);   % Effective vertical stress at the element bottom [kPa]
su_top   = element.cu(i,1);              % Undrained shear strength in the top node[kPa]
su_bot   = element.cu(i,2);              % Undrained shear strength in the bottom node[kPa]

% -------- Initial modulus -------------------------------------------------
k_top_n = element.PISA_param(i,11); % Normalized initial stiffness
k_bot_n = element.PISA_param(i,12); % Normalized initial stiffness

k_top = k_top_n*D^2*G_top; % Non-normalized initial stiffness
k_bot = k_bot_n*D^2*G_bot;  % Non-normalized initial stiffness

% -------- Curvature ------------------------------------------------------
ntop = element.PISA_param(i,13);
nbot = element.PISA_param(i,14);

% -------- Normalized ultimate nodal rotation -----------------------------
teta_u_n  = element.PISA_param(i,8); % Sane for Cowden and Bothkennar clay

% --------- Normalized nodal rotation -------------------------------------
% Parameter teta1
tetatop_n1 = tetatop.*(G_top./su_top);
tetabot_n1 = tetabot.*(G_bot./su_bot);

% Parameter teta2
tetatop_n2 = (tetatop+0.00001).*(G_top./su_top);
tetabot_n2 = (tetabot+0.00001).*(G_bot./su_bot);

% -------- Calculation of parameters a, b and c ---------------------------
a = 1-2*ntop;
atop = 1-2*ntop;
abot = 1-2*nbot;
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
mtop_n1 = mutop_n.*((2*c_top_n1)./(-b_top_n1+sqrt(b_top_n1.^2-4*atop*c_top_n1))); %Normalized distr. moment  
mtop1 = mtop_n1.*su_top*D^2; %Non-normalized distr. moment  
else
mtop_n1 = mutop_n; %Normalized ultimate distr. moment 
mtop1 = mtop_n1.*su_top*D^2; %Non-normalized ultimate distr. moment  
end
% top, calculating m2 at the top of the element, normalized and non-normalized
if tetatop_n2 <= teta_u_n
mtop_n2 = mutop_n.*((2*c_top_n2)./(-b_top_n2+sqrt(b_top_n2.^2-4*atop*c_top_n2))); %Normalized distr. moment 
mtop2 = mtop_n2.*su_top*D^2; %Non-normalized distr. moment 
else
mtop_n2 = mutop_n; %Normalized ultimate distr. moment
mtop2 = mtop_n2.*su_top*D^2; %Non-normalized ultimate distr. moment 
end

mtop=mtop2-mtop1;

% bottom, calculating m1 at the bottom of the element, normalized and non-normalized
if tetabot_n1 <= teta_u_n
mbot_n1 = mubot_n.*((2*c_bot_n1)./(-b_bot_n1+sqrt(b_bot_n1.^2-4*abot*c_bot_n1))); %Normalized distr. moment 
mbot1 = mbot_n1.*su_bot*D^2; %Non-normalized distr. moment 
else
mbot_n1 = mubot_n; %Normalized ultimate distr. moment 
mbot1 = mbot_n1.*su_bot*D^2; %Non-normalized ultimate distr. moment 
end
% bottom, calculating m2 at the bottom of the element, normalized and non-normalized
if tetabot_n2 <= teta_u_n
mbot_n2 = mubot_n.*((2*c_bot_n2)./(-b_bot_n2+sqrt(b_bot_n2.^2-4*abot*c_bot_n2))); %Normalized distr. moment 
mbot2 = mbot_n2.*su_bot*D^2; %Non-normalized distr. moment 
else
mbot_n2 = mubot_n; %Normalized ultimate distr. moment 
mbot2 = mbot_n2.*su_bot*D^2; %Non-normalized ultimate distr. moment 
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



 