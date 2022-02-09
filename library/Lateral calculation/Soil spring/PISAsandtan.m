function [kttop ktbot] = PISAsandtan(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] and pile resistance p 
% in the top and the bottom of each pile segment by applying p-v curves 
% according to PISA sand model. 
% 
% INPUT:  
%         pu (_n)       : Normalized ultimate resistance [-] in top and bottom of element, determined in PISAlayer.m
%         heqv          : Depth from the seabed
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kttop   : Tangent soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Tangent soil stiffness at the bottom of the element [kN/m/m] 
%
% Log: 
% EVVA    09.08.2016  Programming
% FKMV    17.08.2020  Update to use database
%--------------------------------------------------------------------------
% -------- Initializing pile and soil parameters --------------------------

ytop    = abs(y_topbottom(1));      % Horizontal disp. at the top of the
                                    % pile segment [m]

ybot    = abs(y_topbottom(2));      % Horizontal disp. at the bottom of the
                                    % pile segment [m]

xtop    = element.heqv(i,1);        % Distance from the seabed to the top
                                    % of the considered pile segment [m].

xbot    = element.heqv(i,2);        % Distance from the seabed to the bottom
                                    % of the considered pile segment [m].

D       = pile.diameter;            % Outer diameter [m]
G_top   = element.G0(i,1);             % Shear modulus at the top of the element [kPa]
G_bot   = element.G0(i,2);             % Shear modulus at the bottom of the element [kPa]
putop_n   = element.pu(i,1);        % Normalized ultimate resistance [-] in top
pubot_n   = element.pu(i,2);        % and bottom of pile element, determined in PISAlayer.m
sigma_v_top= element.sigma_v_eff(i,1);            % Effective vertical stress at the element top [kPa]
sigma_v_bot= element.sigma_v_eff(i,2);            % Effective vertical stress at the element bottom [kPa]

% % % ---------Calculating G0 ----------------------------------------------------------
% % This is only to replicate PISA results
% B = 875; % Fitted parameter in PISA report [-]
% e = 0.628; % Void ratio  [-]
% K0 = 0.4; % Lateral eath pressure coefficient [-]
% p_ref = 101.325; % Reference pressure in PISA report [kPa]
% p_top = (sigma_v_top+2*K0*sigma_v_top)/3; % Effective mean stress at the element top [kPa]
% p_bot = (sigma_v_bot+2*K0*sigma_v_bot)/3; % Effective mean stress at the element bottom [kPa]

% G_top = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_top/p_ref);
% G_bot = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_bot/p_ref);

% -------- Initial modulus -------------------------------------------------
k_top_n = element.PISA_param(i,4); % Normalized initial stiffness
k_bot_n = element.PISA_param(i,5); % Normalized initial stiffness
% k_top_n = -0.85*(xtop/D)+7.46; % Normalized initial stiffness
% k_bot_n = -0.85*(xbot/D)+7.46; % Normalized initial stiffness
% k_top_n = -0.82*(xtop/D)+7.34; % Normalized initial stiffness - According PISA 2 valid for Dr=45%
% k_bot_n = -0.82*(xbot/D)+7.34; % Normalized initial stiffness - According PISA 2 valid for Dr=45%
% k_top_n = -0.82*(xtop/D)+7.42; % Normalized initial stiffness - According PISA 2 valid for Dr=60%
% k_bot_n = -0.82*(xbot/D)+7.42; % Normalized initial stiffness - According PISA 2 valid for Dr=60%
% k_top_n = -0.83*(xtop/D)+7.31; % Normalized initial stiffness - According PISA 2 valid for Dr=90%
% k_bot_n = -0.83*(xbot/D)+7.31; % Normalized initial stiffness - According PISA 2 valid for Dr=90%
k_top = k_top_n.*G_top; % Non-normalized initial stiffness
k_bot = k_bot_n.*G_bot; % Non-normalized initial stiffness

% -------- Curvature -------------------------------------------------------
n_top = element.PISA_param(i,6);
n_bot = element.PISA_param(i,7);
% n = 0.944;%- According PISA 2 valid for Dr=75%
% n = 0.940;%- According PISA 2 valid for Dr=45%
% n = 0.950;%- According PISA 2 valid for Dr=60%
% n = 0.962;%- According PISA 2 valid for Dr=90%

% -------- Normalized ultimate lateral displacement -----------------------------------
yu_n  = element.PISA_param(i,1); % Same for Cowden clay and Bothkennar clay
% yu_n  = 53.1;%- According PISA 2 valid for Dr=75%
% yu_n  = 102.4;%- According PISA 2 valid for Dr=45%
% yu_n  = 75.8;%- According PISA 2 valid for Dr=60%
% yu_n  = 58.93;%- According PISA 2 valid for Dr=90%

% --------- Normalized horizontal displacement ------------------------------------
% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dp/dy [kN/m^2] for y >0. Real expression 
% for dp/dy is very complex; thus, dy=0.00001m, p2=p(y+0.00001) p1=p(y) and kt=(p(y+0.00001)-p(y))/0.00001.

% Parameter y1
ytop_n = ytop/D*(G_top./sigma_v_top);
ybot_n = ybot/D*(G_bot./sigma_v_bot);
% Parameter y2
ytop_n2 = (ytop+0.00001)/D*(G_top./sigma_v_top);
ybot_n2 = (ybot+0.00001)/D*(G_bot./sigma_v_bot);


% -------- Calculation of parameters a, b and c ----------------------------
% a = 1-2*n;
a_top = 1-2*n_top;
a_bot = 1-2*n_bot;

% Parameter b1
b_top_n = 2*n_top*(ytop_n/yu_n)-(1-n_top)*(1+(ytop_n.*k_top_n./putop_n));
b_bot_n = 2*n_bot*(ybot_n/yu_n)-(1-n_bot)*(1+(ybot_n.*k_bot_n./pubot_n));
% Parameter b2
b_top_n2 = 2*n_top*(ytop_n2/yu_n)-(1-n_top)*(1+(ytop_n2.*k_top_n./putop_n));
b_bot_n2 = 2*n_bot*(ybot_n2/yu_n)-(1-n_bot)*(1+(ybot_n2.*k_bot_n./pubot_n));

% Parameter c1
c_top_n = (ytop_n.*k_top_n./putop_n)*(1-n_top)-n_top*(ytop_n/yu_n)^2;
c_bot_n = (ybot_n.*k_bot_n./pubot_n)*(1-n_bot)-n_bot*(ybot_n/yu_n)^2;
% Parameter c2
c_top_n2 = (ytop_n2.*k_top_n./putop_n)*(1-n_top)-n_top*(ytop_n2/yu_n)^2;
c_bot_n2 = (ybot_n2.*k_bot_n./pubot_n)*(1-n_bot)-n_bot*(ybot_n2/yu_n)^2;


% -------- Pile resistance -------------------------------------------------
% top
% calculating p1, normalized and non-normalized
if ytop_n <= yu_n
ptop_n1 = putop_n.*((2*c_top_n)./(-b_top_n+sqrt(b_top_n.^2-4*a_top*c_top_n))); %Normalized pile resistance 
ptop1 = ptop_n1.*sigma_v_top*D; %Non-normalized pile resistance 
else
ptop_n1 = putop_n; %Normalized ultimate pile resistance 
ptop1 = ptop_n1.*sigma_v_top*D; %Non-normalized ultimate pile resistance 
end
% calculating p2, normalized and non-normalized
if ytop_n2 <= yu_n
ptop_n2 = putop_n.*((2*c_top_n2)./(-b_top_n2+sqrt(b_top_n2.^2-4*a_top*c_top_n2))); %Normalized pile resistance 
ptop2 = ptop_n2.*sigma_v_top*D; %Non-normalized pile resistance 
else
ptop_n2 = putop_n; %Normalized ultimate pile resistance 
ptop2 = ptop_n2.*sigma_v_top*D; %Non-normalized ultimate pile resistance 
end

ptop=ptop2-ptop1;


% bottom
% calculating p1, normalized and non-normalized
if ybot_n <= yu_n
pbot_n1 = pubot_n.*((2*c_bot_n)./(-b_bot_n+sqrt(b_bot_n.^2-4*a_bot*c_bot_n))); %Normalized pile resistance 
pbot1 = pbot_n1.*sigma_v_bot*D; %Non-normalized pile resistance 
else
pbot_n1 = pubot_n; %Normalized ultimate pile resistance 
pbot1 = pbot_n1.*sigma_v_bot*D; %Non-normalized ultimate pile resistance 
end
% calculating p2, normalized and non-normalized
if ybot_n2 <= yu_n
pbot_n2 = pubot_n.*((2*c_bot_n2)./(-b_bot_n2+sqrt(b_bot_n2.^2-4*a_bot*c_bot_n2))); %Normalized pile resistance 
pbot2 = pbot_n2.*sigma_v_bot*D; %Non-normalized pile resistance 
else
pbot_n2 = pubot_n; %Normalized ultimate pile resistance 
pbot2 = pbot_n2.*sigma_v_bot*D; %Non-normalized ultimate pile resistance 
end

pbot=pbot2-pbot1;

% -------- Tangent stiffness, kttop and ktbot -------------------------------

if xtop == 0;
    % This statement is included because Kttop = dp/dy cannot be evaluated 
    % at xtop = 0 since putop_n in this case turns zero. However, the stiff-
    % ness at the soil surface is zero/low
    kttop = 0; 
else
	kttop = ptop/0.00001;
end

	ktbot = pbot/0.00001;

