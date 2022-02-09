function [kstop ksbot] = PISAsandsec(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] and pile resistance p 
% in the top and the bottom of each pile segment by applying p-v curves 
% according to PISA sand model. 
% 
% INPUT:  
%         putop_n       : Normalized ultimate resistance [-] in top of element, determined in PISAlayer.m
%         pubot_n       : Normalized ultimate resistance [-] in bot of element, determined in PISAlayer.m
%         heqv          : Depth from the seabed [m]
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
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
sigma_v_top = element.sigma_v_eff(i,1);           % Effective vertical stress at the element top [kPa]
sigma_v_bot = element.sigma_v_eff(i,2);           % Effective vertical stress at the element bottom [kPa]

% % % ---------Calculating G0 ----------------------------------------------------------
% % % This is only to replicate PISA results
% B = 875;        % Fitted parameter in PISA report [-]
% e = 0.628;      % Void ratio  [-]
% K0 = 0.4;       % Lateral eath pressure coefficient [-]
% p_ref = 101.325;  % Reference pressure in PISA report [kPa]
% p_top = (sigma_v_top+2*K0*sigma_v_top)/3;   % Effective mean stress at the element top [kPa]
% p_bot = (sigma_v_bot+2*K0*sigma_v_bot)/3;   % Effective mean stress at the element bottom [kPa]

% G_top = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_top/p_ref); % PISA report pg. 83/269
% G_bot = (B*p_ref)/(0.3+0.7*e^2)*sqrt(p_bot/p_ref);


% -------- Initial modulus -------------------------------------------------
k_top_n = element.PISA_param(i,4); % Normalized initial stiffness
k_bot_n = element.PISA_param(i,5); % Normalized initial stiffness
% k_top_n = -0.85*(xtop/D)+7.46; % Normalized initial stiffness - According PISA 2 valid for Dr=75%
% k_bot_n = -0.85*(xbot/D)+7.46; % Normalized initial stiffness - According PISA 2 valid for Dr=75%
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
ytop_n = ytop/D*(G_top./sigma_v_top); 
ybot_n = ybot/D*(G_bot./sigma_v_bot); 


% -------- Calculation of parameters a, b and c ----------------------------
% a = 1-2*n;
a_top = 1-2*n_top;
a_bot = 1-2*n_bot;
b_top_n = 2*n_top*(ytop_n/yu_n)-(1-n_top)*(1+(ytop_n.*k_top_n./putop_n));
b_bot_n = 2*n_bot*(ybot_n/yu_n)-(1-n_bot)*(1+(ybot_n.*k_bot_n./pubot_n));
c_top_n = (ytop_n.*k_top_n./putop_n)*(1-n_top)-n_top*(ytop_n/yu_n)^2;
c_bot_n = (ybot_n.*k_bot_n./pubot_n)*(1-n_bot)-n_bot*(ybot_n/yu_n)^2;

% -------- Pile resistance -------------------------------------------------
%top
if ytop_n <= yu_n
ptop_n = putop_n.*((2*c_top_n)./(-b_top_n+sqrt(b_top_n.^2-4*a_top*c_top_n))); %Normalized pile resistance 
ptop = ptop_n.*sigma_v_top*D; %Non-normalized pile resistance 
else
ptop_n = putop_n; %Normalized ultimate pile resistance 
ptop = ptop_n.*sigma_v_top*D; %Non-normalized ultimate pile resistance 
end
%bottom
if ybot_n <= yu_n
pbot_n = pubot_n.*((2*c_bot_n)./(-b_bot_n+sqrt(b_bot_n.^2-4*a_bot*c_bot_n))); %Normalized pile resistance 
pbot = pbot_n.*sigma_v_bot*D; %Non-normalized pile resistance 
else
pbot_n = pubot_n; %Normalized ultimate pile resistance 
pbot = pbot_n.*sigma_v_bot*D; %Non-normalized ultimate pile resistance 
end

% -------- Secant stiffness, kstop and ksbot -------------------------------
% Determination of spring stiffness at the top and bottom of each pile
% segment, kstop is the secant stiffness p/y [kN/m^2] for y >0. For y = 0
% kstop is the tangent stiffness, i.e. kstop = dp/dy evaluated at y = 0. 
% Similar for kbot

if ytop == 0;
    kstop = k_top.*xtop; 
elseif xtop == 0;
    kstop = 0;
else
	kstop = ptop./ytop;
end

if ybot == 0;
    ksbot = k_bot.*xbot;
else
	ksbot = pbot./ybot;
end
