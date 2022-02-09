function [kstop ksbot] = PISAclaysec(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] and pile resistance p 
% in the top and the bottom of each pile segment by applying p-v curves 
% according to PISA clay model. 
% 
% INPUT:  springinput   : [D G0 su]
%         putop         : Normalized ultimate resistance [-] in top of element, determined in PISAlayer.m
%         pubot         : Normalized ultimate resistance [-] in bot of element, determined in PISAlayer.m
%         heqv          : Equivalent length for top and bottom of element
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% Log:
% EVVA    23.08.2016  Programming
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
su_top  = element.cu(i,1);     		% Undrained shear strength in the top node[kPa]
su_bot  = element.cu(i,2);   		% Undrained shear strength in the bottom node[kPa]
G_top   = element.G0(i,1);          % Shear modulus at the top of the element [kPa]
G_bot   = element.G0(i,2);          % Shear modulus at the bottom of the element [kPa]
putop_n   = element.pu(i,1);        % Normalized ultimate resistance [-] in top
pubot_n   = element.pu(i,2);        % and bottom of pile element, determined in layer.m


% -------- Initial modulus -------------------------------------------------
k_top_n = element.PISA_param(i,4); % Normalized initial stiffness
k_bot_n = element.PISA_param(i,5); % Normalized initial stiffness

k_top = k_top_n.*G_top; % Non-normalized initial stiffness 
k_bot = k_bot_n.*G_bot; % Non-normalized initial stiffness 

% -------- Curvature -------------------------------------------------------
n_top = element.PISA_param(i,6);
n_bot = element.PISA_param(i,7);

% -------- Normalized ultimate lateral displacement -----------------------------------
yu_n  = element.PISA_param(i,1); % Same for Cowden clay and Bothkennar clay

% --------- Normalized horizontal displacement ------------------------------------
ytop_n = ytop/D*(G_top./su_top); %
ybot_n = ybot/D*(G_bot./su_bot); %

% -------- Calculation of parameters a, b and c ----------------------------
a_top_n = 1-2*n_top;
a_bot_n = 1-2*n_bot;
b_top_n = 2*n_top*(ytop_n/yu_n)-(1-n_top)*(1+(ytop_n.*k_top_n./putop_n));%
b_bot_n = 2*n_bot*(ybot_n/yu_n)-(1-n_bot)*(1+(ybot_n.*k_bot_n./pubot_n));%
c_top_n = (ytop_n.*k_top_n./putop_n)*(1-n_top)-n_top*(ytop_n/yu_n)^2;%
c_bot_n = (ybot_n.*k_bot_n./pubot_n)*(1-n_bot)-n_bot*(ybot_n/yu_n)^2;%

% -------- Pile resistance -------------------------------------------------
%top
if ytop_n <= yu_n
ptop_n = putop_n.*((2*c_top_n)./(-b_top_n+sqrt(b_top_n.^2-4*a_top_n*c_top_n))); %Normalized pile resistance %
ptop = ptop_n.*su_top*D; %Non-normalized pile resistance %
else
ptop_n = putop_n; %Normalized ultimate pile resistance 
ptop = ptop_n.*su_top*D; %Non-normalized ultimate pile resistance %
end
%bottom
if ybot_n <= yu_n
pbot_n = pubot_n.*((2*c_bot_n)./(-b_bot_n+sqrt(b_bot_n.^2-4*a_bot_n*c_bot_n))); %Normalized pile resistance %!!
pbot = pbot_n.*su_bot*D; %Non-normalized pile resistance %
else
pbot_n = pubot_n; %Normalized ultimate pile resistance 
pbot = pbot_n.*su_bot*D; %Non-normalized ultimate pile resistance %
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

