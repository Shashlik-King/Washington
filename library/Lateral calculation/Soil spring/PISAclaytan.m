function [kttop ktbot] = PISAclaytan(element,pile,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] and pile resistance p 
% in the top and the bottom of each pile segment by applying p-v curves 
% according to PISA sand model. 
% 
% INPUT:  springinput   : [D G0 su]
%         pu (_n)       : Normalized ultimate resistance [-] in top bottom of element, determined in PISAlayer.m
%         heqv          : Embedment
%         u             : Global displacement vector
%         i             : Counter referring to element number  
%
% OUTPUT: kttop   : Tangent soil stiffness at the top of the element [kN/m/m]
%         ktbot   : Tangent soil stiffness at the bottom of the element [kN/m/m]
%
% Log:
%  EVVA    10.08.2016  Programming
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
pubot_n   = element.pu(i,2);        % and bottom of pile element, determined in PISAlayer.m

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
% Determination of spring stiffness at the top and bottom of each pile
% segment, kttop is the tangent stiffness dp/dy [kN/m^2] for y >0. Real expression 
% for dp/dy is very complex; thus, dy=0.00001m, p2=p(y+0.00001) p1=p(y) and kt=(p(y+0.00001)-p(y))/0.00001.
% Parameter y1
ytop_n = ytop/D*(G_top./su_top);
ybot_n = ybot/D*(G_bot./su_bot);

% Parameter y2
ytop_n2 = (ytop+0.00001)/D*(G_top./su_top);
ybot_n2 = (ybot+0.00001)/D*(G_bot./su_bot);

% -------- Calculation of parameters a, b and c ----------------------------
a_top_n = 1-2*n_top;
a_bot_n = 1-2*n_bot;

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
ptop_n1 = putop_n.*((2*c_top_n)./(-b_top_n+sqrt(b_top_n.^2-4*a_top_n*c_top_n))); %Normalized pile resistance 
ptop1 = ptop_n1.*su_top*D; %Non-normalized pile resistance 
else
ptop_n1 = putop_n; %Normalized ultimate pile resistance 
ptop1 = ptop_n1.*su_top*D; %Non-normalized ultimate pile resistance 
end
% calculating p2, normalized and non-normalized
if ytop_n2 <= yu_n
ptop_n2 = putop_n.*((2*c_top_n2)./(-b_top_n2+sqrt(b_top_n2.^2-4*a_top_n*c_top_n2))); %Normalized pile resistance 
ptop2 = ptop_n2.*su_top*D; %Non-normalized pile resistance 
else
ptop_n2 = putop_n; %Normalized ultimate pile resistance 
ptop2 = ptop_n2.*su_top*D; %Non-normalized ultimate pile resistance 
end

ptop=ptop2-ptop1;


% bottom
% calculating p1, normalized and non-normalized
if ybot_n <= yu_n
pbot_n1 = pubot_n.*((2*c_bot_n)./(-b_bot_n+sqrt(b_bot_n.^2-4*a_bot_n*c_bot_n))); %Normalized pile resistance 
pbot1 = pbot_n1.*su_bot*D; %Non-normalized pile resistance 
else
pbot_n1 = pubot_n; %Normalized ultimate pile resistance 
pbot1 = pbot_n1.*su_bot*D; %Non-normalized ultimate pile resistance 
end
% calculating p2, normalized and non-normalized
if ybot_n2 <= yu_n
pbot_n2 = pubot_n.*((2*c_bot_n2)./(-b_bot_n2+sqrt(b_bot_n2.^2-4*a_bot_n*c_bot_n2))); %Normalized pile resistance 
pbot2 = pbot_n2.*su_bot*D; %Non-normalized pile resistance 
else
pbot_n2 = pubot_n; %Normalized ultimate pile resistance 
pbot2 = pbot_n2.*su_bot*D; %Non-normalized ultimate pile resistance 
end

pbot=pbot2-pbot1;

% -------- Tangent stiffness, kttop and ktbot -------------------------------

if xtop == 0;
    % This statement is included because Kttop = dp/dy cannot be evaluated 
    % at xtop = 0 since putop in this case turns zero. However, the stiff-
    % ness at the soil surface is zero
    kttop = 0; 
else
	kttop = ptop/0.00001;

end

	ktbot = pbot/0.00001;
