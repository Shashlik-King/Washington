function ts = toe_shear(element,node,pile,ybot,i)
%--------------------------------------------------------------------------
% PURPOSE
% Calculate toe shear
%
% LAST MODIFIED   : MMOL   ??           Programming
%                   JALY   03.05.2011   Streamlining code
%--------------------------------------------------------------------------
D   = pile.diameter;
A   = pi*(D/2)^2; % Cross sectional area. Used for calculating the toe shear
R   = D/2; % Radius of foundation [m]

msgbox('You have activated the function for toe shear - please check that the toe shear model is applicable for your specific project conditions, see toe_shear.m')

if strcmp(element.model_py(i),'API sand') || strcmp(element.model_py(i),'Kirsch sand') || strcmp(element.model_py(i),'Kallehave sand')
    phi         = element.phi(i);
    sigma_toe   = node.sigma_v_eff(i+1); % Effective vertical stress at pile bottom
    m           = 1000*tand(phi); % acc. to DNV-OS-J101
    G0          = m*sqrt(100*sigma_toe)/(2*(1+0.3)); % acc. to DNV-OS-J101
    pu_cyc      = 0.8*A*sigma_toe*tand(phi);
    p_ts        = 1.1*G0*R*ybot;
elseif strcmp(element.model_py(i),'API clay') || strcmp(element.model_py(i),'Stiff clay w/o free water') || strcmp(element.model_py(i),'Kirsch soft clay') || strcmp(element.model_py(i),'Kirsch stiff clay') || strcmp(element.model_py(i),'Reese stiff clay')
    su          = element.cu(i,2);
    OCR         = 2;
    G0          = 600*su-170*su*sqrt(OCR-1);
    pu_cyc      = 0.5*A*su;
    p_ts        = 2.5*G0*R*ybot;
else
	disp(['Toe shear is not supported for ',element.model_py(i)]);
	p_ts 		= 0;
end

if p_ts<0
    ts     = max(-pu_cyc,p_ts);
else
    ts     = min(pu_cyc,p_ts);
end