function [kttop ktbot] = stiffclay_w_o_water_tan(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to stiff clay w/o
% water according to Reese & Van Impe (2001).
%
% INPUT:  springinput   : cf. mainfile [D Su gamma IDcyc/sta eps50 J]
%                         1 = static behaviour (IDcyc/sta)
%                         2 = cyclic behaviour (IDcyc/sta)
%         putot         : capacity [kN/m] in top and bottom of element
%         heqv          : Equivalent length for top and bottom of element
%         zrtot         : Transition depth [m] - moderate to deep
%         u             : Global displacement vector
%         i             : Counter referring to element number
%
% OUTPUT: kttop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : MMOL
% APPROVED        : CTL

% LAST MODIFIED   : MMOL   05.07.2012   Programming
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Initializing pile and soil parameters
% -------------------------------------------------------------------------
D       = pile.diameter;            % Outer diameter [m]
eps50   = element.epsilon50(i);     % strain which occurs at 50% of failure

% -------------------------------------------------------------------------
% Tangent stiffness, kt
% -------------------------------------------------------------------------
% According to the API standard the p-y (or p/pu-y/yc) curves are
% piece-wise linear. The limits are determined according to the API
% standard p.64.

kt = NaN(2,1); %preallocation
for top_bottom = 1:2 %first the top node spring, thereafter the bottom node spring
    %Initializing parameters
    pu = element.pu(i,top_bottom);          % [kN/m] Ultimate resistance, determined in layer.m
    y  = y_topbottom(top_bottom);           % [m] horizontal displacement
    yc = 2.5*eps50*D; % also sometimes referred to as 'y50'
    
    % Initial stiffness
    w_sta = 0.1; % weight to be applied on static initial part to define initial area
    w_cyc = 1; % weight to be applied on cyclic initial part to define initial area
    
    yini_sta = w_sta*yc;
    pini_sta = 0.5*pu*(yini_sta/yc)^0.25;
    kini_sta = pini_sta/yini_sta;
    
    % Determination of tangent stiffness in top of an element according to the API standard p.64.
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
        yini = yini_sta;
        kini = kini_sta;
		loads.n_cycles = 1;
    % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')    
        pini = pini_sta; % MMOL
        yini = w_sta*yc+w_cyc*9.6*(pini/pu)^4*yc*log10(loads.n_cycles); % MMOL
        kini = pini/yini; % MMOL
    end
        
    ylim = 16*yc+9.6*yc*log10(loads.n_cycles); % limit between curved and horisontal part of curve
        
    if y <= yini;
        kt(top_bottom) = kini;          
    elseif y <= ylim
        kt(top_bottom) = 0.25*pu/(yc*(16+9.6*log10(loads.n_cycles))^0.25*(y/yc)^0.75);
    else
        kt(top_bottom) = 0;
    end
end

kttop = kt(1);
ktbot = kt(2);