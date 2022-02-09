function [kstop ksbot] = kirschstiffclaysec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves for stiff clay without free water 
% according to Kirsch (2014).
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
% OUTPUT: kstop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : MUOE
% APPROVED        : MMOL
%
% LAST MODIFIED   : MUOE   30.09.2015   Programming
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Initializing pile and soil parameters
% -------------------------------------------------------------------------
D       = pile.diameter;            % Outer diameter [m]
eps50   = element.epsilon50(i);     % strain which occurs at 50% of failure

% -------------------------------------------------------------------------
% Secant stiffness, ks
% -------------------------------------------------------------------------

ks = NaN(2,1); %preallocation
for top_bottom = 1:2 %first the top node spring, thereafter the bottom node spring
    %Initializing parameters
    pu = element.pu(i,top_bottom);  % [kN/m] Ultimate resistance, determined in layer.m
    y  = y_topbottom(top_bottom);   % [m] horizontal displacement
    Es = element.Es(i,1)/1E3; % [MPa] Large-strain constrained modulus 
    Esd =10^(-0.42*log10(0.0006*Es))*Es; % [MPa] Small-strain constrained modulus
    yc = 2.5*eps50*D;
    
    % Initial stiffness
    w_sta = 0.1; % weight to be applied on static initial part to define initial area
    w_cyc = 1; % weight to be applied on cyclic initial part to define initial area
    
    tol = 0.01*eps50; %SPSO modified
    
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
        
        loads.n_cycles = 1;
        
    % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')
        
    end
        
    % iteration to find eps50mod
    e50mod = eps50;
    p = pu*(y/yc)^0.25/(16+9.6*log10(loads.n_cycles))^0.25; %SPSO new
    eps50mod = eps50*(1+(1-p/pu)*(Es/Esd-1)); %SPSO new
    while abs(e50mod-eps50mod) > tol
        e50mod = (e50mod+eps50mod)/2;
        ycmod = 2.5*e50mod*D;
        p = pu*(y/ycmod)^0.25/(16+9.6*log10(loads.n_cycles))^0.25; %SPSO new
        eps50mod = eps50*(1+(1-p/pu)*(Es/Esd-1)); %SPSO new
    end  
    ycmod = 2.5*eps50mod*D;
    yini_sta = w_sta*ycmod;
    pini = 0.5*pu*(yini_sta/ycmod)^0.25; % MMOL
    eps50_ini = eps50*(1+(1-pini/pu)*(Es/Esd-1)); %SPSO
    yc50_ini = eps50_ini*2.5*D; %SPSO
    yini = w_sta*yc50_ini+w_cyc*9.6*(pini/pu)^4*yc50_ini*log10(loads.n_cycles); % MMOL + SPSO modified
    kini = pini/yini; % MMOL

    ylim = 16*ycmod+9.6*ycmod*log10(loads.n_cycles); % limit between curved and horisontal part of curve

    if y <= yini;
        ks(top_bottom) = kini; 
    elseif y <= ylim
        p = pu*(y/ycmod)^0.25/(16+9.6*log10(loads.n_cycles))^0.25;
        ks(top_bottom) = p/y;
    else
        p = pu;
        ks(top_bottom) = p/y;
    end
end

kstop = ks(1);
ksbot = ks(2);