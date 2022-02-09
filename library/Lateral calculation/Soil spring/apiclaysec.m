function [kstop ksbot]   = apiclaysec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to API(1993)- soft
% clay.
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
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    19.11.2007   Programming
%                   AHA    22.05.2008   Introduction of layered soil.
%                                       Correction line 149 putop=>pubot
%                   AHA    26.05.2008   Cyclic behaviour associated with
%                                       soft clay criterion
%                   CTL    05.11.2009   Change to curve description of p-y formulation
%                   MMOL   04.11.2010   Integrating into calculation routine
%                   JALY   29.04.2011   Streamlining code
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Initializing pile and soil parameters
% -------------------------------------------------------------------------
D       = pile.diameter;            % Outer diameter [m]
eps50   = element.epsilon50(i);     % strain which occurs at 50% of failure

xi      = 10;                       % Factor used for determining initial stiffness

% -------------------------------------------------------------------------
% Secant stiffness, ks
% -------------------------------------------------------------------------
% According to the API standard the p-y (or p/pu-y/yc) curves are
% piece-wise linear. The limits are determined according to the API
% standard p.64.

yc      = 2.5*eps50*D;

yyslim1 = 8.0;
yyclim1 = 3.0;
yyclim2 = 15.0;

ks=NaN(2,1); %preallocation
for top_bottom=1:2 %first the top node spring, thereafter the bottom node spring
    %Initializing parameters
    pu = element.pu(i,top_bottom);          % [kN/m] Ultimate resistance, determined in layer.m
    zr = element.hr(i,top_bottom);          % [m] Transition depth
    x  = element.heqv(i,top_bottom);        % [m] Equivalent depth, according to Georgiadis' principle
    y  = y_topbottom(top_bottom);           % [m] horizontal displacement
    
%     kini = xi*pu/(D*(eps50^0.25));        % initial stiffness, other version:
    kini=(0.23*pu)/(0.1*yc);
    
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
        % Determination of secant stiffness in top of an element according to the
        % API standard p.64.
        
        if y == 0
            ks(top_bottom) = kini;
        elseif y/yc <= yyslim1
            p = min(kini*y , pu/2*(y/yc)^(1/3));
            ks(top_bottom) = p/y;
        else
            p      = pu;
            ks(top_bottom) = p/y;
        end
        
        % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')
        
        % Determination of secant stiffness in top of an element according to the
        % API standard p.64. Parts are given according to the static criterion.
        if y == 0
            ks(top_bottom) = kini;
        elseif y/yc <= yyclim1
            p = min(kini*y , pu/2*(y/yc)^(1/3));
            ks(top_bottom) = p/y;
        elseif x > zr && y/yc > yyclim1
            p      = 0.72*pu;
            ks(top_bottom) = p/y;
        elseif x <= zr && y/yc > yyclim1 && y/yc <= yyclim2
            p        = 0.72*pu*(1-(1-x/zr)*(y-yyclim1*yc)/((yyclim2-yyclim1)*yc));
            ks(top_bottom) = p/y;
        elseif x <= zr && y/yc > yyclim2
            p      = 0.72*pu*x/zr;
            ks(top_bottom) = p/y;
        end
    end
end
kstop=ks(1);
ksbot=ks(2);