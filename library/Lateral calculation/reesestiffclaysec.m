function [kstop ksbot]   = reesestiffclaysec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to Reese et al. (1975)
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
%
% LAST MODIFIED   : AHA    19.11.2007   Programming
%                   AHA    22.05.2008   Introduction of layered soil.
%                                       Correction line 149 putop=>pubot
%                   AHA    26.05.2008   Cyclic behaviour associated with
%                                       soft clay criterion
%
%                   MMOL   ??           Programming
%                   CTL    03.11.2009   debug
%                   JALY   29.04.2011   Streamlining code
%--------------------------------------------------------------------------
%disp('Warning - Reese stiff clay: Initial slope of p-y curve modified in accordance with recommendation in DNV-OS-J101. For some choices of soil parameters the static p-y curve will therefore not be stiffest at y=0.');
% -------------------------------------------------------------------------
% Initializing pile and soil parameters
% -------------------------------------------------------------------------

D       = pile.diameter;            % Outer diameter [m]
eps50   = element.epsilon50(i);     % strain which occurs at 50% of failure

xi=30;                              % Factor used for determining initial stiffness

%preallocation
ks=NaN(2,1);

% -------------------------------------------------------------------------
% Secant stiffness, ks
% -------------------------------------------------------------------------

for top_bottom=1:2 %first the top node spring, thereafter the bottom node spring    
    %Initializing parameters
    pu = element.pu(i,top_bottom);          % [kN/m] Ultimate resistance, determined in layer.m
    zr = element.hr(i,top_bottom);          % [m] Transition depth
    x  = element.heqv(i,top_bottom);        % [m] Equivalent depth, according to Georgiadis' principle
    y  = y_topbottom(top_bottom);           % [m] horizontal displacement
    
    %Curve fitting parameters
	A  = min(0.0106*(x/D)^3-0.1006*(x/D)^2+0.3323*(x/D)+0.2042,0.6); 
    B   = min(0.0229*(x/D)^3-0.1038*(x/D)^2+0.1652*(x/D)+0.202,0.3);
	yc      = eps50*D;
    yp      = 4.1*A*yc;
    kini   = 5*pu*0.1^0.5/yc;%          SPSO added 15-12-2015. Initial slope up to 0.1yc
    %kini   = xi*pu/(D*(eps50^0.25));
    
    
    
    
    
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
    % Determination of secant stiffness in top of an element according to
    % Reese et al. 1975
    yyslim1=A*yc;
    yyslim2=6*A*yc;
    yyslim3=18*A*yc;
        
        if y <= yyslim1
            if y == 0
                ks(top_bottom) = kini;   % CTL modification 03.11.2009
            else
                p1     = kini*y;
                p2     = 0.5*pu*(y/yc)^0.5;
                p      = min(p1,p2);
                ks(top_bottom)  = p/y;
            end
        elseif y > yyslim1 && y <= yyslim2
            p      = 0.5*pu*(y/yc)^0.5-0.055*pu*((y-A*yc)/(A*yc))^1.25;
            ks(top_bottom)  = p/y;
        elseif y > yyslim2 && y <= yyslim3
            p      = 0.5*pu*(6*A)^0.5-0.411*pu-(0.0625/yc)*pu*(y-6*A*yc);
            ks(top_bottom)  = p/y;
        else
            p      = 0.5*pu*(6*A)^0.5-0.411*pu-0.75*pu*A;
            ks(top_bottom)  = p/y;
        end
    
    % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')
    % Determination of secant stiffness in top of an element according to
    % Reese et al. Parts are given according to the static criterion.
    yyclim1=0.6*yp;
    yyclim2=1.8*yp;
        if y <= yyclim1
            if y <= 0.1*yc
                pyc=B*pu*(1-(abs((0.1*yc-0.45*yp)/(0.45*yp)))^2.5);%SPSO modified 15-12-2015
                ks(top_bottom)  = pyc/(0.1*yc);%SPSO added 16-12-2015
            else
                p     = B*pu*(1-(abs((y-0.45*yp)/(0.45*yp)))^2.5);%SPSO modified 15-12-2015
                ks(top_bottom)  = p/y;  % MMOL&CTL modification 03.11.2009
            end
        elseif y > yyclim1 && y <= yyclim2
            p      = 0.936*B*pu-(0.085/yc)*pu*(y-0.6*yp);
            ks(top_bottom)  = p/y;   % CTL modification 03.11.2009
        else
            p      = 0.936*B*pu-(0.102/yc)*pu*yp;
            ks(top_bottom)  = p/y;   % CTL modification 03.11.2009
        end
    end
end

kstop=ks(1);
ksbot=ks(2);