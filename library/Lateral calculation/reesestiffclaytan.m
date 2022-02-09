function [kttop ktbot]   = reesestiffclaytan(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to Reese et al. (1975)
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
% OUTPUT: kttop   : Soil stiffness at the top of the element [kN/m/m]
%         ksbot   : Soil stiffness at the bottom of the element [kN/m/m]
%
% CODE            : MMOL
% APPROVED        : CTL
%
% LAST MODIFIED   : MMOL   22.10.2009   Programming
%                   CTL    03.11.2009   Debug
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
kt=NaN(2,1);

% -------------------------------------------------------------------------
% Tangent stiffness, kt
% -------------------------------------------------------------------------

for top_bottom=1:2 %first the top node spring, thereafter the bottom node spring
    %Initializing parameters
    pu = element.pu(i,top_bottom);          % [kN/m] Ultimate resistance, determined in layer.m
    zr = element.hr(i,top_bottom);          % [m] Transition depth
    x  = element.heqv(i,top_bottom);        % [m] Equivalent depth, according to Georgiadis' principle
    y = y_topbottom(top_bottom);
    
    %Curve fitting parameters
    A   = min(0.0106*(x/D)^3-0.1006*(x/D)^2+0.3323*(x/D)+0.2042,0.6);
    B   = min(0.0229*(x/D)^3-0.1038*(x/D)^2+0.1652*(x/D)+0.202,0.3);
	yc      = eps50*D;
    yp      = 4.1*A*yc;
	k   = 5*pu*0.1^0.5/yc;%          SPSO added 15-12-2015. Initial slope up to 0.1yc
    %k=xi*pu/(D*(eps50)^0.25);      % CTL add 03.11.2009
   
    
    

    
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
        % Determination of secant stiffness in top of an element according to
        % Reese et al. 1975
        
        yyslim1=A*yc;
        yyslim2=6*A*yc;
        yyslim3=18*A*yc;
        
        if y == 0
            kt(top_bottom)  = k;      % CTL add 03.11.2009
        elseif y <= yyslim1
            kt(top_bottom)  = 0.25*(pu/((y/yc)^0.5*yc));
            kt(top_bottom)  = min(kt(top_bottom),k);        % CTL add 03.11.2009
        elseif y > yyslim1 && y <= yyslim2
            kt(top_bottom)  =0.25*(pu/((y/yc)^0.5*yc))-0.06875*pu*(((y-A*yc)/(A*yc))^0.25)/(A*yc);       % CTL modify 03.11.2009
        elseif y > yyslim2 && y <= yyslim3
            kt(top_bottom)  =  (-0.0625/yc)*pu;
        else
            kt(top_bottom)  = 0;
        end
        
        % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')
        % Determination of secant stiffness in top of an element according to
        % Reese et al. Parts are given according to the static criterion.
        yyclim0=0.1*yc;
        yyclim1=0.6*yp;
        yyclim2=1.8*yp;
        
        if y<=yyclim0% SPSO modificaion 15-12-2015
            kt(top_bottom)  = -18.4*B*pu*((0.1*yc-0.45*yp)^2/yp^2)^0.25*((0.1*yc-0.45*yp)/yp^2);
            
        elseif y<= yyclim1 && y>yyclim0      % CTL modification 03.11.2009
            kt(top_bottom)  = -18.4*B*pu*((y-0.45*yp)^2/yp^2)^0.25*((y-0.45*yp)/yp^2);
        elseif y > yyclim1 && y <= yyclim2
            kt(top_bottom)  = (-0.085/yc)*pu;
        elseif y > yyclim2
            kt(top_bottom)  = 0;
        end
    end
end

kttop=kt(1);
ktbot=kt(2);