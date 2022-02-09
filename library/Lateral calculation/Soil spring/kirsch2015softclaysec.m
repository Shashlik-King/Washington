function [kstop ksbot] = kirschsoftclaysec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves for soft clay according to Kirsch (2014).
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
    % Initializing parameters
    pu = element.pu(i,top_bottom); % [kN/m] Ultimate resistance, determined in layer.m
    zr = element.hr(i,top_bottom); % [m] Transition depth
    x  = element.heqv(i,top_bottom); % [m] Equivalent depth, according to Georgiadis' principle
    y  = y_topbottom(top_bottom); % [m] horizontal displacement
    Es = element.Es(i)/1E3; % [MPa] Large-strain constrained modulus 
    Esd =10^(-0.42*log10(0.0006*Es))*Es; % [MPa] Small-strain constrained modulus
    yyslim1 = 8.0;
    yyclim1 = 3.0;
    yyclim2 = 15.0;
    
    tol = 0.001*eps50; %SPSO modified
    
    e50mod = eps50;
    yc = 2.5*e50mod*D;
    
    % -------- Static criterion -----------------------------------------------
    if strcmp(loads.static_cyclic,'static')
        
        % iteration to find eps50mod
        eps50mod = min(eps50*(1+(1-0.50*(y/yc)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
        while abs(e50mod-eps50mod) > tol && abs(eps50mod-eps50)>0.001*eps50 %SPSO added part from &&
            e50mod = (e50mod+eps50mod)/2;
            ycmod = 2.5*e50mod*D;
            eps50mod = min(eps50*(1+(1-0.50*(y/ycmod)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
        end  
        ycmod = 2.5*eps50mod*D;
        
        %iteration to find yini
        yini=0.1*ycmod;
        delta_yini=abs(y-yini);
        while delta_yini>0.0001
            eps50mod_ini = min(eps50*(1+(1-0.50*(yini/ycmod)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
            while abs(e50mod-eps50mod_ini) > tol
                e50mod = (e50mod+eps50mod_ini)/2;
                ycmod_ini = 2.5*e50mod*D;
                eps50mod_ini = min(eps50*(1+(1-0.50*(yini/ycmod_ini)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
            end  
            ycmod_ini = 2.5*eps50mod_ini*D;
            delta_yini=abs(yini-0.1*ycmod_ini);
            yini=0.1*ycmod_ini;
        end
        
        kini = (0.23208*pu)/(yini); %SPSO modified
        
        if y == 0
            ks(top_bottom) = kini;
        elseif y/ycmod <= yyslim1
            p = min(kini*y , pu/2*(y/ycmod)^(1/3));
            ks(top_bottom) = p/y;
        else
            p      = pu;
            ks(top_bottom) = p/y;
        end
        
    % -------- Cyclic criterion -----------------------------------------------
    elseif strcmp(loads.static_cyclic,'cyclic')   
        
        % iteration to find eps50mod
        eps50mod = min(eps50*(1+(1-0.50/0.72*(y/yc)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
        while abs(e50mod-eps50mod) > tol && abs(eps50mod-eps50)>0.001*eps50 %SPSO added part from &&
            e50mod = (e50mod+eps50mod)/2;
            ycmod = 2.5*e50mod*D;
            eps50mod = min(eps50*(1+(1-0.50/0.72*(y/ycmod)^(1/3))*(Es/Esd-1)),eps50); %SPSO modified
        end  
        ycmod = 2.5*eps50mod*D;
        
        %iteration to find yini
        yini=0.1*ycmod; %new
        delta_yini=abs(y-yini);%new
        while delta_yini>0.0001%new
            eps50mod_ini = min(eps50*(1+(1-0.50/0.72*(yini/ycmod)^(1/3))*(Es/Esd-1)),eps50); %new
            while abs(e50mod-eps50mod_ini) > tol%new
                e50mod = (e50mod+eps50mod_ini)/2;%new
                ycmod_ini = 2.5*e50mod*D;%new
                eps50mod_ini = min(eps50*(1+(1-0.50/0.72*(yini/ycmod_ini)^(1/3))*(Es/Esd-1)),eps50); %new
            end  %new
            ycmod_ini = 2.5*eps50mod_ini*D;%new
            delta_yini=abs(yini-0.1*ycmod_ini);%new
            yini=0.1*ycmod_ini;%new
        end%new
        
        kini = (0.23208*pu)/(yini); %SPSO modified
            
        if y == 0
            ks(top_bottom) = kini;
        elseif y/ycmod <= yyclim1
            p = min(kini*y , pu/2*(y/ycmod)^(1/3));
            ks(top_bottom) = p/y;
        elseif x > zr && y/ycmod > yyclim1
            p      = 0.72*pu;
            ks(top_bottom) = p/y;
        elseif x <= zr && y/ycmod > yyclim1 && y/ycmod <= yyclim2
            p        = 0.72*pu*(1-(1-x/zr)*(y-yyclim1*ycmod)/((yyclim2-yyclim1)*ycmod));
            ks(top_bottom) = p/y;
        elseif x <= zr && y/ycmod > yyclim2
            p      = 0.72*pu*x/zr;
            ks(top_bottom) = p/y;
        end
    end
end
kstop=ks(1);
ksbot=ks(2);