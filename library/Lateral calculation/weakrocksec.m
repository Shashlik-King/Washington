function [kstop ksbot]   = weakrocksec(element,pile,loads,y_topbottom,i)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the secant spring stiffness [kN/m/m] in the top and the bottom of
% each pile segment by applying p-y curves according to the Weak Rock
% model.
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
% APPROVED        : 
%
% LAST MODIFIED   : MUOE   26.07.2013   Programming
%
%--------------------------------------------------------------------------
disp('Warning - Weak rock: Project specific change has been applied to the initial slope of the p-y curve.');
% -------------------------------------------------------------------------
% Initializing pile and soil parameters
% -------------------------------------------------------------------------
D = pile.diameter;            % Outer diameter [m]

% -------------------------------------------------------------------------
% Secant stiffness, ks
% -------------------------------------------------------------------------

ks = NaN(2,1); % Preallocation
for top_bottom = 1:2 %first the top node spring, thereafter the bottom node spring    
    %Initializing parameters
    pu = element.pu(i,top_bottom);                  % [kN/m] Ultimate resistance, determined in layer.m
    hr = abs(element.level(i,top_bottom)-element.level(i,3));    % [m] Depth below rock surface
    y  = y_topbottom(top_bottom);                   % [m] horizontal displacement
    
    % No difference between cyclic and static
    % Determination of secant stiffness in top of an element

%% Wikinger approach    
    
    E_max = element.G0(i,top_bottom)*2*(1+element.poisson(i));
%     % Determination of initial modulus
%     if hr <= 3*D
%         k_ir = 100+400*hr/(3*D);
%     else
%         k_ir = 500;
%     end
%     K_ir = k_ir*E_max; 
    
%% SNA & COU approach

K_ir=E_max;

    % Determination of boundaries for piece wise function
    y_rm = element.k_rm(i)*D;
    y_A = (pu/(2*y_rm^(1/4)*K_ir))^(4/3); 
    
    if y <= y_A
        ks(top_bottom)  = K_ir;
    elseif y > y_A 
        p               = (1/2)*pu*(y/y_rm)^(1/4);
        if p > pu
            p           = pu;
        end
        ks(top_bottom)  = p/y;
    end
end

kstop = ks(1);
ksbot = ks(2);