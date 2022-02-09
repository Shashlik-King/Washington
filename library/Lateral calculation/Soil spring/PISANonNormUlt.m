function [punon munon ptop pbot] = PISANonNormUlt(element,pile,y_topbottom,ks,c)
%--------------------------------------------------------------------------
% PURPOSE
% Compute non-normalized ultimated parameters pu_non, mu_non which are needed for
% UR plots. 
%
% INPUT:  D   : pile outer diameter
%         L   : pile length
%         nelem  : number of elements
%
% OUTPUT: pu_non      : Normalized ultimate horisontal capacity at the top and the bottom of the
%                   segment [-]  [pu_top pu_bot]. Size=(nelem,2)
%         mu_non      : Normalized ultimate moment capacity at the top and the bottom of the

%
% Log:
% 17.08.2016    EVVA    Setting up the function
%--------------------------------------------------------------------------

nelem=length(element.model_py);

%Preallocation
element.pu_non=NaN(nelem,2);
element.mu_non=NaN(nelem,2);

% Input
D = pile.diameter;          % Outer diameter [m]
L = pile.L;      % Pile length [m]


su_top   = element.cu(c,1);     	% Undrained shear strength in the top node[kPa]
su_bot   = element.cu(c,2);         % Undrained shear strength in the bottom node[kPa]
ytop    = abs(y_topbottom(1));      % Horizontal disp. at the top of the pile element [m]
ybot    = abs(y_topbottom(2));      % Horizontal disp. at the bottom of the pile element [m]
kstop   = ks(1);
ksbot   = ks(2);
sigma_v_top = element.sigma_v_eff(c,1);
sigma_v_bot = element.sigma_v_eff(c,2);
pu_n_top = element.pu(c,1);         % Ultimate distributed load at the top of the p.el. (normalized) [-]
pu_n_bot = element.pu(c,2);         % Ultimate distributed load at the bot of the p.el. (normalized) [-]
mu_n_top = element.mu(c,1);         % Ultimate distr. moment at the top of the p.el. (normalized) [-]
mu_n_bot = element.mu(c,2);         % Ultimate distr. moment at the bot of the p.el. (normalized) [-]

         %% -------------------------------------------------------------------------
         % Free standing length, zero soil
         % -------------------------------------------------------------------------
     if strcmp(element.model_py(c),'Zero soil')
            pu_non_top = 0;          % Ult. dist. load at the top of the p.el.[kN/m]
            pu_non_bot = 0;          % Ult. dist. load at the top of the p.el.[kN/m]
            ptop = 0;                % Distributed load  at the top of the pile element [m]                        
            pbot = 0;                % Distributed load  at the bot of the pile element [m]  
            mu_non_top = 0;          % Ult. dist. load at the top of the p.el.[kNm/m]
            mu_non_bot = 0;          % Ult. dist. load at the top of the p.el.[kNm/m]   
		%% -------------------------------------------------------------------------
        % PISA sand
        % --------------------------------------------------------------------------
     elseif strcmp(element.type{c},'Sand')
 %for top_bottom = 1:2
            % -------- Determination of non-normalized ultimate values
           % pu_non = element.pu(c,top_bottom)*element.sigma_v_eff(c,top_bottom)*D;   % Ult. dist. load at the top of the p.el.[kN/m]
            pu_non_top = pu_n_top.*sigma_v_top*D;   % Ult. dist. load at the top of the p.el.[kN/m]
            pu_non_bot = pu_n_bot.*sigma_v_bot*D;   % Ult. dist. load at the top of the p.el.[kN/m]
            ptop = kstop.*ytop;                     % Distributed load  at the top of the pile element [m]                        
            pbot = ksbot.*ybot;                     % Distributed load  at the bot of the pile element [m]  
            mu_non_top = mu_n_top.*ptop*D;          % Ult. dist. load at the top of the p.el.[kNm/m]
            mu_non_bot = mu_n_bot.*pbot*D;          % Ult. dist. load at the top of the p.el.[kNm/m]
           % mu_non = element.mu(c,top_bottom)*element.sigma_v_eff(c,top_bottom)*D; % Ult. dist. load at the top of the p.el.[kNm/m]
 %end		
 		%% -------------------------------------------------------------------------
        % PISA clay
        % --------------------------------------------------------------------------
     elseif strcmp(element.type{c},'Clay')
            pu_non_top = pu_n_top.*su_top*D;   % Ult. dist. load at the top of the p.el.[kN/m]
            pu_non_bot = pu_n_bot.*su_bot*D;   % Ult. dist. load at the top of the p.el.[kN/m]
            ptop = kstop.*ytop;                % Distributed load  at the top of the pile element [m]                        
            pbot = ksbot.*ybot;                % Distributed load  at the bot of the pile element [m]  
            mu_non_top = mu_n_top.*su_top*D^2; % Ult. dist. load at the top of the p.el.[kNm/m]
            mu_non_bot = mu_n_bot.*su_bot*D^2; % Ult. dist. load at the top of the p.el.[kNm/m]
       
    
        %% -------------------------------------------------------------------------
        % Not recognised soil model
        % -------------------------------------------------------------------------
    else
        error('The specified soil model is not supported in PISANonNormUlt.m')
     end   
    punon=[pu_non_top pu_non_bot];
    munon=[mu_non_top mu_non_bot];
    ptop;
    pbot;
end