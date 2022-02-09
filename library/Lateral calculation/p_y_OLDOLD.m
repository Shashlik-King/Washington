function [p y y_tot output toe_plot_u] = p_y(node,pile,element,loads,output,settings)
%% FORMULATION FOR p-y CURVES
%--------------------------------------------------------------------------
% CHANGE LOG
% 01.08.2013    MUOE - SEPARATING FROM WINKLER_PY.M
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% pile.xx:      pile.diameter is used to determine yc
% element.xx:   model.py is used to determine which p-y curve to use
%               epsilon50 is used to determine yc
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% p:        lateral resistance
% y:        lateral displacement
%--------------------------------------------------------------------------
%% Initial
%--------------------------------------------------------------------------
nelem = length(element.model_py);
y_tot = zeros(nelem,2);
%--------------------------------------------------------------------------
%% Calculation routine
%--------------------------------------------------------------------------
for c = 1:(nelem-1)
    if strcmp(element.model_py(c),'Zero soil')
        upy_top(c,:) = zeros(1,17);
        upy_bot(c,:) = upy_top(c,:);
        
    elseif strcmp(element.model_py(c),'API sand') || strcmp(element.model_py(c),'Kirsch sand') || strcmp(element.model_py(c),'Kallehave sand')|| strcmp(element.model_py(c),'PISA sand')|| strcmp(element.model_py(c),'PISA clay') || strcmp(element.model_py(c),'PISA Bothkennar clay') || strcmp(element.model_py(c),'PISA Cowden clay') || strcmp(element.model_py(c),'PISA Upper clay') || strcmp(element.model_py(c),'PISA Till') || strcmp(element.model_py(c),'PISA Chalk')
%         yc = 0.10;
%         upy_top(c,:) = [0 0.02 0.05 0.10 0.15 0.2 0.25 0.30 0.37 0.45 0.5 0.65 0.85 1 1.5]*yc;
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4 0.6 0.8 1];
        upy_bot(c,:) = upy_top(c,:);
%         upy_top(c,:) = [0 0.00000005 0.0001 0.0003 0.0005 0.00075 0.001 0.002 0.003 0.004 0.005 0.0075 0.01 0.01125 0.0125 0.01375 0.015 0.0175 0.02 0.025 0.03 0.035 0.04 0.05 0.06 0.07 0.080 0.1 0.12 0.16 0.2 0.3 0.4];
%         upy_bot(c,:) = upy_top(c,:);
    elseif element.PISA_switch==1 && strcmp(element.type{c}, 'Sand') 
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4 0.6 0.8 1];
        upy_bot(c,:) = upy_top(c,:);    
        
    elseif element.PISA_switch==1 && strcmp(element.type{c}, 'Clay') 
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4 0.6 0.8 1];
        upy_bot(c,:) = upy_top(c,:);   
        
    elseif strcmp(element.model_py(c),'API clay') || strcmp(element.model_py(c),'Stiff clay w/o free water') || strcmp(element.model_py(c),'Kirsch soft clay') || strcmp(element.model_py(c),'Kirsch stiff clay')
        yc_top(c) = 2.5*element.epsilon50(c,1)*element.diameter(c);
        yc_bot(c) = 2.5*element.epsilon50(c,2)*element.diameter(c);
%         if strcmp(element.model_py(c),'API clay') || strcmp(element.model_py(c),'Kirsch soft clay')
%             ylimit_top(c) = 15*yc_top(c);
%             ylimit_bot(c) = 15*yc_bot(c);
%             upy_top(c,:) = [0 0.00667 0.02 0.04 0.08 0.1 0.12 0.15 0.20 0.30 0.40 0.50 0.70 1.00 1.50]*ylimit_top(c);
%             upy_bot(c,:) = [0 0.00667 0.02 0.04 0.08 0.1 0.12 0.15 0.20 0.30 0.40 0.50 0.70 1.00 1.50]*ylimit_bot(c);
%         elseif strcmp(element.model_py(c),'Stiff clay w/o free water') || strcmp(element.model_py(c),'Kirsch stiff clay')
%             ylimit_top(c) = 16*yc_top(c);
%             ylimit_bot(c) = 16*yc_bot(c);
%             upy_top(c,:) = [0 0.00625 0.02 0.04 0.08 0.10 0.12 0.15 0.20 0.30 0.40 0.5 0.70 1.00 1.20]*ylimit_top(c);
%             upy_bot(c,:) = [0 0.00625 0.02 0.04 0.08 0.10 0.12 0.15 0.20 0.30 0.40 0.5 0.70 1.00 1.20]*ylimit_bot(c);
%         end
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4 0.6 0.8 1];
        upy_bot(c,:) = upy_top(c,:);
%         upy_top(c,:) = [0 0.00000005 0.0001 0.0003 0.0005 0.00075 0.001 0.002 0.003 0.004 0.005 0.0075 0.01 0.01125 0.0125 0.01375 0.015 0.0175 0.02 0.025 0.03 0.035 0.04 0.05 0.06 0.07 0.080 0.1 0.12 0.16 0.2 0.3 0.4];
%         upy_bot(c,:) = upy_top(c,:);
        y_tot(c,1) = yc_top(c);
        y_tot(c,2) = yc_bot(c);
        
    elseif strcmp(element.model_py(c),'Reese stiff clay')
        A_top = min(0.0106*(element.heqv(c,1)/element.diameter(c))^3-0.1006*(element.heqv(c,1)/element.diameter(c))^2+0.3323*(element.heqv(c,1)/element.diameter(c))+0.2042,0.6) ;
        A_bot = min(0.0106*(element.heqv(c,2)/element.diameter(c))^3-0.1006*(element.heqv(c,2)/element.diameter(c))^2+0.3323*(element.heqv(c,2)/element.diameter(c))+0.2042,0.6);
        yc_top(c) = element.epsilon50(c,1)*element.diameter(c);
        yc_bot(c) = element.epsilon50(c,2)*element.diameter(c);
        yp_top(c) = 4.1*A_top*yc_top(c);
        yp_bot(c) = 4.1*A_bot*yc_bot(c);
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4];
        upy_bot(c,:) = upy_top(c,:);
%         upy_top(c,:) = [0 0.00000005 0.0001 0.0003 0.0005 0.00075 0.001 0.002 0.003 0.004 0.005 0.0075 0.01 0.01125 0.0125 0.01375 0.015 0.0175 0.02 0.025 0.03 0.035 0.04 0.05 0.06 0.07 0.080 0.1 0.12 0.16 0.2 0.3 0.4];
%         upy_bot(c,:) = upy_top(c,:);
        y_tot(c,1) = yp_top(c);
        y_tot(c,2) = yp_bot(c);
        
    elseif strcmp(element.model_py(c),'Modified Weak rock')
        yrm(c) = element.k_rm(c)*element.diameter(c);
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4];
        upy_bot(c,:) = upy_top(c,:);
%         upy_top(c,:) = [0 0.00000005 0.0001 0.0003 0.0005 0.00075 0.001 0.002 0.003 0.004 0.005 0.0075 0.01 0.01125 0.0125 0.01375 0.015 0.0175 0.02 0.025 0.03 0.035 0.04 0.05 0.06 0.07 0.080 0.1 0.12 0.16 0.2 0.3 0.4];
%         upy_bot(c,:) = upy_top(c,:);
        y_tot(c,1) = yrm(c);
        y_tot(c,2) = yrm(c);
    
    end
    
    npoints = size(upy_top,2);
    
    for j = 1:npoints
        [ksppynode_top ksppynode_bot] = secspringstiff(element,pile,loads,[upy_top(c,j) upy_bot(c,j)],c);
        ppynode_top(c,j) = ksppynode_top*upy_top(c,j);  
        ppynode_bot(c,j) = ksppynode_bot*upy_bot(c,j);  
    end
    
    if strcmp(element.model_py(c),'Zero soil')
        upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4];
        upy_bot(c,:) = upy_top(c,:);
    end
    
         %% Toe shear
         for g = nelem % the last element, toe shear
             ybot = linspace(0.0,0.1,100);
           for  j = 1:length(ybot)
                 %
                 if settings.toe_shear == 1
    %                 output.ts = toe_shear(element,node,pile,ybot,c);
               [kHb_top(j) kHb_bot(j)] = secHBstiff(element,pile,[ybot(j) ybot(j)],g); 
                Hb_bot(j) = kHb_bot(j)*ybot(j); 
                output.ts = Hb_bot*(element.level(end,1)-element.level(end,2));
                elseif settings.toe_shear == 0
                     output.ts = 0;
                end
           end
			 toe_plot_u = ybot;
         end
end

% Putting p and y data into same format as t and z
for i = 1:size(upy_top,1)
    p.top(i,:)      = ppynode_top(i,:);
    p.bottom(i,:)   = ppynode_bot(i,:);
    y.top(i,:)      = upy_top(i,:);
    y.bottom(i,:)   = upy_bot(i,:);
end