function [t z] = t_z(fs,pile,element,reduction,plug_unplug,A,G)
%% FORMULATION FOR t-z CURVES
% BASED ON API FORMULATION FOR CLAY
% IS USED FOR ALL COHESIVE MATERIALS 
%--------------------------------------------------------------------------
% CHANGE LOG
% 19.07.2013    MUOE - PROGRAMMING
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% fs:           is used to calculate tmax
% pile.xx:      pile.diameter is used directly in the computation of t-z curves
% element.xx:   level is used to get depth in which the t-z curves are
%               computed
%               model.tz is used to get which model of t-z curves to use
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% t:        vertical resistance
% z:        vertical displacement
%--------------------------------------------------------------------------
%% t-z formulation
%--------------------------------------------------------------------------
% values used to describe t-z curve
for i = 1:(length(element.level)-1)
    fs.G_soil(i) = G.soil(i)/(pi*(pile.diameter-2*element.thickness(i))*element.height(i));
    
    % In SACS the t-z curves need to be read from negative to positive
    % - from left to right. SACS reads positive values as compression,
    % and negative values as tension. That means we must have the t-z
    % curves defined as going from: 
    % maximum tension - tension - zero - compression - maximum compression
    
    if strcmp(element.model_axial(i),'API clay')
        %% Cohesive soils (API)
        z_peak = 0.01*pile.diameter;
        z_res = 0.02*pile.diameter; % FKMV
        z_disp = [-10 -2 -1 -0.80 -0.57 -0.31 -0.16 0 0.16 0.31 0.57 0.80 1 2 10]*z_peak;
        t_res = 0.8;
        t_tens_unplug = [-t_res -t_res -1.00 -0.90 -0.75 -0.50 -0.30]; 
        t_tens_plug = [-1.00 -1.00 -1.00 -0.90 -0.75 -0.50 -0.30];
        t_comp_unplug = [0.00 0.30 0.50 0.75 0.90 1.00 t_res t_res];
        t_comp_plug = [0.00 0.30 0.50 0.75 0.90 1.00 1.00 1.00];
        
    elseif strcmp(element.model_axial(i),'API sand')
        %% Friction soils (API)
        z_peak = 0.01*pile.diameter;
        z_disp = [-10 -2 -1 -0.80 -0.57 -0.31 -0.16 0 0.16 0.31 0.57 0.80 1 2 10]*z_peak;
        t_res = 1.0;
        t_tens_unplug = [-t_res -t_res -1.00 -0.90 -0.75 -0.50 -0.30]; 
        t_tens_plug = [-1.00 -1.00 -1.00 -0.90 -0.75 -0.50 -0.30];
        t_comp_unplug = [0.00 0.30 0.50 0.75 0.90 1.00 t_res t_res];
        t_comp_plug = [0.00 0.30 0.50 0.75 0.90 1.00 1.00 1.00];
        
    elseif strcmp(element.model_axial(i),'Zero soil')
        %% Zero soil
        z_peak = 0.01*pile.diameter;
        z_disp = [-10 -2 -1 -0.80 -0.57 -0.31 -0.16 0 0.16 0.31 0.57 0.80 1 2 10]*z_peak;
        t_tens_unplug = zeros(1,7);
        t_tens_plug = t_tens_unplug;
        t_comp_unplug = zeros(1,8);
        t_comp_plug = t_comp_unplug;
    else
        error('Please select a t-z curve that has been implemented')
    end
    
    %% Common calculations
    for j = 1:length(z_disp) % related to z_D
        z.top(i,j) = z_disp(j)*element.degradation_tz_z(i);
        z.bottom(i,j) = z.top(i,j);
    end
    for jt = 1:(length(z_disp)-1)/2 % related to tension part of z_D (not including 0.00)
        % Tension
        if strcmp(plug_unplug.tens_index,'Unplugged')
            t_tmax_tens_i = reduction.skin_tension*t_tens_unplug; 
            t.top_tens_i(i,jt) = t_tmax_tens_i(jt)*fs.i(i,1);
            t.bottom_tens_i(i,jt) = t_tmax_tens_i(jt)*fs.i(i,2);
        else % if pile is plugged in tension
            t_tmax_tens_i = t_tens_plug;
            t.top_tens_i(i,jt) = t_tmax_tens_i(jt)*fs.G_soil(i);
            t.bottom_tens_i(i,jt) = t_tmax_tens_i(jt)*fs.G_soil(i);
        end
        t_tmax_tens_o = reduction.skin_tension*t_tens_unplug;
        t.top_tens_o(i,jt) = t_tmax_tens_o(jt)*fs.o(i,1);
        t.bottom_tens_o(i,jt) = t_tmax_tens_o(jt)*fs.o(i,2);
        
        t.top_tens_w(i,jt) = t.top_tens_o(i,jt) + A.si(i)/A.so(i)*t.top_tens_i(i,jt);
        t.bottom_tens_w(i,jt) = t.bottom_tens_o(i,jt) + A.si(i)/A.so(i)*t.bottom_tens_i(i,jt);
    end     
    for jc = 1:(length(z_disp)+1)/2 % related to compression part of z_D (including 0.00)
        % Compression
        if strcmp(plug_unplug.comp_index,'Unplugged')
            t_tmax_comp_i = t_comp_unplug;
            t.top_comp_i(i,jc) = t_tmax_comp_i(jc)*fs.i(i,1);
            t.bottom_comp_i(i,jc) = t_tmax_comp_i(jc)*fs.i(i,2);
        else % if pile is plugged in compression
            t_tmax_comp_i = t_comp_plug;
            t.top_comp_i(i,jc) = t_tmax_comp_i(jc)*0;
            t.bottom_comp_i(i,jc) = t_tmax_comp_i(jc)*0;
        end
        t_tmax_comp_o = t_comp_unplug;
        t.top_comp_o(i,jc) = t_tmax_comp_o(jc)*fs.o(i,1);
        t.bottom_comp_o(i,jc) = t_tmax_comp_o(jc)*fs.o(i,2);
        
        t.top_comp_w(i,jc) = t.top_comp_o(i,jc) + A.si(i)/A.so(i)*t.top_comp_i(i,jc);
        t.bottom_comp_w(i,jc) = t.bottom_comp_o(i,jc) + A.si(i)/A.so(i)*t.bottom_comp_i(i,jc);
    end
end
    
%% Putting the tension and compression parts together
t.top = [t.top_tens_w t.top_comp_w];
t.bottom = [t.bottom_tens_w t.bottom_comp_w];