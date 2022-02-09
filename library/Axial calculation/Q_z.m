function [Q zQ] = Q_z(qp,pile,element)
%% FORMULATION FOR Q-z CURVES
% IS USED FOR ALL MATERIALS 
%--------------------------------------------------------------------------
% CHANGE LOG
% 29.07.2014    MUOE - PROGRAMMING
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% qp:           is used to calculate Qmax
% pile.xx:      pile.diameter is used directly in the computation of t-z curves
% element.xx:   level is used to get depth in which the t-z curves are
%               computed
%               model.tz is used to get which model of t-z curves to use
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% Q:        end bearing resistance
% z:        vertical displacement
%--------------------------------------------------------------------------
%% Q-z formulation
%--------------------------------------------------------------------------
% values used to describe Q-z curve
for i = 1:length(element.level)-1
    
    if strcmp(element.model_axial(i),'API clay')
        %% Cohesive soils (API)
        z_disp = [0.000 0.002 0.013 0.042 0.073 0.100 0.200]*pile.diameter;
        Qp = [0.00 0.25 0.50 0.75 0.90 1.00 1.00];
        
    elseif strcmp(element.model_axial(i),'API sand')
        %% Friction soils (API)
        z_disp = [0.000 0.002 0.013 0.042 0.073 0.100 0.200]*pile.diameter;
        Qp = [0.00 0.25 0.50 0.75 0.90 1.00 1.00];
        
    elseif strcmp(element.model_axial(i),'Zero soil')
        %% Zero soil
        z_disp = [0.000 0.002 0.013 0.042 0.073 0.100 0.200]*pile.diameter;
        Qp = zeros(1,7);
    else
        error('Please select a Q-z curve that has been implemented')
    end
    
    %% Common calculations
    for j = 1:length(z_disp) % related to z_D
        zQ.top(i,j) = z_disp(j);
        zQ.bottom(i,j) = zQ.top(i,j);
        
        Q.top(i,j) = Qp(j)*qp(i,1);
        Q.bottom(i,j) = Qp(j)*qp(i,2);
    end
end