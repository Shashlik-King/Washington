function [torsion_stiffness] = torsion_stiffness(pile,element,t)
%% Estimation of torsional pile stiffness
% 
%--------------------------------------------------------------------------
% CHANGE LOG
% 24.08.2015    SPSO - PROGRAMMING
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% t:            top_comp_o and bottom_comp_o are used directly
% pile.xx:      pile.diameter is used directly in the computation of the
%               rotational stiffness
% element.xx:   level is used to get depth in which the t-z curves are
%               computed
%               model.tz is used to get which model of t-z curves to use
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% torsion:      Torsional soil-pile stiffness in kNm/rad
%--------------------------------------------------------------------------
%% Calculation
%--------------------------------------------------------------------------
z_max=0.01*pile.diameter;%m
rad_max=z_max/(pile.diameter/2);%radians
for i=1:(size(element.level,1)-1)
    t.max_average(i)=(max(t.top_comp_o(i,:))+max(t.bottom_comp_o(i,:)))/2;%kN/m^2
    torsional_moment_max(i)=t.max_average(i)*pile.diameter/2*pi*pile.diameter*(element.level(i,1)-element.level(i,2));%kNm
end
torsion_stiffness=sum(torsional_moment_max)/rad_max;%kNm/rad