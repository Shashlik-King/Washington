% Function that combines all PISA spring formulations for respective layers
% into one function to decrease QA time and chances of errors. 
%
% Units: kN, m, s, kPa
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Log of changes------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Date            Initials        Change
%2020.03.25      FKMV            Programming
%2020.08.18      FKMV            Correction of sand formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAM 1 = Normalized ultimate lateral displacement for p-y

function [element]=PISA_formulations_calculation(pile,~,element,~,~,~,i)
% (pile,soil,element,scour,settings,data,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P-Y springs

% Normalized ultimate lateral displacement
element.PISA_param(i,1) = element.PISA_prelim_param.p_y(i,2);

for j = 1:2 % loop over top and bottom of element

% Ultimate lateral load
if strcmp(element.type{i},'Clay')
    if element.PISA_prelim_param.p_y(i,5) == 0
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6);
    elseif element.PISA_prelim_param.p_y(i,5) == 1
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.p_y(i,7);
    elseif element.PISA_prelim_param.p_y(i,5) == 2
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6)+element.PISA_prelim_param.p_y(i,7)*exp(element.PISA_prelim_param.p_y(i,8)*(-element.level(i,j)/pile.diameter));
    else
        error('wrong type of function defined for a spring')
    end
elseif strcmp(element.type{i},'Sand')
    if element.PISA_prelim_param.p_y(i,5) == 0
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6);
    elseif element.PISA_prelim_param.p_y(i,5) == 1
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6)*(-element.level(i,j)/pile.L)+element.PISA_prelim_param.p_y(i,7);
    elseif element.PISA_prelim_param.p_y(i,5) == 2
        element.PISA_param(i,1+j) = element.PISA_prelim_param.p_y(i,6)+element.PISA_prelim_param.p_y(i,7)*exp(element.PISA_prelim_param.p_y(i,8)*(-element.level(i,j)/pile.diameter));
    else
        error('wrong type of function defined for a spring')
    end
end															 

% Normalized initial stiffness
if element.PISA_prelim_param.p_y(i,9) == 0
    element.PISA_param(i,3+j) = element.PISA_prelim_param.p_y(i,10);
elseif element.PISA_prelim_param.p_y(i,9) == 1
    element.PISA_param(i,3+j) = element.PISA_prelim_param.p_y(i,10)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.p_y(i,11);
elseif element.PISA_prelim_param.p_y(i,9) == 2
    element.PISA_param(i,3+j) = element.PISA_prelim_param.p_y(i,10)+element.PISA_prelim_param.p_y(i,11)*exp(element.PISA_prelim_param.p_y(i,12)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,3+j) = max(element.PISA_param(i,3+j) , element.PISA_param(i,1+j)/element.PISA_param(i,1));

% Curvature
if element.PISA_prelim_param.p_y(i,13) == 0
    element.PISA_param(i,5+j) = element.PISA_prelim_param.p_y(i,14);
elseif element.PISA_prelim_param.p_y(i,13) == 1
    element.PISA_param(i,5+j) = element.PISA_prelim_param.p_y(i,14)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.p_y(i,15);
elseif element.PISA_prelim_param.p_y(i,13) == 2
    element.PISA_param(i,5+j) = element.PISA_prelim_param.p_y(i,14)+element.PISA_prelim_param.p_y(i,15)*exp(element.PISA_prelim_param.p_y(i,16)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,5+j) = max(element.PISA_param(i,5+j),0);
element.PISA_param(i,5+j) = min(element.PISA_param(i,5+j),0.999);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%																								  
%% M-T springs

% Normalized ultimate lateral displacement
element.PISA_param(i,8) = element.PISA_prelim_param.m_t(i,2);

for j = 1:2 % loop over top and bottom of element

% Ultimate lateral load
if strcmp(element.type{i},'Clay')
    if element.PISA_prelim_param.m_t(i,5) == 0
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6);
    elseif element.PISA_prelim_param.m_t(i,5) == 1
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.m_t(i,7);
    elseif element.PISA_prelim_param.m_t(i,5) == 2
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6)+element.PISA_prelim_param.m_t(i,7)*exp(element.PISA_prelim_param.m_t(i,8)*(-element.level(i,j)/pile.diameter));
    else
        error('wrong type of function defined for a spring')
    end
elseif strcmp(element.type{i},'Sand')
    if element.PISA_prelim_param.m_t(i,5) == 0
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6);
    elseif element.PISA_prelim_param.m_t(i,5) == 1
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6)*(-element.level(i,j)/pile.L)+element.PISA_prelim_param.m_t(i,7);
    elseif element.PISA_prelim_param.m_t(i,5) == 2
        element.PISA_param(i,8+j) = element.PISA_prelim_param.m_t(i,6)+element.PISA_prelim_param.m_t(i,7)*exp(element.PISA_prelim_param.m_t(i,8)*(-element.level(i,j)/pile.diameter));
    else
        error('wrong type of function defined for a spring')
    end
end

% element.PISA_param(i,8+j) = max(element.PISA_param(i,8+j) , 0.03);

% Normalized initial stiffness
if element.PISA_prelim_param.m_t(i,9) == 0
    element.PISA_param(i,10+j) = element.PISA_prelim_param.m_t(i,10);
elseif element.PISA_prelim_param.m_t(i,9) == 1
    element.PISA_param(i,10+j) = element.PISA_prelim_param.m_t(i,10)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.m_t(i,11);
elseif element.PISA_prelim_param.m_t(i,9) == 2
    element.PISA_param(i,10+j) = element.PISA_prelim_param.m_t(i,10)+element.PISA_prelim_param.m_t(i,11)*exp(element.PISA_prelim_param.m_t(i,12)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

%element.PISA_param(i,10+j) = max(element.PISA_param(i,10+j) , element.PISA_param(i,8+j)/element.PISA_param(i,8));

% Normalized ultimate rotation 
if strcmp(element.type(i),'Clay')
    element.PISA_param(i,8) = element.PISA_prelim_param.m_t(i,2);
elseif strcmp(element.type(i),'Sand')
    element.PISA_param(i,8) = element.PISA_param(i,10+j)/element.PISA_param(i,8+j);
end
% Curvature
if element.PISA_prelim_param.m_t(i,13) == 0
    element.PISA_param(i,12+j) = element.PISA_prelim_param.m_t(i,14);
elseif element.PISA_prelim_param.m_t(i,13) == 1
    element.PISA_param(i,12+j) = element.PISA_prelim_param.m_t(i,14)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.m_t(i,15);
elseif element.PISA_prelim_param.m_t(i,13) == 2
    element.PISA_param(i,12+j) = element.PISA_prelim_param.m_t(i,14)+element.PISA_prelim_param.m_t(i,15)*exp(element.PISA_prelim_param.m_t(i,16)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,12+j) = max(element.PISA_param(i,12+j),0);
element.PISA_param(i,12+j) = min(element.PISA_param(i,12+j),0.999);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P-Y toe springs

% Normalized ultimate lateral displacement
if element.PISA_prelim_param.Hb(i,1) == 0
    element.PISA_param(i,15) = element.PISA_prelim_param.Hb(i,2);
elseif element.PISA_prelim_param.Hb(i,1) == 1
    element.PISA_param(i,15) = element.PISA_prelim_param.Hb(i,2)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Hb(i,3);
elseif element.PISA_prelim_param.Hb(i,1) == 2
    element.PISA_param(i,15) = element.PISA_prelim_param.Hb(i,2)+element.PISA_prelim_param.Hb(i,3)*exp(element.PISA_prelim_param.Hb(i,4)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end
element.PISA_param(i,15) = max(element.PISA_param(i,15),0.001);															   

% Ultimate lateral load
if element.PISA_prelim_param.Hb(i,5) == 0
    element.PISA_param(i,16) = element.PISA_prelim_param.Hb(i,6);
elseif element.PISA_prelim_param.Hb(i,5) == 1
    element.PISA_param(i,16) = element.PISA_prelim_param.Hb(i,6)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Hb(i,7);
elseif element.PISA_prelim_param.Hb(i,5) == 2
    element.PISA_param(i,16) = element.PISA_prelim_param.Hb(i,6)+element.PISA_prelim_param.Hb(i,7)*exp(element.PISA_prelim_param.Hb(i,8)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end
element.PISA_param(i,16) = max(element.PISA_param(i,16),0.001);															   

% Normalized initial stiffness
if element.PISA_prelim_param.Hb(i,9) == 0
    element.PISA_param(i,17) = element.PISA_prelim_param.Hb(i,10);
elseif element.PISA_prelim_param.Hb(i,9) == 1
    element.PISA_param(i,17) = element.PISA_prelim_param.Hb(i,10)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Hb(i,11);
elseif element.PISA_prelim_param.Hb(i,9) == 2
    element.PISA_param(i,17) = element.PISA_prelim_param.Hb(i,10)+element.PISA_prelim_param.Hb(i,11)*exp(element.PISA_prelim_param.Hb(i,12)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end
element.PISA_param(i,17) = max(element.PISA_param(i,17),0.001);															   

% Curvature
if element.PISA_prelim_param.Hb(i,13) == 0
    element.PISA_param(i,18) = element.PISA_prelim_param.Hb(i,14);
elseif element.PISA_prelim_param.Hb(i,13) == 1
    element.PISA_param(i,18) = element.PISA_prelim_param.Hb(i,14)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Hb(i,15);
elseif element.PISA_prelim_param.Hb(i,13) == 2
    element.PISA_param(i,18) = element.PISA_prelim_param.Hb(i,14)+element.PISA_prelim_param.Hb(i,15)*exp(element.PISA_prelim_param.Hb(i,16)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,18) = max(element.PISA_param(i,18),0);
element.PISA_param(i,18) = min(element.PISA_param(i,18),0.999);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M-T toe springs

% Normalized ultimate lateral displacement
if element.PISA_prelim_param.Mb(i,1) == 0
    element.PISA_param(i,19) = element.PISA_prelim_param.Mb(i,2);
elseif element.PISA_prelim_param.Mb(i,5) == 1
    element.PISA_param(i,19) = element.PISA_prelim_param.Mb(i,2)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Mb(i,3);
elseif element.PISA_prelim_param.Mb(i,5) == 2
    element.PISA_param(i,19) = element.PISA_prelim_param.Mb(i,2)+element.PISA_prelim_param.Mb(i,3)*exp(element.PISA_prelim_param.Mb(i,4)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end														   
element.PISA_param(i,19) = max(element.PISA_param(i,19),0.001);															   

% Ultimate lateral load
if element.PISA_prelim_param.Mb(i,5) == 0
    element.PISA_param(i,20) = element.PISA_prelim_param.Mb(i,6);
elseif element.PISA_prelim_param.Mb(i,5) == 1
    element.PISA_param(i,20) = element.PISA_prelim_param.Mb(i,6)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Mb(i,7);
elseif element.PISA_prelim_param.Mb(i,5) == 2
    element.PISA_param(i,20) = element.PISA_prelim_param.Mb(i,6)+element.PISA_prelim_param.Mb(i,7)*exp(element.PISA_prelim_param.Mb(i,8)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end
element.PISA_param(i,20) = max(element.PISA_param(i,20),0.001);															   

% Normalized initial stiffness
if element.PISA_prelim_param.Mb(i,9) == 0
    element.PISA_param(i,21) = element.PISA_prelim_param.Mb(i,10);
elseif element.PISA_prelim_param.Mb(i,9) == 1
    element.PISA_param(i,21) = element.PISA_prelim_param.Mb(i,10)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Mb(i,11);
elseif element.PISA_prelim_param.Mb(i,9) == 2
    element.PISA_param(i,21) = element.PISA_prelim_param.Mb(i,10)+element.PISA_prelim_param.Mb(i,11)*exp(element.PISA_prelim_param.Mb(i,12)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,21) = max(element.PISA_param(i,21) , element.PISA_param(i,20)/element.PISA_param(i,19));																											 
% Curvature
if element.PISA_prelim_param.Mb(i,13) == 0
    element.PISA_param(i,22) = element.PISA_prelim_param.Mb(i,14);
elseif element.PISA_prelim_param.Mb(i,13) == 1
    element.PISA_param(i,22) = element.PISA_prelim_param.Mb(i,14)*(-element.level(i,j)/pile.diameter)+element.PISA_prelim_param.Mb(i,15);
elseif element.PISA_prelim_param.Mb(i,13) == 2
    element.PISA_param(i,22) = element.PISA_prelim_param.Mb(i,14)+element.PISA_prelim_param.Mb(i,15)*exp(element.PISA_prelim_param.Mb(i,16)*(-element.level(i,j)/pile.diameter));
else
    error('wrong type of function defined for a spring')
end

element.PISA_param(i,22) = max(element.PISA_param(i,22),0);
element.PISA_param(i,22) = min(element.PISA_param(i,22),0.999);

% element.element.PISA_param = element.PISA_param; % save whole matrix into element structure
end
