function [element] = A_factor(element,pile,loads,i)

%%%This is for permanent rotations with N>100%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                    


for top_bottom=1:2 % A corresponding to considered pile segment [-]
    if strcmp(loads.static_cyclic,'cyclic') % cyclic case
        if strcmp(loads.A,'API')
            element.A(i,top_bottom)    =   0.9; % acc. API   
        elseif strcmp(loads.A,'TUHH')
            element.A(i,top_bottom)    =   min(0.343*element.heqv(i,top_bottom)/pile.diameter,0.9);  % for more than 100 load cycles, acc. EA-Piles (p.443, Equ. D3.6)
        end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(loads.static_cyclic,'static') % static case
        element.A(i,top_bottom) = max(0.9 , 3.0 - 0.8*element.heqv(i,top_bottom)/pile.diameter); %acc. API
    end
end