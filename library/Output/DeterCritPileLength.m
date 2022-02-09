function [Critpilelength, output] = DeterCritPileLength(change,output,pile,node) %SPSO 01-02-2019 having a function in the bottom of such a file here is different than what we typically have done in COSPIN. Typically we have one file for each subroutine. Both of course works, but maybe to keep it uniform throughout we should make a seperate file for this subroutine?
%created 16/01/2019 by GMME as update to COSPIN
    if change(1,1) < pile.criterion/100
        disp('--------------------------------------')
        disp('Critical pile length due to 10% criterium can''t be determined.')   %SPSO: 01-02-2019 Changed to write 10% and not 1%
        error('Investigated minimum pile length is too large.')  					%SPSO: 01-02-2019 Do we want the program to report an error and stop? We could consider to just 1) exit function + 2) report that critical pile length cannot be determined due to too large minimum pile length + 3) make typical critical pile length plot but without the marker for the critical pile length 
        
    else
        for kl = 1:pile.criterion    % loop over several percentage criterion (10% ... 1%)
            ij = 2;     % counter
            crit_percent(kl) = (pile.criterion - (kl-1))/100;   % determine the decimal fraction corresponding to percentage (0.1 0.09 ... 0.01) %SPSO 01-02-2019 Changed decimal fraction in the brackets of explanation
            while change(1,ij) > crit_percent(kl)   % look for pile length until %-criterium is reached
                ij = ij + 1;
            end
            % interpolate critical pile length 
            output.pile_length_deter.crit_length(kl) = pile.length(ij-1) + (pile.length(ij) - pile.length(ij-1))...
                / (change(1,ij) - change (1,ij-1)) * (crit_percent(kl) - change(1,ij-1));
            % interpolate pile head deflection 
            output.pile_length_deter.pile_head_defl(kl) = output.hor_defl(1,ij-1) + (output.hor_defl(1,ij) - output.hor_defl(1,ij-1))...
                / (change(1,ij) - change (1,ij-1)) * (crit_percent(kl) - change(1,ij-1));
            % interpolate pile toe deflection 
            m = 1; % counter
            while -node.level(m) < pile.length(ij-1)
                m = m+1;
            end
            n = 1; % counter
            while -node.level(n) < pile.length(ij)
                n = n+1;
            end
            output.pile_length_deter.pile_toe_defl(kl) = output.hor_defl(m,ij-1) + (output.hor_defl(n,ij) - output.hor_defl(m,ij-1))...
                / (change(1,ij) - change (1,ij-1)) * (crit_percent(kl) - change(1,ij-1));
        end
    end

    kl = pile.criterion-(pile.criterion-1); % 1%-> kl = 10; 10% -> kl = 1; 5% -> kl = 6
    Critpilelength.value = output.pile_length_deter.crit_length(kl);
    Critpilelength.text.x = output.pile_length_deter.crit_length(kl)+0.1;
    Critpilelength.text.y = (pile.criterion/100+0.005)*100;
    Critpilelength.text.value = [' (',num2str(round(crit_percent(kl)*100)),'% : L = ',num2str(ceil(output.pile_length_deter.crit_length(kl)*10)/10),' m)'];	
end