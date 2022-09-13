function [plug_unplug G output pile] = axial_uls(data,A,G,fs,qp,element,pile,SF,reduction,id,output,i,plots,loads) 
%% CALCULATION OF STATIC AXIAL CAPACITY
%--------------------------------------------------------------------------
% CHANGE LOG
% 16.07.2013    MUOE - PROGRAMMING
% 08.09.2013    MUOE - SEPARATING CALCULATION OF SKIN FRICTION AND TIP 
%                      RESISTANCE FROM THIS CALCULATION OF AXIAL ULS
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% Rd:   total design axial capacity of pile             [kN]
%--------------------------------------------------------------------------
%% Pre-allocating
%--------------------------------------------------------------------------
G.pile  = zeros(length(element.model_axial),1);
R_so    = zeros(length(element.model_axial),2);
R_si    = zeros(length(element.model_axial),2);
%--------------------------------------------------------------------------
%% Supplementary calculations
%--------------------------------------------------------------------------
for j = 1:(element.nelem-1)
    G.pile(j) = element.ep(j,2)*element.height(j)*(pile.density-10); % effective weight of the pile
    for top_bottom = 1:2
        R_so(j,top_bottom) = fs.o(j,top_bottom)*A.so(j,top_bottom); % outer skin friction force on the element
        R_si(j,top_bottom) = fs.i(j,top_bottom)*A.si(j,top_bottom); % inner skin friction force on the element
    end
end
R_so_total = sum(sum(R_so));
R_si_total = sum(sum(R_si));
%--------------------------------------------------------------------------
%% Determination of plugged/unplugged pile
%--------------------------------------------------------------------------
% Determine plugged/unplugged
    % Used to 1) determine correct axial capacity
    %         2) determine correct t-z curves for inside of pile
    
%% Plug/unplug control
% Tension
if pile.plug_unplug.tens == 0 % Automatic determination
    if R_si_total > sum(G.soil)
        plug_unplug.tens_index = 'Plugged';
    else
        plug_unplug.tens_index = 'Unplugged';
    end
elseif pile.plug_unplug.tens == 1 % Force plugged
    plug_unplug.tens_index = 'Plugged';
elseif pile.plug_unplug.tens == 2 % Force unplugged
    plug_unplug.tens_index = 'Unplugged';
else
    disp('Choose plug_unplug_tens = [0,1,2]')
end

% Compression
if pile.plug_unplug.comp == 0 % Automatic determination
    if R_si_total+qp(end)*pile.cross_section.endarea+sum(G.soil) > qp(end)*pi/4*pile.diameter^2
        plug_unplug.comp_index = 'Plugged';
    else
        plug_unplug.comp_index = 'Unplugged';
    end
elseif pile.plug_unplug.comp == 1 % Force plugged
    plug_unplug.comp_index = 'Plugged';
elseif pile.plug_unplug.comp == 2 % Force unplugged
    plug_unplug.comp_index = 'Unplugged';
else
    disp('Choose plug_unplug_comp = [0,1,2]')
end
        
%% Bearing capacity
% Tension
if strcmp(plug_unplug.tens_index,'Plugged') % Plugged
    Rsk.tens = reduction.skin_tension*R_so_total; % characteristic skin resistance in tension
    Rbk.tens = 0; % characteristic tip resistance in tension
    Rbd.tens = Rbk.tens/SF.tens; % design tip resistance in tension
    Rsd.tens = Rsk.tens/SF.tens; % design skin resistance in tension
    Rd.tens  = Rsd.tens + Rbd.tens + sum(G.pile)*SF.weight_advan + sum(G.soil)*SF.weight_advan; % summing up the total (design) axial forces acting on the pile in tension
else % Unplugged
    Rsk.tens = reduction.skin_tension*R_so_total+reduction.skin_tension*R_si_total; % characteristic skin resistance in tension
    Rbk.tens = 0; % characteristic tip resistance in tension
    Rbd.tens = Rbk.tens/SF.tens; % design tip resistance in tension
    Rsd.tens = Rsk.tens/SF.tens; % design skin resistance in tension
    Rd.tens  = Rsd.tens + Rbd.tens + sum(G.pile)*SF.weight_advan; % summing up the total (design) axial forces acting on the pile in tension
end
% Compression
if strcmp(plug_unplug.comp_index,'Plugged') % Plugged
    Rbk.comp = qp(end)*pi/4*pile.diameter^2;
    Rsk.comp = R_so_total;
    Rbd.comp = Rbk.comp/SF.comp;
    Rsd.comp = Rsk.comp/SF.comp;
    Rd.comp  = Rsd.comp + Rbd.comp - sum(G.pile) - sum(G.soil); 
    R_so(end,:) = [];
    G.pile(end,:) = [];
%     Rd.cumulative=R_so(:,1)/SF.comp + (qp(:,1)*pi/4*pile.diameter^2)/SF.comp - (G.pile) - (G.soil'); 
else % Unplugged
    Rbk.comp = qp(end)*pile.cross_section.endarea;
    Rsk.comp = R_so_total+R_si_total; 
    Rbd.comp = Rbk.comp/SF.comp;
    Rsd.comp = Rsk.comp/SF.comp;
    Rd.comp  = Rsd.comp + Rbd.comp - sum(G.pile);  
     
end

     
for j = 1:(element.nelem-1)
    
    Rsd.cumulative_comp(j)=(sum(R_so(1:j,1)+R_so(1:j,2))+sum(R_si(1:j,1)+R_si(1:j,2)))/SF.comp;
    Rbd.cumulative_comp(j)=(qp(j)*pile.cross_section.endarea)/SF.comp;
    Rd.cumulative_comp(j)= Rsd.cumulative_comp(j)+Rbd.cumulative_comp(j)- sum(G.pile(1:j));
    
end

% Save for plotting
output.Rd(i,1:2)  = [Rd.tens Rd.comp];
output.Rbd(i,1:2) = [Rbd.tens Rbd.comp];
output.Rsd(i,1:2) = [Rsd.tens Rsd.comp];
output.Rd_cum  = Rd.cumulative_comp; 
output.Rbd_cum = Rbd.cumulative_comp; 
output.Rsd_cum = Rsd.cumulative_comp;

if plots.res_vs_pilelength         == 1
    Min=(abs(abs(element.level(:,1))-(pile.L-pile.extra_L)));
%     Min=(abs(abs(element.level(:,1))-(pile.L)));
    pile.axial_elem=find(Min==min(Min))-1;
    if size(pile.axial_elem,1) > 1
        pile.axial_elem = pile.axial_elem(1);
    end
end
   
if plots.res_vs_pilelength         == 1
    output.axial_util_ratio=(loads.Vc/Rd.cumulative_comp((pile.axial_elem)));
    axial_capacity=Rd.cumulative_comp((pile.axial_elem));
else
	pile.axial_elem = size(element.level(:,1),1);
	output.axial_util_ratio = (loads.Vc/Rd.cumulative_comp((pile.axial_elem-1)));
    axial_capacity=Rd.comp;
end