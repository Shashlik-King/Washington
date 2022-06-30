function [output] = plot_functions_axial(element,pile,output,Coord,settings,plots,loads,Ex,Ey,Es,i,node,soil,L,G,ii,data,path)
%--------------------------------------------------------------------------
% PURPOSE
% Compute plots for axial calculation.
% 
% INPUT:  from main file (run_COSPIN.m)
%
% CODE            : EBGX
% APPROVED        : 
%
% MODIFIED        : EBGX    24.11.2016   Programming
%                 : FKMV    06.04.2020   Adjusting for AO1
%--------------------------------------------------------------------------
%% Initial
%--------------------------------------------------------------------------




output.plot_node  = node.level(1:end-1)<=soil.toplevel(1,1); % Determines which nodes are below seabed (embedded)
output.plot_level = node.level(output.plot_node);   % Sorts out the nodes below seabed
if (~isdeployed)
    data.save_path2 = [pwd,'/AppendixGenerationFiles/ProjectLocation'];
    data.save_path = [pwd,'/output'];
    pathname2 = [pwd,'/library/Output/Temporary_plots']; % temporary folder for plot to be loaded into database by perl, can be specified arbitrarily

else
    data.save_path2 = [path.output,'/AppendixGenerationFiles/ProjectLocation'];
    data.save_path = [path.output,'/output'];
    pathname2 = [path.output,'/library/Output/Temporary_plots']; % temporary folder for plot to be loaded into database by perl, can be specified arbitrarily

end
filename = [data.location,'_',soil.type,'_']; % general part of each filename

if ~exist(data.save_path2,'dir')
    mkdir (data.save_path2);
end
if ~exist(data.save_path,'dir')
    mkdir (data.save_path);
end
if ~exist(pathname2,'dir')
    mkdir (pathname2);
end
if ~exist([pathname2,'/Axial_capacity'],'dir')
    mkdir ([pathname2,'/Axial_capacity']);
end
if ~exist([pathname2,'/Axial_UR'],'dir')
    mkdir ([pathname2,'/Axial_UR']);
end
if ~exist([pathname2,'/Normal_force'],'dir')
    mkdir ([pathname2,'/Normal_force']);
end
if ~exist([pathname2,'/Axial_load_displacement'],'dir')
    mkdir ([pathname2,'/Axial_load_displacement']);
end
if ~exist([pathname2,'/Axial_displacement'],'dir')
    mkdir ([pathname2,'/Axial_displacement']);
end
%--------------------------------------------------------------------------
%% Plots
% -------------------------------------------------------------------------
%% Vertical displacement in last load step as a function of depth 
if plots.res_vs_pilelength == 1 && i == length(pile.length) && soil.psf==0
    clf
    figure(1)
    subplot(1,2,1)
    hold on
    plot(output.Rd_cum(:)/1000,-element.level(1:end-1,1))%,'-o',output.Rd(:,2),pile.length,'-o')
    plot([loads.Vc/1000 loads.Vc/1000],[-element.level(1,1) -element.level(end-1,1)*2], '--k')
    plot([0 output.Rd_cum(pile.axial_elem+1)/1000],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
    depth_util = interp1(loads.Vc./output.Rd_cum, -element.level(1:end-1,1),0.6);
    index = find(-element.level(1:end-1,1) > depth_util);
    plot([0 output.Rd_cum(index(1))/1000],[depth_util depth_util], '--b')
    plot(output.Rd_cum(index)/1000,-element.level(index,1),'b','Linewidth',4)
%     plot([0 output.Rd_cum(pile.axial_elem)/1000],[pile.L pile.L], '--r')
    hold off
    legend('Axial capacity', 'Factored axial load', 'Selected emb. length','Utilisation ratio - 0.6')
    set(gca,'YDir','rev')
    xlabel('Axial Compression Capacity [MN]')
    ylabel('Depth below mudline [m]')
    xlim([0 100])
    ylim([0 70])
    xticks([0 25 50 75 100])
	% grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    subplot(1,2,2)
    output.Rd_cum(find(output.Rd_cum<0))=0;
    plot(loads.Vc./output.Rd_cum, -element.level(1:end-1,1))
    hold on
    plot([0 loads.Vc./output.Rd_cum(pile.axial_elem+1)],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
    plot([0 0.6],[depth_util depth_util], '--b')
    plot([0.6 0.6],[depth_util depth_util*2], '--b')
    plot(loads.Vc./output.Rd_cum(index), -element.level(index,1),'b','Linewidth',4)
%     plot([0 loads.Vc./output.Rd_cum(pile.axial_elem)],[pile.L pile.L], '--r')
    set(gca,'YDir','rev')
    xlabel('Axial utilisation ratio')
    ylabel('Depth below mudline [m]')
    xlim([0.2 1])
    ylim([0 70])
    xticks([0.2 0.4 0.6 0.8 1])
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    legend('Axial utilisation ratio','Selected emb. length','Utilisation ratio - 0.6')
    if settings.save_plots
        saveas(gcf,[data.save_path,'/axial_capacity_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/axial_capacity.png'])												   
    end
    save_name2{1} = [pathname2,'/Axial_capacity/',filename,'axial_capacity.png']; % temporary file for plot to be loaded into database by perl
    print(figure(1),save_name2{1}, '-r300','-dpng'); % temporary file for plot to be loaded into database by perl
    %close(1)
    
    % save for table
    output.Rd_total = output.Rd_cum(pile.axial_elem)/1000;
    output.Rd_skin = output.Rbd_cum(pile.axial_elem)/1000;
    output.Rd_end_bearing = output.Rsd_cum(pile.axial_elem)/1000;
end

if plots.res_vs_pilelength == 1 && i == length(pile.length) && soil.psf==1
        clf
    figure(1)
    subplot(1,2,1)
    hold on
    plot(output.Rd_cum(:)/1000,-element.level(1:end-1,1))%,'-o',output.Rd(:,2),pile.length,'-o')
    plot([loads.Vc/1000 loads.Vc/1000],[-element.level(1,1) -element.level(end-1,1)*2], '--k')
    plot([0 output.Rd_cum(pile.axial_elem)/1000],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
    depth_util = interp1(loads.Vc./output.Rd_cum, -element.level(1:end-1,1),0.6);
    index = find(-element.level(1:end-1,1) > depth_util);
    plot([0 output.Rd_cum(index(1))/1000],[depth_util depth_util], '--b')
    plot(output.Rd_cum(index)/1000,-element.level(index,1),'b','Linewidth',4)
%     plot([0 output.Rd_cum(pile.axial_elem)/1000],[pile.L pile.L], '--r')
    hold off
    legend('Axial capacity', 'Factored axial load', 'Selected emb. length','Utilisation ratio - 0.6')
    set(gca,'YDir','rev')
    xlabel('Axial Compression Capacity [MN]')
    ylabel('Depth below mudline [m]')
    xlim([0 100])
    ylim([0 70])
    xticks([0 25 50 75 100])
	% grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    subplot(1,2,2)
    output.Rd_cum(find(output.Rd_cum<0))=0;
    plot(loads.Vc./output.Rd_cum, -element.level(1:end-1,1))
    hold on
    plot([0 loads.Vc./output.Rd_cum(pile.axial_elem)],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
    plot([0 0.6],[depth_util depth_util], '--b')
    plot([0.6 0.6],[depth_util depth_util*2], '--b')
    plot(loads.Vc./output.Rd_cum(index), -element.level(index,1),'b','Linewidth',4)
%     plot([0 loads.Vc./output.Rd_cum(pile.axial_elem)],[pile.L pile.L], '--r')
    set(gca,'YDir','rev')
    xlabel('Axial utilisation ratio')
    ylabel('Depth below mudline [m]')
    xlim([0.2 1])
    ylim([0 70])
    xticks([0.2 0.4 0.6 0.8 1])
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    legend('Axial utilisation ratio','Selected emb. length','Utilisation ratio - 0.6')
    if settings.save_plots
        saveas(gcf,[data.save_path,'/axial_capacity_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/axial_capacity.png'])												   
    end
    save_name2{1} = [pathname2,'/Axial_capacity/',filename,'axial_capacity.png']; % temporary file for plot to be loaded into database by perl
    print(figure(1),save_name2{1}, '-r300','-dpng'); % temporary file for plot to be loaded into database by perl
    %close(1)
end

%% Utilization ratio 
if plots.UR_axial == 1
    figure(2)
    output.UR_axial_plot = output.UR_axial(output.plot_node(1:end-1)); 
    plot(abs(output.UR_axial_plot),output.plot_level(1:end-1,1)-pile.head)   
%     eval(['title(''Axial utilization plot, L = ' num2str(L) ' m'');'])
    xlabel('Utilization ratio [-]')
    ylabel('Embedded length [m]')
    ylim([-pile.length(i) 0])
    %xlim([0.40 0.46])
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    title(['Axial UR - ',settings.analysis_loading{ii,1},' L = ' num2str(pile.length(i)) ' m'])
    if settings.save_plots
       saveas(gcf,[data.save_path,'/Axial_UR_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/Axial_UR.png'])											
    end
    save_name2{2} = [pathname2,'/Axial_UR/',filename,'Axial_UR.png'];
	print(figure(2),'-dpng',save_name2{2}, '-r300')
end

%% Normal force plot 
if plots.normal_force_axial == 1
    N = Es{end}(:,1,1);
%     for x=1:element.nelem
%     eldia2(Ex(x,:),Ey(x,:),N(x,:),[1 4],1.0);
%     end 
%     title(['Normal force plot, L = ' num2str(pile.length(i)) ' m'])
%     xlabel('Normal force [kN]')
%     ylabel('Depth [m]')
    figure(3)
    plot([Ey(:,1);Ey(end,2)],[N(:);0],'b')
    view([90 -90]) 
    ylabel('Normal force [kN]')
    xlabel('Depth [m]')
    xlim([-pile.length(i) 0])
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    title(['Normal force - ',settings.analysis_loading{ii,1},' L = ' num2str(pile.length(i)) ' m'])
    if settings.save_plots
       saveas(gcf,[data.save_path,'/Normal_force_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/Normal_force.png'])												
    end
    save_name2{3} = [pathname2,'/Normal_force/',filename,'Normal_force.png'];
	print(figure(3),'-dpng',save_name2{3}, '-r300')
end 
    
%% Resistance vs. pile length 
% if plots.res_vs_pilelength == 1 && i == length(pile.length)
%     figure
%     plot(output.Rd(:,1),pile.length,'-o',output.Rd(:,2),pile.length,'-o')
%     legend('Tension','Compression')
%     set(gca,'YDir','rev')
%     xlabel('Design bearing capacity [kN]')
%     ylabel('Pile length [m]')
% %     title('Axial bearing capacity, API method')
%     title(['Axial bearing capacity API - ',settings.analysis_loading{ii,1},' L = ' num2str(pile.length(i)) ' m'])
%     if settings.save_plots
%        saveas(gcf,[data.save_path,'/Axial_resistance_vs_pile_length_',data.location,'.png'])
%     end
% %     save_name2{2} = [pathname2,'/Axial_resistance_vs_pile_length/',filename,'Axial_resistance_vs_pile_length.png'];
% % 	print(figure(2),'-dpng',save_name2{2}, '-r300')
% end
   

%% Load-displacement curve 
if plots.load_deflection_axial == 1
    
    APILE.load = [0.2276E+01 0.2276E+02 0.1138E+03 0.2276E+03 0.1138E+04 0.2276E+04 0.8020E+04 0.1175E+05 0.1183E+05]; 
    APILE.disp = [0.2820E-05 0.2820E-04 0.1410E-03 0.2820E-03 0.1410E-02 0.2820E-02 0.1363E-01 0.2668E-01 0.5186E-01]; 

    if strcmp(settings.analysis_loading(ii),'Comp') == 1
        F = loads.Vc+sum(G.pile); 
        F_axis = loads.Vc;
    elseif strcmp(settings.analysis_loading(ii),'Tens') == 1
        F = -loads.Vt+sum(G.pile); 
        F_axis = -loads.Vt;
    end 
    if settings.n_max < 300
        disp('To get proper results for the load-deflection curve, settings.n_max should be at least 300')
    end
    deflections = output.deflections;
    loads_help = output.loads; 
    for iii = 1:settings.n_max % putting a column of zero deflection in the first column of the deflection matrix. This is done to ensure that the displacement plot passes through origo.
        output.deflections(:,settings.n_max+2-iii) = deflections(:,settings.n_max+1-iii);
        output.loads(:,settings.n_max+2-iii) = loads_help(:,settings.n_max+1-iii);
    end
    output.deflections(:,1) = zeros(size(output.deflections,1),1); % putting in zero displacement for zero load
    output.loads(:,1) = zeros(size(output.loads,1),1); % Putting in zero load for zero displacement
    axial_set = interp1((output.loads(2,1:min([settings.n_max output.n_possible])+1)),output.deflections(2,1:min([settings.n_max output.n_possible])+1)*1000,F_axis);
    limit = find((output.loads(2,1:min([settings.n_max output.n_possible])+1)) < F_axis);
    output.axial_set = axial_set;
    
    
    figure(4)
    hold on
%     plot(output.deflections(2,1:min([settings.n_max output.n_possible])+1)*1000,(output.loads(2,1:min([settings.n_max output.n_possible])+1))/(F),'-bx') % 
    plot(axial_set,loads.Vc,'or')    
%     plot(output.deflections(2,1:min([settings.n_max output.n_possible])+1)*1000,(output.loads(2,1:min([settings.n_max output.n_possible])+1)),'-kx') % 
    plot([output.deflections(2,1:limit(end))*1000 axial_set],[output.loads(2,1:limit(end)) F_axis],'-kx') % 
%     plot(APILE.disp*1000,APILE.load,'-bx') % 
    plot([axial_set axial_set],[0 loads.Vc],':r')
    plot([0 axial_set],[loads.Vc loads.Vc],':r')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    ylabel('Vertical force applied [kN]')
    xlabel('Vertical pile head deflection [mm]')
%     ylim([0 round(F_axis, -4)])
    ylim([0 F_axis])
    legend('Axial load applied','Location','SouthEast')
    title(['Pile head load-displacement - ',settings.analysis_loading{ii,1},' L = ' num2str(pile.length(i)) ' m'])
    if strcmp(settings.analysis_loading(ii),'Tens') == 1
        set(gca,'XDir','Reverse')
    end 
    
    if settings.save_plots
       saveas(gcf,[data.save_path,'/Axial_load_displacement_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/Axial_load_displacement.png'])														   
    end
    save_name2{4} = [pathname2,'/Axial_load_displacement/',filename,'Axial_load_displacement.png'];
	print(figure(4),'-dpng',save_name2{4}, '-r300')
end


%% Displacement curve 

if plots.deflection_plot_axial == 1
for j = 1:element.nelem
    output.ver_defl(j,i) = -output.deflections(j*3-1,min(settings.n_max,output.n_possible)); % save vertical displacement for plotting at final load step for each pile length 
end
    figure(5)
    subplot(1,2,ii)
    %plot(roundp(output.ver_defl(:,i),4)*10^3,Coord(:,2))
    plot(output.ver_defl(:,i)*10^3,Coord(:,2));
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('Vertical displacement [mm]')
    ylabel('Depth [m]')
    title(['Settlement - ',settings.analysis_loading{ii,1},' L = ' num2str(pile.length(i)) ' m'])
    ylim([-pile.length(i) 0]) 
    if settings.save_plots
       saveas(gcf,[data.save_path,'/Axial_displacement_',data.location,'.png'])
		saveas(gcf,[data.save_path2,'/Axial_displacement.png'])													  
    end
    save_name2{5} = [pathname2,'/Axial_displacement/',filename,'Axial_displacement.png'];
	print(figure(5),'-dpng',save_name2{5}, '-r300')
end 

end 