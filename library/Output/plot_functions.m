function [output] = plot_functions(element,pile,node,soil,plots,output,settings,i, loads, data,SF)
%% PLOTTING FUNCTION
% RECEIVES INPUT FROM MAIN.M TO MAKE PLOTS
%--------------------------------------------------------------------------
% CHANGE LOG
% 02.08.2013    MUOE - OPTIMIZING PLOT, SINCE CERTAIN DATA IS REDUNDANT
%--------------------------------------------------------------------------
%% Initial
%--------------------------------------------------------------------------
output.plot_node  = node.level<=soil.toplevel(1,1); % Determines which nodes are below seabed (embedded)
output.plot_level = node.level(output.plot_node);   % Sorts out the nodes below seabed
plots.node_rot_def= 1;                              % number of node to plot permanent rotations for, 1 = node at pile head

for k = 1:6
    save_name2{k}     = '0';                     % initialise filenames of plots, temporary plot to be inserted into database
end

data.save_path2 = [pwd,'\AppendixGenerationFiles\ProjectLocation'];
data.save_path = [pwd,'\output'];
pathname2 = [pwd,'\library\Output\Temporary_plots']; % temporary folder for plot to be loaded into database by perl, can be specified arbitrarily
filename = [data.location,'_',soil.type,'_']; % general part of each filename

%% write documentation of revision settings for each position
% by use of the revision no. the applied loads, geometry and soil
% parameters can be backtracked within the MySQL database
% F = fopen([pathname,data.location,'_AppliedRevisions_',date,'.dat'],'wt');
% fprintf(F,'**********************************************************************************************************************************************\n');
% fprintf(F,'**************************************** Revisions applied for all checks at position %5s **************************************************\n',data.location);
% fprintf(F,'************************************************************* %s ******************************************************************\n',date);
%     fprintf(F,'   \n');
%     fprintf(F,'Global revision: %2d\n',data.revision.global);
%     fprintf(F,'   \n');
%     fprintf(F,'Soil revision: %2d\n',data.revision.soil);
%     fprintf(F,'   \n');
%     fprintf(F,'Structural revision: %2d\n',data.revision.structure);
%     fprintf(F,'   \n');
%     fprintf(F,'Load revision: %2d\n',data.revision.loads);
%     fprintf(F,'   \n');
% fprintf(F,'**********************************************************************************************************************************************\n');
% fclose(F);

%--------------------------------------------------------------------------
%%  Resistance vs. pile depth plots
%--------------------------------------------------------------------------

% if plots.res_vs_pilelength == 1 && i == length(pile.length) && soil.psf==0
%     clf
%     figure(1)
%     subplot(1,2,1)
%     hold on
%     plot(Rd.cumulative_comp(:)/1000,-element.level(1:end-1,1))%,'-o',output.Rd(:,2),pile.length,'-o')
%     plot([loads.Vc/1000 loads.Vc/1000],[-element.level(1,1) -element.level(end-1,1)], '--k')
%     plot([0 Rd.cumulative_comp(pile.axial_elem)/1000],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
%     hold off
%     legend('Axial capacity', 'Factored axial load', 'Selected emb. length')
%     set(gca,'YDir','rev')
%     xlabel('Axial Compression Capacity [MN]')
%     ylabel('Depth below mudline [m]')
%     xlim([0 100])
%     ylim([0 70])
%     xticks([0 25 50 75 100])
%     subplot(1,2,2)
%     Rd.cumulative_comp(find(Rd.cumulative_comp<0))=0;
%     plot(loads.Vc./Rd.cumulative_comp, -element.level(1:end-1,1))
%     hold on
%     plot([0 loads.Vc./Rd.cumulative_comp(pile.axial_elem)],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
%     set(gca,'YDir','rev')
%     xlabel('Axial utilisation ratio')
%     ylabel('Depth below mudline [m]')
%     xlim([0.2 1])
%     ylim([0 70])
%     xticks([0.2 0.4 0.6 0.8 1])
%     legend('Axial utilisation ratio','Selected emb. length')
%     if settings.save_plots
%         saveas(gcf,[data.save_path,'\axial_capacity_',data.location,'.png'])
%     end
%     save_name2{1} = [pathname2,'\Axial_capacity\',filename,'axial_capacity.png']; % temporary file for plot to be loaded into database by perl
%     print(figure(1),save_name2{1}, '-r300','-dpng'); % temporary file for plot to be loaded into database by perl
%     %close(1)
% end

if plots.res_vs_pilelength == 1 && i == length(pile.length) && soil.psf==1
    clf
    figure(1)
    subplot(1,2,1)
    hold on
    plot(output.Rd_cum(:)/1000,-element.level(1:end-1,1))%,'-o',output.Rd(:,2),pile.length,'-o')
    plot([loads.Vc/1000 loads.Vc/1000],[-element.level(1,1) -element.level(end-1,1)], '--k')
    plot([0 output.Rd_cum(pile.axial_elem)/1000],[pile.L-pile.extra_L pile.L-pile.extra_L], '--r')
    hold off
    legend('Axial capacity', 'Factored axial load', 'Selected emb. length')
    set(gca,'YDir','rev')
    xlabel('Axial comp. capacity [MN]')
    ylabel('Depth below mudline [m]')
    xlim([0 100])
    ylim([0 ceil(-element.level(end-1,1))])
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
    set(gca,'YDir','rev')
    xlabel('Axial utilisation ratio')
    ylabel('Depth below mudline [m]')
    xlim([0 1])
    ylim([0 ceil(-element.level(end-1,1))])
    xticks([0.2 0.4 0.6 0.8 1])
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    legend('Axial utilisation ratio','Selected emb. length')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\axial_capacity_',data.location,'.png'])
        saveas(gcf,[data.save_path2,'\axial_capacity.png'])
    end
    save_name2{1} = [pathname2,'\Axial_capacity\',filename,'axial_capacity.png']; % temporary file for plot to be loaded into database by perl
    print(figure(1),save_name2{1}, '-r300','-dpng'); % temporary file for plot to be loaded into database by perl
    %close(1)
end

%--------------------------------------------------------------------------
%%  Critical pile length plots
%--------------------------------------------------------------------------

if plots.pilehead_vs_length == 1 && i == length(pile.length) && soil.psf==0
    
    asym = min(output.pilehead_rotation(1,:));
    change = (output.pilehead_rotation-asym)/asym;
    
    if pile.fixed_lenght_switch
        %%% Interpolation to find relative rotation that corresponds to critical L
        D_L=abs(pile.length_start_or-pile.length);
        L_pos=find(D_L==min(D_L));
        if length(L_pos)==2
            L_pos=min(L_pos);
        end
        
        if pile.length(L_pos)>=pile.length_start_or
            Lx1=L_pos-1;
            Lx2=L_pos;
        else
            Lx1=L_pos;
            Lx2=L_pos+1;
        end
        Ly1=output.pilehead_rotation(find(pile.length==pile.length(Lx1)));
        Ly2=output.pilehead_rotation(find(pile.length==pile.length(Lx2)));
        Lyx=interp1([pile.length(Lx1) pile.length(Lx2)],[Ly1 Ly2],pile.length_start_or);
        
        Ry1=change(Lx1);
        Ry2=change(Lx2);
        
        Ryx=interp1([pile.length(Lx1) pile.length(Lx2)],[Ry1 Ry2],pile.length_start_or);
    end
    %%%
    figure(2)
    clf;
    [Critpilelength, output] = DeterCritPileLength(change,output,pile,node); %function to get interpolated critical pile length
    subplot(2,1,1)
    hold on
    plot(pile.length,output.pilehead_rotation)
    if pile.fixed_lenght_switch
        plot(pile.length_start_or,Lyx, 's','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
        text(pile.length_start_or+0.5,Lyx+0.01,'Embedded length')
    end
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('Pile length [m]')
    ylabel('Pile head rotation [{\circ}]')
    hold off
    subplot(2,1,2)
    hold on
    plot(pile.length,change*100)
    if pile.fixed_lenght_switch
        plot(pile.length_start_or,Ryx*100, 's','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
        text(pile.length_start_or+0.5,Ryx*100+4,'Embedded length')
    else
        plot(Critpilelength.value,pile.criterion, 's','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
        text(Critpilelength.text.x,Critpilelength.text.y,Critpilelength.text.value)
    end
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('Pile length [m]')
    ylabel('Relative pile head rotation [%]')
    hold off
    if settings.save_plots
        saveas(gcf,[data.save_path,'\pilehead_vs_length_',data.location,'.png'])
        saveas(gcf,[data.save_path2,'\pilehead_vs_length.png'])
    end
    save_name2{2} = [pathname2,'\pilehead_vs_length\',filename,'pilehead_vs_length.png'];
    print(figure(2),'-dpng',save_name2{2}, '-r300')
    %close(2)
end

% if plots.pilehead_vs_length == 1 && i == length(pile.length) && soil.psf==1
%
%     asym = min(output.pilehead_rotation(1,:));
%     change = (output.pilehead_rotation-asym)/asym;
%
%     if pile.fixed_lenght_switch
%          %%% Interpolation to find relative rotation that corresponds to critical L
%         D_L=abs(pile.length_start_or-pile.length);
%         L_pos=find(D_L==min(D_L));
%         if length(L_pos)==2
%             L_pos=min(L_pos);
%         end
%
%         if pile.length(L_pos)>=pile.length_start_or
%             Lx1=L_pos-1;
%             Lx2=L_pos;
%         else
%             Lx1=L_pos;
%             Lx2=L_pos+1;
%         end
%         Ly1=output.pilehead_rotation(find(pile.length==pile.length(Lx1)));
%         Ly2=output.pilehead_rotation(find(pile.length==pile.length(Lx2)));
%         Lyx=interp1([pile.length(Lx1) pile.length(Lx2)],[Ly1 Ly2],pile.length_start_or);
%
%         Ry1=change(Lx1);
%         Ry2=change(Lx2);
%
%         Ryx=interp1([pile.length(Lx1) pile.length(Lx2)],[Ry1 Ry2],pile.length_start_or);
%     end
%   %%%
%     figure(2)
%     clf;
%
% 	[Critpilelength, output] = DeterCritPileLength(change,output,pile,node); %function to get interpolated critical pile length
%     subplot(2,1,1)
%     hold on
%     plot(pile.length,output.pilehead_rotation)
%     if pile.fixed_lenght_switch
%         plot(pile.length_start_or,Lyx, 's','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
%         text(pile.length_start_or+0.5,Lyx+0.01,'Embedded length')
%     end
%     % grid setup
%     grid
%     ax.GridLineStyle = '-';
%     ax.GridAlpha = 0.5;
%     grid minor
%     ax.MinorGridLineStyle = ':';
%     ax.MinorGridAlpha = 0.2;
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     xlabel('Pile length [m]')
%     ylabel('Pile head rotation [{\circ}]')
%     hold off
%     subplot(2,1,2)
%     hold on
%     plot(pile.length,change*100)
%     if pile.fixed_lenght_switch
%         plot(pile.length_start_or,Ryx*100, 's','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
%         text(pile.length_start_or+0.5,Ryx*100+4,'Embedded length')
%     end
%     % grid setup
%     grid
%     ax.GridLineStyle = '-';
%     ax.GridAlpha = 0.5;
%     grid minor
%     ax.MinorGridLineStyle = ':';
%     ax.MinorGridAlpha = 0.2;
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     xlabel('Pile length [m]')
%     ylabel('Relative pile head rotation [%]')
%     hold off
%     if settings.save_plots
%         saveas(gcf,[data.save_path,'\pilehead_vs_length_',data.location,'.png'])
% 		saveas(gcf,[data.save_path2,'\pilehead_vs_length.png'])
%     end
%     save_name2{2} = [pathname2,'\pilehead_vs_length\',filename,'pilehead_vs_length.png'];
% 	print(figure(2),'-dpng',save_name2{2}, '-r300')
%     %close(2)
% end
%--------------------------------------------------------------------------
%%  Deflection plots
%--------------------------------------------------------------------------

if plots.deflection_plot == 1 && soil.psf == 0
    Coord = output.Coord;
    figure(3)
    plot(output.hor_defl(1:end-1,end),Coord(1:end-1,2))
    xlabel('Horizontal deflection [m]')
    ylabel('Level [m vref]')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\deflection_plot_',data.location,'_analysisID=',data.analysis_id,'.png'])
        saveas(gcf,[data.save_path2,'\deflection_plot.png'])
    end
    save_name2{3} = [pathname2,'\Deflection\',filename,'deflection.png'];
    print(figure(3),save_name2{3},'-dpng', '-r300');
    %close(3)
end

if plots.deflection_plot == 1 && soil.psf == 1
    Coord = output.Coord;
    figure(3)
    plot(output.hor_defl(1:end-1,end),Coord(1:end-1,2))
    xlabel('Horizontal deflection [m]')
    ylabel('Level [m vref]')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\deflection_plot_',data.location,'_analysisID=',data.analysis_id,'.png'])
        saveas(gcf,[data.save_path2,'\deflection_plot.png'])
    end
    save_name2{3} = [pathname2,'\Deflection\',filename,'deflection.png'];
    print(figure(3),save_name2{3}, '-r300','-dpng');
    %close(3)
end

%--------------------------------------------------------------------------
%%  UR along pile plot
%--------------------------------------------------------------------------

if plots.utilization_ratio == 1 && soil.psf==0
    figure(4)
    output.UR_plot = output.UR(output.plot_node(1:end-2),2);
    plot(output.UR_plot,output.plot_level(2:end-1,1)-pile.head)
    eval(['title(''Utilization plot, L = ' num2str(pile.L) ' '');'])
    xlabel('Utilization ratio [-]')
    ylabel('Embedded length [m]')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\utilization_ratio_',data.location,'.png'])
    end
    save_name2{4} = [pathname2,'\Utilization_ratio\',filename,'utilization_ratio.png'];
    print(figure(4),'-dpng',save_name2{4}, '-r300')
end

if plots.utilization_ratio == 1 && soil.psf==1
    figure(4)
    output.UR_plot = output.UR(output.plot_node(1:end-2),2);
    plot(output.UR_plot,output.plot_level(2:end-1,1)-pile.head)
    eval(['title(''Utilization plot, L = ' num2str(pile.L) ' '');'])
    xlabel('Utilization ratio [-]')
    ylabel('Embedded length [m]')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\utilization_ratio_',data.location,'.png'])
        saveas(gcf,[data.save_path2,'\utilization_ratio.png'])
    end
    save_name2{4} = [pathname2,'\Utilization_ratio\',filename,'utilization_ratio.png'];
    print(figure(4),'-dpng',save_name2{4}, '-r300')
end

if plots.utilization_ratio == 1
    figure(44)
    output.p_UR_plot = output.p_UR(output.plot_node(1:end-2),2);
    output.pu_UR_plot = output.pu_UR(output.plot_node(1:end-2),2);
    hold on
    plot(output.p_UR_plot,output.plot_level(2:end-1,1)-pile.head)
    plot(output.pu_UR_plot,output.plot_level(2:end-1,1)-pile.head)
    eval(['title(''Soil utilization, L = ' num2str(pile.L) ' '');'])
    xlabel('Soil resistance [kN/m]')
    ylabel('Embedded length [m]')
    legend('p','p_u_l_t')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    if settings.save_plots
        saveas(gcf,[data.save_path,'\utilization_resistance_',data.location,'.png'])
        saveas(gcf,[data.save_path2,'\utilization_resistance.png'])
    end
    save_name2{44} = [pathname2,'\Utilization_resistance\',filename,'utilization_resistance.png'];
    print(figure(44),'-dpng',save_name2{4}, '-r300')
end



%--------------------------------------------------------------------------
%%  Deflection bundle plots
%--------------------------------------------------------------------------

if plots.deflection_bundle == 1
    figure(100)
    hold all
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    plot(output.hor_defl(:,i),output.plot_level(:,1)-pile.head,'k')
    title('Deflection plot, bundle')
    xlabel('Deflection [m]')
    ylabel('Depth below pile head [m]')
end

if settings.toe_shear == 1
    if plots.toe_shear_graph == 1
        toe_plot_u = output.toe_plot_u;
        figure
        % ts = toe_shear(element,node,pile,abs(output.hor_defl(end)),length(element.model_py));
        [~, kHb_bot] = secHBstiff(element,pile,[abs(output.hor_defl(end)) abs(output.hor_defl(end))],length(element.model_py));
        ts = kHb_bot*abs(output.hor_defl(end))*(element.level(end,1)-element.level(end,2));
        hold all
        % grid setup
        grid
        ax.GridLineStyle = '-';
        ax.GridAlpha = 0.5;
        grid minor
        ax.MinorGridLineStyle = ':';
        ax.MinorGridAlpha = 0.2;
        set(gca,'XMinorTick','on','YMinorTick','on')
        plot(toe_plot_u*1000,output.ts,'b-')
        plot(abs(output.hor_defl(end))*1000,ts,'ro')
        ylabel('Toe shear force [kN]')
        xlabel('Toe kick [mm]')
    end
end

%--------------------------------------------------------------------------
%%  Permanent deformtion plots
%--------------------------------------------------------------------------

if plots.permanent_rot_def == 1 && soil.psf==0
    if settings.n_max < 10
        disp('To get proper results for the permanent rotation, settings.n_max should be at least 10')
    end
    deflections = output.deflections;
    for ii = 1:settings.n_max % putting a column of zero deflection in the first column of the deflection matrix. This is done to ensure that the rotation plot passes through origo.
        output.deflections(:,settings.n_max+2-ii) = deflections(:,settings.n_max+1-ii);
    end
    output.deflections(:,1) = zeros(size(output.deflections,1),1); % putting in zero displacement for zero load
    F = linspace(0,1,settings.n_max+1); % this is valid because the load is applied in equally sized steps - the magnitude of the load doesn't matter, only the fact that it is applied in equally sized steps
    output.perm_rot = output.deflections(3*plots.node_rot_def,end)-F(end)/output.elasticstiff; % the permanent rotation
    figure(5)
    clf;
    hold on
    scatter(output.perm_rot*180/pi,F(1),'ro')
    plot(output.deflections(3*plots.node_rot_def,:)*180/pi,F,'-bx')
    plot([output.perm_rot output.deflections(3*plots.node_rot_def,end)]*180/pi,[F(1) F(end)],'-r')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    eval(['title(''Permanent rotation plot, L = ' num2str(pile.L) ' '');'])
    eval(['xlabel(''Rotation of node ' num2str(plots.node_rot_def) ' [^o] '');'])
    ylabel('Load ratio [-]')
    eval(['legend(''Permanent rotation is ' num2str(output.perm_rot*180/pi,'%.3f') '^o '');'])
    if settings.save_plots
        saveas(gcf,[data.save_path,'\permanent_rot_def_',data.location,'_',num2str(loads.n_cycles) ,'_analysisID=', data.analysis_id, '.png'])
        saveas(gcf,[data.save_path2,'\permanent_rot_def_',data.location,'_',num2str(loads.n_cycles),'_analysisID=', data.analysis_id, '.png'])
    end
    % % %     save_name2{5} = [pathname2,'\Perm_rot\',filename,'perm_rotation.png'];
    % % %     print(figure(5),'-dpng',save_name2{5}, '-r300');
    %close(5)
end

if plots.permanent_rot_def == 1 && soil.psf==1
    if settings.n_max < 10
        disp('To get proper results for the permanent rotation, settings.n_max should be at least 10')
    end
    deflections = output.deflections;
    for ii = 1:settings.n_max % putting a column of zero deflection in the first column of the deflection matrix. This is done to ensure that the rotation plot passes through origo.
        output.deflections(:,settings.n_max+2-ii) = deflections(:,settings.n_max+1-ii);
    end
    output.deflections(:,1) = zeros(size(output.deflections,1),1); % putting in zero displacement for zero load
    F = linspace(0,1,settings.n_max+1); % this is valid because the load is applied in equally sized steps - the magnitude of the load doesn't matter, only the fact that it is applied in equally sized steps
    output.perm_rot = output.deflections(3*plots.node_rot_def,end)-F(end)/output.elasticstiff; % the permanent rotation
    figure(5)
    clf;
    hold on
    scatter(output.perm_rot*180/pi,F(1),'ro')
    plot(output.deflections(3*plots.node_rot_def,:)*180/pi,F,'-bx')
    plot([output.perm_rot output.deflections(3*plots.node_rot_def,end)]*180/pi,[F(1) F(end)],'-r')
    % grid setup
    grid
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.5;
    grid minor
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    set(gca,'XMinorTick','on','YMinorTick','on')
    eval(['title(''Permanent rotation plot, L = ' num2str(pile.L) ' '');'])
    eval(['xlabel(''Rotation of node ' num2str(plots.node_rot_def) ' [^o] '');'])
    ylabel('Load ratio [-]')
    eval(['legend(''Permanent rotation is ' num2str(output.perm_rot*180/pi,'%.3f') '^o '');'])
    if settings.save_plots
        saveas(gcf,[data.save_path,'\permanent_rot_def_',data.location,'_analysisID=', data.analysis_id,'.png'])
        saveas(gcf,[data.save_path2,'\permanent_rot_def.png'])
    end
    save_name2{5} = [pathname2,'\Perm_rot\',filename,'perm_rotation.png'];
    print(figure(5),'-dpng',save_name2{5}, '-r300');
    %close(5)
end
%--------------------------------------------------------------------------
%%  Load deflection plots
%--------------------------------------------------------------------------

if plots.load_deflection == 1 % && strcmp(settings.interface,'FAC')
    if settings.n_max < 50
        disp('To get proper results for the load-deflection and UR curve, settings.n_max should be at least 50; this is AO1 project specific')
    end
    deflections = output.deflections;
    for ii = 1:settings.n_max % putting a column of zero deflection in the first column of the deflection matrix. This is done to ensure that the displacement plot passes through origo.
        output.deflections(:,settings.n_max+2-ii) = deflections(:,settings.n_max+1-ii);
    end
    
    output.deflections(:,1) = zeros(size(output.deflections,1),1); % putting in zero displacement for zero load
    
    %  this is valid because the load is applied in equally sized steps - the magnitude of the load doesn't matter, only the fact that it is applied in equally sized steps
    
    F = linspace(0,settings.max_load_ratio,settings.n_max+1)*loads.H; % DNV calc, material factors are applied
    
    if ~isfield(output,'n_possible')
        error('Increase max_load_ratio in order to plot the load-deflection curve!')
    end
    
    output.defl_plot=output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000;
    output.load_interp=F(1:min([settings.n_max output.n_possible])+1);
    
    %%%% finding horizontal load that that corresponds to displacement=10%D
    value=find(abs(output.defl_plot-pile.diameter*100)==(min(abs(output.defl_plot-pile.diameter*100))));
    
    if length(value)==2
        value=value(1);
    end
    
    if output.defl_plot(value)<=pile.diameter*100
        value1=value;
        value2=value+1;
    elseif output.defl_plot(value)>pile.diameter*100
        value1=value-1;
        value2=value;
    end
    
    valuex=interp1([output.defl_plot(value1) output.defl_plot(value2)],[output.load_interp(value1) output.load_interp(value2)],pile.diameter*100);
    FF_interp = F/loads.H; % factor multiplied to horizontal load for load-deflection plot
    H_load_10_ratio=interp1([output.load_interp(value1) output.load_interp(value2)],[FF_interp(value1) FF_interp(value2)],valuex);
    H_load_10 = H_load_10_ratio*loads.H; % load that corresponds to displacement=10%D
    % matrix_name=['H_load_10_pos_',data.location,'.mat'];
    % save(['output\rev0.1\mat_files/',matrix_name], 'H_load_10');
    
    % interpolation to find the displacement that corresponds to loads.H UR
    load_10=find(abs(output.load_interp-loads.H)==(min(abs(output.load_interp-loads.H))));
    if length(load_10)==2
        load_10=load_10(1);
    end
    
    if output.load_interp(load_10)<=loads.H
        load1=load_10;
        load2=load1+1;
    elseif output.load_interp(load_10)>loads.H
        load1=load_10-1;
        load2=load_10;
    end
    
    
    %%%% DNV-GL plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if  plots.load_deflection_type == 0
        
        disp_10=interp1([output.load_interp(load1) output.load_interp(load2)],[output.defl_plot(load1) output.defl_plot(load2)],loads.H);
        
        figure(4)
        clf;
        output.H_load_10_dnv=H_load_10;
        output.UR_lat_dnv=loads.H/output.H_load_10_dnv;
        
        x1=[0 pile.diameter*100];
        y1=[output.H_load_10_dnv/1000 output.H_load_10_dnv/1000];
        x2=[pile.diameter*100 pile.diameter*100];
        y2=[0 output.H_load_10_dnv/1000];
        
        plot(pile.diameter*100,output.H_load_10_dnv/1000,'or')
        hold on
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %disp_10
        
        plot(disp_10,loads.H/1000,'og')
        plot(output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000,F(1:min([settings.n_max output.n_possible])+1)/1000,'-kx')
        plot(x1,y1,':r')
        plot(x2,y2,':r')
        % grid setup
        grid
        ax.GridLineStyle = '-';
        ax.GridAlpha = 0.5;
        grid minor
        ax.MinorGridLineStyle = ':';
        ax.MinorGridAlpha = 0.2;
        set(gca,'XMinorTick','on','YMinorTick','on')
        eval(['title(''Pile head load-displacement curve, UR_D_N_V_G_L = ' num2str(output.UR_lat_dnv,'%0.2f') ' '');'])
        %     xlim([0 1000])
        legend('Lateral Capacity','ALS Unfactored Load','Location','SouthEast')
        %legend('Lateral Capacity','ULS Factored Load','Location','SouthEast')
        xlim([0 pile.diameter*120])
        ylabel('Horizontal load at mudline [MN]')
        xlabel('Pile head deflection [mm]')
        output.def_calibration = output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000;
        output.force_calibration = F(1:min([settings.n_max output.n_possible])+1);
        if settings.save_plots
            saveas(gcf,[data.save_path,'\utilization_ratio_',data.location,'_DNV', '_analysisID=', data.analysis_id,'.png'])
            saveas(gcf,[data.save_path2,'\utilization_ratio_DNV.png'])
        end
        % % % %         save_name2{4} = [pathname2,'\Utilization_ratio\',filename,'utilization_ratio_DNV.png'];
        % % % %         print(figure(4),'-dpng',save_name2{4}, '-r300');
        
        %%%% BSH plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % % % % % % % % % % % % % % % % % % % % % %     elseif soil.psf==0 && plots.load_deflection_type == 0
        % % % % % % % % % % % % % % % % % % % % % %         disp_10=interp1([output.load_interp(load1) output.load_interp(load2)],[output.defl_plot(load1) output.defl_plot(load2)],loads.H);
        % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % %         figure(4)
        % % % % % % % % % % % % % % % % % % % % % %         clf;
        % % % % % % % % % % % % % % % % % % % % % %         output.H_load_10_bsh=H_load_10;
        % % % % % % % % % % % % % % % % % % % % % %         output.UR_lat_bsh=loads.H/output.H_load_10_bsh;
        % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % %         x1=[0 pile.diameter*100];
        % % % % % % % % % % % % % % % % % % % % % %         y1=[output.H_load_10_bsh/1000 output.H_load_10_bsh/1000];
        % % % % % % % % % % % % % % % % % % % % % %         x2=[pile.diameter*100 pile.diameter*100];
        % % % % % % % % % % % % % % % % % % % % % %         y2=[0 output.H_load_10_bsh/1000];
        % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % %         plot(pile.diameter*100,output.H_load_10_bsh/1000,'or')
        % % % % % % % % % % % % % % % % % % % % % %         hold on
        % % % % % % % % % % % % % % % % % % % % % %         plot(disp_10,loads.H/1000,'og')
        % % % % % % % % % % % % % % % % % % % % % %         plot(output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000,F(1:min([settings.n_max output.n_possible])+1)/1000,'-kx')
        % % % % % % % % % % % % % % % % % % % % % %         plot(x1,y1,':r')
        % % % % % % % % % % % % % % % % % % % % % %         plot(x2,y2,':r')
        % % % % % % % % % % % % % % % % % % % % % %         % grid setup
        % % % % % % % % % % % % % % % % % % % % % %         grid
        % % % % % % % % % % % % % % % % % % % % % %         ax.GridLineStyle = '-';
        % % % % % % % % % % % % % % % % % % % % % %         ax.GridAlpha = 0.5;
        % % % % % % % % % % % % % % % % % % % % % %         grid minor
        % % % % % % % % % % % % % % % % % % % % % %         ax.MinorGridLineStyle = ':';
        % % % % % % % % % % % % % % % % % % % % % %         ax.MinorGridAlpha = 0.2;
        % % % % % % % % % % % % % % % % % % % % % %         set(gca,'XMinorTick','on','YMinorTick','on')
        % % % % % % % % % % % % % % % % % % % % % %         eval(['title(''Pile head load-displacement curve, UR_B_S_H = ' num2str(output.UR_lat_bsh,'%0.2f') ' '');'])
        % % % % % % % % % % % % % % % % % % % % % %         xlim([0 1000])
        % % % % % % % % % % % % % % % % % % % % % %         legend('Lateral Capacity','ULS Factored Load','Location','SouthEast')
        % % % % % % % % % % % % % % % % % % % % % %         ylabel('Horizontal load at mudline [MN]')
        % % % % % % % % % % % % % % % % % % % % % %         xlabel('Pile head deflection [mm]')
        % % % % % % % % % % % % % % % % % % % % % %         output.def_calibration = output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000;
        % % % % % % % % % % % % % % % % % % % % % %         output.force_calibration = F(1:min([settings.n_max output.n_possible])+1);
        % % % % % % % % % % % % % % % % % % % % % %         if settings.save_plots
        % % % % % % % % % % % % % % % % % % % % % %             saveas(gcf,[data.save_path,'\utilization_ratio_',data.location,'_BSH.png'])
        % % % % % % % % % % % % % % % % % % % % % %             saveas(gcf,[data.save_path2,'\utilization_ratio_BSH.png'])
        % % % % % % % % % % % % % % % % % % % % % %         end
        % % % % % % % % % % % % % % % % % % % % % %         save_name2{4} = [pathname2,'\Utilization_ratio\',filename,'utilization_ratio_BSH.png'];
        % % % % % % % % % % % % % % % % % % % % % %         print(figure(4),'-dpng',save_name2{4}, '-r300');
    end
    
    %%%% Load-disp plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plots.load_deflection == 1
        if plots.load_deflection_type == 1
            
            if settings.n_max < 50
                disp('To get proper results for the load-deflection curve, settings.n_max should be at least 50')
            end
            output.deflections(:,1) = zeros(size(output.deflections,1),1); % putting in zero displacement for zero load
            F = linspace(0,settings.max_load_ratio,settings.n_max+1)*loads.H; % this is valid because the load is applied in equally sized steps - the magnitude of the load doesn't matter, only the fact that it is applied in equally sized steps
            figure(6)
            hold on
            plot(output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000,F(1:min([settings.n_max output.n_possible])+1)/1000,'-bx')
            grid on
            eval(['title(''Pile head load-displacement curve plot, L = ' num2str(pile.L) ' ' data.location ' '');'])
            xlim([0 1000])
            legend('1D Model','3D Model')
            ylabel('Pile Head Force [MN]')
            xlabel('Pile head deflection [mm]')
            
            output.def_calibration = output.deflections(1,1:min([settings.n_max output.n_possible])+1)*1000;
            output.force_calibration = F(1:min([settings.n_max output.n_possible])+1);
            if settings.save_plots
                saveas(gcf,[data.save_path,'\load_deflection_',data.location,'.png'])
                saveas(gcf,[data.save_path2,'\load_deflection.png'])
            end
            save_name2{6} = [pathname2,'\load_deflection\',filename,'load_deflection.png'];
            print(figure(6),'-dpng',save_name2{6}, '-r300');
        end
    end
    
    %% Internal forces plots
    
    if plots.moment_distribution == 1
        figure('visible','on')
        subplot(1,3,1)
        plot(Es{1,1}(:,1,3),element.level(:,1)) %FKMV
        %     plot(soil.Mxy*1000, soil.levelPLAX, '-rx')
        ylabel('Depth [m]')
        xlabel('Moment [kNm]')
        grid on
        subplot(1,3,2)
        plot(Es{1,1}(:,1,2),element.level(:,1)) %FKMV
        %     plot(soil.Mxy*1000, soil.levelPLAX, '-rx')
        ylabel('Depth [m]')
        xlabel('Shear [kNm]')
        grid on
        subplot(1,3,3)
        plot(output.hor_defl(1:end-1,end),element.level(:,1)) %FKMV
        %     plot(soil.Mxy*1000, soil.levelPLAX, '-rx')
        ylabel('Depth [m]')
        xlabel('Displacement [kNm]')
        grid on
        if settings.save_plots
            saveas(gcf,[data.save_path,'\struct_forces',data.location,'.png'])
        end
        save_name2{6} = [pathname2,'\struct_forces\',filename,'struct_forces.png'];
        print(figure(6),'-dpng',save_name2{6}, '-r300');
    end
    
end
if settings.update_db, database_write(output,settings,loads,data,plots,save_name2,pile,i,soil); end


end