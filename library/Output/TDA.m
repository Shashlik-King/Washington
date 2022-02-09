function [p_TDA, t_TDA, level_TDA]   = TDA(node,p,t,data,scour,pile,y,z,settings)
%%calculate springs at 1 m interval

node.level_top(1)=0;
for i=2:length(node.level)-1
    node.level_top(i)=node.level(i);
    node.level_bot(i-1)=node.level_top(i);
end
node.level_bot(length(node.level)-1)=min(node.level);

if pile.diameter>9
	element_size=1.5; % distance between soil springs. Must not be less than 10% of pile diameter
elseif pile.diameter>4
	element_size=1.0; % distance between soil springs. Must not be less than 10% of pile diameter
elseif settings.nelem_factor==3  && settings.GeniE
	element_size=0.5; % distance between soil springs. Must not be less than 10% of pile diameter
elseif settings.nelem_factor>=4  && settings.GeniE
	element_size=0.25; % distance between soil springs. Must not be less than 10% of pile diameter
else
	element_size=1.0; % distance between soil springs. Must not be less than 10% of pile diameter
end

level_TDA=0:-element_size:min(node.level);
level_TDAtop(1)=0;
level_TDAtop(2:length(level_TDA))=(level_TDA(2:length(level_TDA))+level_TDA(1:length(level_TDA)-1))/2;
level_TDAbot(1:length(level_TDA)-1)=level_TDAtop(2:length(level_TDA));
level_TDAbot(length(level_TDA))=min(node.level);

p_TDA=zeros(length(level_TDA),size(p.top,2));
t_TDA=zeros(length(level_TDA),size(t.top,2));
for i=1:length(level_TDA)
    for j=1:length(node.level)-2
        if node.level_top(j)<=level_TDAtop(i) && node.level_top(j)>level_TDAbot(i)
            if node.level_bot(j)>level_TDAbot(i)
                p_TDA(i,:)=p_TDA(i,:)+(p.top(j,:)+p.bottom(j,:))/2*(node.level_top(j)-node.level_bot(j))/(level_TDAtop(i)-level_TDAbot(i));
                t_TDA(i,:)=t_TDA(i,:)+(t.top(j,:)+t.bottom(j,:))/2*(node.level_top(j)-node.level_bot(j))/(level_TDAtop(i)-level_TDAbot(i));
            else
                p_TDA(i,:)=p_TDA(i,:)+(p.top(j,:)*2+(p.bottom(j,:)-p.top(j,:))*(node.level_top(j)-level_TDAbot(i))/(node.level_top(j)-node.level_bot(j)))/2*(node.level_top(j)-level_TDAbot(i))/(level_TDAtop(i)-level_TDAbot(i));
                t_TDA(i,:)=t_TDA(i,:)+(t.top(j,:)*2+(t.bottom(j,:)-t.top(j,:))*(node.level_top(j)-level_TDAbot(i))/(node.level_top(j)-node.level_bot(j)))/2*(node.level_top(j)-level_TDAbot(i))/(level_TDAtop(i)-level_TDAbot(i));
            end
        elseif node.level_bot(j)<level_TDAtop(i) && node.level_bot(j)>=level_TDAbot(i)
                p_TDA(i,:)=p_TDA(i,:)+(p.top(j,:)+(p.bottom(j,:)-p.top(j,:))*(level_TDAtop(i)-node.level_bot(j))/(node.level_top(j)-node.level_bot(j))+p.bottom(j,:))/2*(level_TDAtop(i)-node.level_bot(j))/(level_TDAtop(i)-level_TDAbot(i));
                t_TDA(i,:)=t_TDA(i,:)+(t.top(j,:)+(t.bottom(j,:)-t.top(j,:))*(level_TDAtop(i)-node.level_bot(j))/(node.level_top(j)-node.level_bot(j))+t.bottom(j,:))/2*(level_TDAtop(i)-node.level_bot(j))/(level_TDAtop(i)-level_TDAbot(i));
        end
    end
end

% avg_TDA2=mean(p_TDA(:,2))
% avg_ptop2=mean(p.top(:,2))
% avg_pbot2=mean(p.bottom(:,2))
% avgp2=(avg_ptop2+avg_pbot2)/2
% 
% figure(101)
% plot(p.top(:,2),node.level_top(1:end-1),'b.',p.bottom(:,2),node.level_bot(1:end-1),'r.',p_TDA(:,2),level_TDA,'k-')
% ylabel('Level [m Vref]')
% xlabel('p [kN/m]')
% legend('p.top','p.bottom','p_T_D_A')
% 
% % avg_TDA4=mean(p_TDA(:,4))
% % avg_ptop4=mean(p.top(:,4))
% % avg_pbot4=mean(p.bottom(:,4))
% % avgp4=(avg_ptop4+avg_pbot4)/2
% 
% figure(102)
% plot(p.top(:,4),node.level_top(1:end-1),'b.',p.bottom(:,4),node.level_bot(1:end-1),'r.',p_TDA(:,4),level_TDA,'k-')
% ylabel('Level [m Vref]')
% xlabel('p [kN/m]')
% legend('p.top','p.bottom','p_T_D_A')
% 
% figure(103)
% plot(t.top(:,2),node.level_top(1:end-1),'b.',t.bottom(:,2),node.level_bot(1:end-1),'r.',t_TDA(:,2),level_TDA,'k-')
% ylabel('Level [m Vref]')
% xlabel('p [kN/m]')
% legend('t.top','t.bottom','t_T_D_A','location','Southeast')
% 
% figure(104)
% plot(t.top(:,4),node.level_top(1:end-1),'b.',t.bottom(:,4),node.level_bot(1:end-1),'r.',t_TDA(:,4),level_TDA,'k-')
% ylabel('Level [m Vref]')
% xlabel('t [kN/m^2]')
% legend('t.top','t.bottom','t_T_D_A','location','Southeast')

% %generate txt file to TDA (similar to IBDAS format
% name = sprintf('ANSYSinput.%s-%s_%s%s_%s-GEOrev%s',data.project,data.location,data.type,data.curves,settings.lateralmultipliers,data.soil_rev);
% F = fopen(name,'wt');
% 
% gamma_tot=1.0;
% 
% %%% Header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(F,'%s %s %s %s','* SOIL PROFILE FOR LOCATION:', data.location ,'at', data.project ,'');
% fprintf(F,'\n');
% fprintf(F,'********************************************************************************\n');
% fprintf(F,'%s %0.2f', '* LOCAL SCOUR [m]:',scour.local);
% fprintf(F,'\n');
% fprintf(F,'%s %0.2f', '* WATER DEPTH [m]:',scour.water_depth);
% fprintf(F,'\n');
% fprintf(F,'%s %.2f', '* PILE LENGTH [m]:',pile.length(end));
% fprintf(F,'\n');
% fprintf(F,'%s %.2f', '* PILE DIAMETER [m]:',pile.diameter(end));
% fprintf(F,'\n');
% fprintf(F,'%s %.2f', '* PILE HEAD [m VREF]:',pile.head);
% fprintf(F,'\n');
% fprintf(F,'%s %.2f', '* PILE STICK-UP [m]:',pile.stick_up);
% fprintf(F,'\n');
% fprintf(F,'%s %s', '* FILE CREATED BY:',data.prepared_by,',',datestr(now));
% fprintf(F,'\n');
% fprintf(F,'%s','* p[kN/m], y[m], t[kN/m] and z[m]');
% fprintf(F,'\n');
% fprintf(F,'%s','SOIL');
% fprintf(F,'\n');
% fprintf(F,'%s %f','  MUDD',-scour.water_depth);
% fprintf(F,'\n');
% 
% %%%% p-y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(F,'%s','                              Displacement [m]');
% fprintf(F,'\n');
% fprintf(F,'%s %s','Depth','Level');
% for ii=1:size(y.top,2)
%     fprintf(F,'%5.4f   ',y.top(1,ii));
% end
% fprintf(F,'\n');
% fprintf(F,'%s  %s  %s','[m]','[m]','                          Lateral resistance [kN/m]');
% fprintf(F,'\n');
% for jj=1:size(p_TDA,1)
%     fprintf(F,'%8.1f  %8.1f  ',-level_TDA(1,jj),-scour.water_depth+level_TDA(1,jj));
%     for kk=1:size(y.top,2)
%         fprintf(F,'%8.1f  ',p_TDA(jj,kk));
%     end
%     fprintf(F,'\n');
% end
% fprintf(F,'\n');
% 
% %%%% t-z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(F,'%s','                              Settlement [m]');
% fprintf(F,'\n');
% fprintf(F,'%s %s','Depth','Level');
% out_circ=pi*pile.diameter;
% for ii=1:size(z.top,2)
%     fprintf(F,'%5.3f   ',z.top(1,ii));
% end
% fprintf(F,'\n');
% fprintf(F,'%s  %s  %s','[m]','[m]','                          Axial resistance [kN/m]');
% fprintf(F,'\n');
% for jj=1:size(t_TDA,1)
%     fprintf(F,'%5.1f  %5.1f  ',-level_TDA(1,jj),-scour.water_depth+level_TDA(1,jj));
%     for kk=1:size(z.top,2)
%         fprintf(F,'%5.1f  ',t_TDA(jj,kk)*out_circ);
%     end
%     fprintf(F,'\n');
% end
% 
% %%%% Q-z
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT READY
% fprintf(F,'%s','                              Settlement [m]');
% fprintf(F,'\n');
% fprintf(F,'%s %s','Depth','Level');
% out_circ=pi*pile.diameter;
% for ii=1:size(z.top,2)
%     fprintf(F,'%5.3f   ',zQ.bottom(1,ii)/5);
% end
% fprintf(F,'\n');
% fprintf(F,'%s  %s  %s','[m]','[m]','                          Axial resistance [kN/m]');
% fprintf(F,'\n');
% for jj=1:size(t_TDA,1)
%     fprintf(F,'%5.1f  %5.1f  ',-level_TDA(1,jj),-scour.water_depth+level_TDA(1,jj));
%     for kk=1:size(z.top,2)
%         fprintf(F,'%5.1f  ',t_TDA(jj,kk)*out_circ);
%     end
%     fprintf(F,'\n');
% end
% 
