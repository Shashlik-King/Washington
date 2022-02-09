function ANSYS_ASAS_TDA(p,t,level,y,z,scour,pile,data,iii,Q,zQ,plug_unplug,element)
%% FILE FOR GENERATING ANSYS-INPUT
%--------------------------------------------------------------------------
% CHANGE LOG
% 2014.01.22        MMOL        Programming
%--------------------------------------------------------------------------

name = ['ANSYSinp.',data.location,'.txt'];
F = fopen(name,'wt');

gamma_tot=1.0;
iii;
%%%% Header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(F,'%s %s %s %s','* SOIL PROFILE FOR LOCATION:', data.location ,'at', data.project ,'');
fprintf(F,'\n');
fprintf(F,'********************************************************************************\n');
fprintf(F,'%s %0.2f', '* LOCAL SCOUR [m]:',scour.local);
fprintf(F,'\n');
fprintf(F,'%s %0.2f', '* WATER DEPTH [m]:',scour.water_depth);
fprintf(F,'\n');
fprintf(F,'%s %.2f', '* PILE LENGTH [m]:',pile.length(iii));
fprintf(F,'\n');
fprintf(F,'%s %.2f', '* PILE DIAMETER [m]:',pile.diameter(iii));
fprintf(F,'\n');
fprintf(F,'%s %.2f', '* PILE HEAD [m VREF]:',pile.head);
fprintf(F,'\n');
fprintf(F,'%s %.2f', '* PILE STICK-UP [m]:',pile.stick_up);
fprintf(F,'\n');
fprintf(F,'%s %s', '* FILE CREATED BY:',data.prepared_by,',',datestr(now));
fprintf(F,'\n');
fprintf(F,'%s','* p[kN/m], y[m], t[kN/m] and z[m]');
fprintf(F,'\n');
fprintf(F,'%s','SOIL');
fprintf(F,'\n');
fprintf(F,'%s %f','  MUDD',-scour.water_depth);
fprintf(F,'\n');
if size(y,2)>17
	disp('p-y curves for ANSYS file given for 17 points. Check if modifications are needed')
end
%%%% p-y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(p,1)-1
    if i<size(p,1)
        fprintf(F,'%s %f %f', '1 P-Y', -level(i), -level(i+1));
    else
        fprintf(F,'%s %f %f', '1 P-Y', -level(i), pile.length(iii));
    end
    fprintf(F,'\n');
    fprintf(F,'%s %4.1f ','P',p(i,1));
    fprintf(F,'%4.1f ',p(i,2));
    fprintf(F,'%4.1f ',p(i,3));
    fprintf(F,'%4.1f ',p(i,4));
    fprintf(F,'%4.1f ',p(i,5));
    fprintf(F,'%4.1f ',p(i,6));
    fprintf(F,'%4.1f ',p(i,7));
	fprintf(F,'\n');
    fprintf(F,'%s %4.1f ',':',p(i,8));
	fprintf(F,'%4.1f ',p(i,9));
	fprintf(F,'%4.1f ',p(i,10));
	fprintf(F,'%4.1f ',p(i,11));
	fprintf(F,'%4.1f ',p(i,12));
	fprintf(F,'%4.1f ',p(i,13));
	fprintf(F,'%4.1f ',p(i,14));
	fprintf(F,'\n');
    fprintf(F,'%s %4.1f ',':',p(i,15));
	fprintf(F,'%4.1f ',p(i,16));
	fprintf(F,'%4.1f ',p(i,17));
	fprintf(F,'%4.1f ',p(i,17));% same p value applied at last two y-values
    fprintf(F,'\n');
    fprintf(F,'%s %4.4f ','Y',y.top(1,1));
    fprintf(F,'   %4.4f',y.top(1,2));
    fprintf(F,'   %4.4f',y.top(1,3));
    fprintf(F,'   %4.4f ',y.top(1,4));
    fprintf(F,'   %4.4f ',y.top(1,5));
    fprintf(F,'   %4.4f ',y.top(1,6));
    fprintf(F,'   %4.4f ',y.top(1,7));
	fprintf(F,'\n');
    fprintf(F,'%s %4.4f ',':',y.top(1,8));
    fprintf(F,'   %4.4f ',y.top(1,9));
	fprintf(F,'   %4.4f ',y.top(1,10));
	fprintf(F,'   %4.4f ',y.top(1,11));
	fprintf(F,'   %4.4f ',y.top(1,12));
	fprintf(F,'   %4.4f ',y.top(1,13));
	fprintf(F,'   %4.4f ',y.top(1,14));
	fprintf(F,'\n');
    fprintf(F,'%s %4.4f ',':',y.top(1,15));
    fprintf(F,'   %4.4f ',y.top(1,16));
	fprintf(F,'   %4.4f ',y.top(1,17));
	fprintf(F,'   %4.4f ',y.top(1,17)*2);
    fprintf(F,'\n');    
end
%%%% t-z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t[kN/m]
out_circ=pi*pile.diameter;
fprintf(F,'\n'); 
for i=1:size(t,1)-1
    if i<size(p,1)
        fprintf(F,'%s %f %f', '1 T-Z', -level(i), -level(i+1));
    else
        fprintf(F,'%s %f %f', '1 T-Z', -level(i), pile.length(iii));
    end
    fprintf(F,'\n');
    
    fprintf(F,'%s %4.1f ','T',t(i,1)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,2)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,3)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,4)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,5)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,6)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,7)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,8)/gamma_tot*out_circ);
	fprintf(F,'\n');
	fprintf(F,'%s %4.1f ',':',t(i,9)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,10)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,11)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,12)/gamma_tot*out_circ);
	fprintf(F,'%4.1f ',t(i,13)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,14)/gamma_tot*out_circ);
    fprintf(F,'%4.1f ',t(i,15)/gamma_tot*out_circ);
    fprintf(F,'\n');
    fprintf(F,'%s %4.3f ','Z',z.top(1,1)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,2)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,3)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,4)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,5)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,6)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,7)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,8)/gamma_tot);
	fprintf(F,'\n');
	fprintf(F,'%s %4.3f ',':',z.top(i,9)/gamma_tot);
	fprintf(F,'%4.3f ',z.top(1,10)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,11)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,12)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,13)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,14)/gamma_tot);
    fprintf(F,'%4.3f ',z.top(1,15)/gamma_tot);
    fprintf(F,'\n'); 
end

%%%% Q-z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(plug_unplug.comp_index,'Unplugged')
        Q_area = (pile.diameter^2-(pile.diameter-2*element.thickness(end))^2);
elseif strcmp(plug_unplug.comp_index,'Plugged')
        Q_area = pile.diameter^2;
end
fprintf(F,'\n');
fprintf(F,'%s','1 ENDB');
fprintf(F,'\n'); 
fprintf(F,'%s %4.1f ','t',Q.bottom(end,1)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,1)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,2)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,3)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,4)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,5)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,6)/gamma_tot*Q_area);
fprintf(F,'%4.1f ',Q.bottom(end,7)/gamma_tot*Q_area);
fprintf(F,'\n'); 
fprintf(F,'%s %4.3f ','z',-0.15);
fprintf(F,'%4.3f ',zQ.bottom(end,1)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,2)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,3)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,4)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,5)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,6)/gamma_tot);
fprintf(F,'%4.3f ',zQ.bottom(end,7)/gamma_tot);
fprintf(F,'\n'); 
fprintf(F,'%s','END');
fclose(F);