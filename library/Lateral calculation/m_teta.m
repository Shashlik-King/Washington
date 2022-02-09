function [m teta output toe_plot_teta] = m_teta(node,pile,element,loads,output,settings)
%% FORMULATION FOR p-y CURVES
%--------------------------------------------------------------------------
% CHANGE LOG
% 01.08.2013    MUOE - SEPARATING FROM WINKLER_PY.M
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% pile.xx:      pile.diameter is used to determine yc
% element.xx:   model.py is used to determine which p-y curve to use
%               epsilon50 is used to determine yc
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% p:        lateral resistance
% y:        lateral displacement
%--------------------------------------------------------------------------
%% Initial
%--------------------------------------------------------------------------
nelem = length(element.model_py);
y_tot = zeros(nelem,2);
%--------------------------------------------------------------------------
%% Calculation routine
%--------------------------------------------------------------------------
for c = 1:(nelem-1)
    teta_top(c,:) = [0 0.000005 0.000025 0.00005 0.00015 0.00025 0.0005 0.000625 0.00075 0.001 0.0015 0.002	0.003 0.004 0.006 0.01 0.02];
    teta_bot(c,:) = teta_top(c,:);
    upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4 0.8];
    upy_bot(c,:) = upy_top(c,:);    

    
    npoints = size(teta_top,2);
	npoints2 = size(upy_top,2);
    
    for j = 1:npoints2
		[ksppynode_top ksppynode_bot] = secspringstiff(element,pile,loads,[upy_top(c,j) upy_bot(c,j)],c);
        ppynode_top(c,j) = ksppynode_top*upy_top(c,j);  
        ppynode_bot(c,j) = ksppynode_bot*upy_bot(c,j); 
		for k=1:npoints
			[ksmtetanode_top ksmtetanode_bot] = secmomstiff(element,pile,[teta_top(c,k) teta_bot(c,k)],[upy_top(c,j) upy_bot(c,j)],ksppynode_top,ksppynode_bot,c);
			mmtetanode_top{1,j}(c,k) = ksmtetanode_top*teta_top(c,k);  
			mmtetanode_bot{1,j}(c,k) = ksmtetanode_bot*teta_bot(c,k);  
		end
    end
    
     
    
         % Base moment shear
         for g = nelem % the last element, toe shear
             tetabot = [0 0.0000003 0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.01 0.1];
           for  j = 1:length(tetabot)
                 %
                 if settings.toe_shear == 1
    %                 output.ts = toe_shear(element,node,pile,ybot,c);
               [kMb_top(j) kMb_bot(j)] = secMBstiff(element,pile,[tetabot(j) tetabot(j)],g); 
                Mb_bot(j) = kMb_bot(j)*tetabot(j); 
                output.Mb = Mb_bot*(element.level(end,1)-element.level(end,2));
                elseif settings.toe_shear == 0
                     output.Mb = 0;
                end
           end
			 toe_plot_teta = tetabot;
         end
end

    
    for i =2:element.nelem-1
        for j=1:length(ppynode_top(1,:))
            if ppynode_top(i,j)<ppynode_bot(i-1,j)||ppynode_top(i,j)>ppynode_bot(i-1,j)
                p_av_M(i,j) = (ppynode_top(i,j).*(abs(element.level(i,2))-...
                    abs(element.level(i,1)))+ppynode_bot(i-1,j).*...
                    (abs(element.level(i,1))-abs(element.level(i-1,1))))...
                    ./((abs(element.level(i,2))-abs(element.level(i,1)))...
                    +(abs(element.level(i,1))-abs(element.level(i-1,1))));
            else
                p_av_M(i,j) = ppynode_top(i,j);
            end
        end
    end
    












% Putting p and y data into same format as t and z
for j=1:npoints2
	for i = 1:size(upy_top,1)
		m.top{1,j}(i,:)      = mmtetanode_top{1,j}(i,:);
		m.bottom{1,j}(i,:)   = mmtetanode_bot{1,j}(i,:);
		teta.top{1,j}(i,:)      = teta_top(i,:);
		teta.bottom{1,j}(i,:)   = teta_bot(i,:);
    end
    
    for i=2:size(upy_top,1)  % start from second element , beacuse frist element is always zero
    %% There would be a table of the cells 
    %% first row is the table of M-T
    %%%from the second row until the end there is P_Y of each node (from 2
    %%%until the last node)
    m.top{1+i,j}(1,:)=uniquetol(p_av_M(i,:),0.0001);
	m.bottom{1+i,j}(1,:)=uniquetol(p_av_M(i,:),0.0001);
    SizeofUnique=size(m.bottom{1+i,j}(1,:),2);
    teta.top{1+i,j} (1,:)=upy_top(i,1:SizeofUnique);
	teta.bottom{1+i,j}(1,:) =upy_bot(i,1:SizeofUnique);
    
    if strcmp(element.model_py(i),'Zero soil')
        m.top{1+i,j}(1,:)=([0:1:size(ppynode_top(i,:),2)-1]/(size(ppynode_top(i,:),2)-1))*1000;
        m.bottom{1+i,j}(1,:)=([0:1:size(ppynode_bot(i,:),2)-1]/(size(ppynode_bot(i,:),2)-1))*1000;
    end 
    end 
    m.top{2,j}=uniquetol(p_av_M(2,:),0.0001); 
    m.bottom{2,j}=uniquetol(p_av_M(2,:),0.0001);
    
end
end 