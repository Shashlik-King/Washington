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
    teta_top(c,:) = [0 0.0000003 0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.01 0.1];
    teta_bot(c,:) = teta_top(c,:);
    upy_top(c,:) = [0 0.0001 0.001 0.01 0.03 0.12 0.4];
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

% Putting p and y data into same format as t and z
for j=1:npoints2
	for i = 1:size(upy_top,1)
		m.top{1,j}(i,:)      = mmtetanode_top{1,j}(i,:);
		m.bottom{1,j}(i,:)   = mmtetanode_bot{1,j}(i,:);
		teta.top{1,j}(i,:)      = teta_top(i,:);
		teta.bottom{1,j}(i,:)   = teta_bot(i,:);
	end
	teta.top{2,j}=upy_top(end,j);
	teta.bottom{2,j}=upy_bot(end,j);
    m.top{2,j}=ppynode_top(end,j);
	m.bottom{2,j}=ppynode_bot(end,j);
end