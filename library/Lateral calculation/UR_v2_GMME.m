function output = UR_v2_GMME(element,settings,pile,loads,data,plots,ustep,y_tot,output,node)
%% UTILIZATION RATIO
% CALCULATES THE UTILIZATION RATIO FOR THE SOIL
%--------------------------------------------------------------------------
% CHANGE LOG
% 06.08.2013    MUOE - MAKING IT WORK WITH MAIN.M
% 30.04.2018 	GMME - Proper UR when using p- and y- mutlipliers
%--------------------------------------------------------------------------
%% Pre-allocation
%--------------------------------------------------------------------------
nelem           = element.nelem;
output.UR       = zeros(nelem,2); % initializing numerical integrator
output.pu_UR    = zeros(nelem,2); % initializing numerical integrator
output.p_UR     = zeros(nelem,2); % initializing numerical integrator
output.Bhc      = 0; % initializing numerical integrator
output.PU       = 0; % initializing numerical integrator
zero_defl       = 0; % maker for point of zero deflection
[SF reduction]  = factors(loads); % include partial safety factors
%--------------------------------------------------------------------------
%% Calculation routine
%--------------------------------------------------------------------------
for c = 1:nelem
    upy_top=abs(ustep(c*3-2,settings.n_max));
    upy_bot=abs(ustep((c+1)*3-2,settings.n_max)); 
    teta_top=abs(ustep(c*3,settings.n_max));
    teta_bot=abs(ustep(c*3+3,settings.n_max));
    %% API Sand
    if strcmp(element.model_py(c),'API sand') || strcmp(element.model_py(c),'Kirsch sand') || strcmp(element.model_py(c),'Kallehave sand')
        if element.pu(c,1) == 0
            element.pu(c,1) = 1;
        elseif element.pu(c,2) == 0
            element.pu(c,2) = 1;
        end       
        % -------- Calculation of UR --------------------------------------------        
        [output.UR(c,1),output.UR(c,2), output.p_UR(c,:), output.pu_UR(c,:)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,upy_top*1e3,upy_bot*1e3); % calculate UR
   %% API clay
    elseif (strcmp(element.model_py(c),'API clay') || strcmp(element.model_py(c),'Kirsch soft clay')) 
        % calculation of utilisation ratio:
        if strcmp(loads.static_cyclic,'static') % static case -> pu equal to pu_static
            [output.UR(c,1),output.UR(c,2), output.p_UR(c,1), output.p_UR(c,2), output.pu_UR(c,1),  output.pu_UR(c,2)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,upy_top*1e3,upy_bot*1e3);
        elseif strcmp(loads.static_cyclic,'cyclic') % cyclic case -> pu equal to 0.72*pu_static
            if upy_top<3*y_tot(c,1) % before peak (0.72*pu_static)
                [output.UR(c,1), ~, output.p_UR(c,1), ~, output.pu_UR(c,1), ~] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,3*y_tot(c,1),3*y_tot(c,2));           
            else % after peak (0.72*pu_static) -> UR = 1.0
                output.UR(c,1) = 1; % for UR-plot UR = 1 if capacity is exeeded
			end
			if upy_bot<3*y_tot(c,2) % before peak (0.72*pu_static)
                [~,output.UR(c,2), ~, output.p_UR(c,2), ~,  output.pu_UR(c,2)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,3*y_tot(c,1),3*y_tot(c,2));           
			else
				output.UR(c,2) = 1; % for UR-plot UR = 1 if capacity is exeeded
			end
        end
    %% stiff clay w/o free water, weak rock 
    elseif (strcmp(element.model_py(c),'Stiff clay w/o free water') || strcmp(element.model_py(c),'Kirsch stiff clay') || strcmp(element.model_py(c),'Modified Weak rock'))
        % calculation of p_mobilised (same method for static and cyclic case)
        [output.UR(c,1),output.UR(c,2), output.p_UR(c,1), output.p_UR(c,2), output.pu_UR(c,1),  output.pu_UR(c,2)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,upy_top*1e3,upy_bot*1e3); % calculate UR
    %% Zero soil    
    elseif strcmp(element.model_py(c),'Zero soil')
        output.UR(c,:)= [1,1];
        output.p_UR(c,:)=0;
        output.pu_UR(c,:)=0;
    %% PISA models 
    elseif (strcmp(element.model_py(c),'PISA sand') || strcmp(element.model_py(c),'PISA clay'))
	     [output.UR(c,1),output.UR(c,2), output.mUR(c,1),output.mUR(c,2), ... % UR values for PISA curves
		  output.p_UR(c,1), output.p_UR(c,2), output.m_UR(c,1), output.m_UR(c,2),... %  mobilised horizontal/rotational forces for PISA curves
		  output.pu_UR(c,1),  output.pu_UR(c,2), output.mu_UR(c,1),  output.mu_UR(c,2)] = Calc_UR_PISA(element,pile,loads,c,upy_top,upy_bot,teta_top,teta_bot,upy_top*1e3,upy_bot*1e3); % calculate UR
    %% Reese stiff clay
    elseif strcmp(element.model_py(c),'Reese stiff clay')
        % calculation of p_mobilised (same method for static and cyclic case)
		% calculation of A-factor
        xtop    = element.heqv(c,1);    % [m] equivalent height of element top node
        xbot    = element.heqv(c,2);    % [m] equivalent height of element bottom node
        D       = pile.diameter;        % [m] pile outer diameter
        % -------- Determination of A ---------------------------------------------
		Atop=min(0.0106*(xtop/D)^3-0.1006*(xtop/D)^2+0.3323*(xtop/D)+0.2042,0.6);
		Abot=min(0.0106*(xbot/D)^3-0.1006*(xbot/D)^2+0.3323*(xbot/D)+0.2042,0.6);
        % calculation of ultimate resistance pu and utilisation ratio:
        if strcmp(loads.static_cyclic,'static') % static case -> pu equal to pu_static
			disp_peak_top = (4.24*Atop^2+1.98*Atop-0.11)*y_tot(c,1)/4.1/Atop;
			disp_peak_bot = (4.24*Abot^2+1.98*Abot-0.11)*y_tot(c,2)/4.1/Abot;			

			if upy_top < disp_peak_top %before peak at (4.24A^2+1.98A-0.11)yc
                [output.UR(c,1), ~, output.p_UR(c,1), ~, output.pu_UR(c,1), ~] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,disp_peak_top,disp_peak_bot);
			else %after peak
				output.UR(c,1)    = 1; 	%For UR-plot UR=1
			end
			if upy_bot < disp_peak_bot %before peak at (4.24A^2+1.98A-0.11)yc
                [~,output.UR(c,2), ~, output.p_UR(c,2), ~,  output.pu_UR(c,2)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,disp_peak_top,disp_peak_bot);
			else
				output.UR(c,2)    = 1; 	%For UR-plot UR=1
			end
        elseif strcmp(loads.static_cyclic,'cyclic') % cyclic case -> pu equal to pu_cyclic at 0.45*yp
            if upy_top<0.45*y_tot(c,1) % before peak at 0.45*yp
                [output.UR(c,1), ~, output.p_UR(c,1), ~, output.pu_UR(c,1), ~] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,0.45*y_tot(c,1),0.45*y_tot(c,2));
            else % after peak at 0.45*yp
                output.UR(c,1)    = 1;       %For UR-plot UR=1
			end
			if upy_bot<0.45*y_tot(c,2) % before peak at 0.45*yp
                [~,output.UR(c,2), ~, output.p_UR(c,2), ~,  output.pu_UR(c,2)] = Calc_UR(element,pile,loads,c,upy_top,upy_bot,0.45*y_tot(c,1),0.45*y_tot(c,2));
			else
				output.UR(c,2)    = 1;       %For UR-plot UR=1
%             output.pu_UR(c,:)=output.p_UR(c,:);   %For integration ultimate resistance equal to mobilized resistance when above peak
            end
        end
    end
	    if output.hor_defl(c)/output.hor_defl(1)>0.0 && zero_defl==0    %Calculation of utilization ratio from seabed to rotation point
            if output.hor_defl(c+1)/output.hor_defl(1)>0.0
                output.Bhc=output.Bhc+(sum(output.p_UR(c,:)*0.5*(node.level(c)-node.level(c+1))));     %Summing up mobilized resistance
                output.PU=output.PU+(sum(output.pu_UR(c,:)*0.5*(node.level(c)-node.level(c+1))));      %Summing up ultimate resistance  
            else
                xrot=node.level(c)+(node.level(c+1)-node.level(c))*abs(output.hor_defl(c))/abs(output.hor_defl(c+1)-output.hor_defl(c));
                output.Bhc=output.Bhc+(sum(output.p_UR(c,:)*0.5*abs(node.level(c)-xrot)));     %Summing up mobilized resistance
                output.PU=output.PU+(sum(output.pu_UR(c,:)*0.5*(node.level(c)-xrot)));      %Summing up ultimate resistance  
            end
    else                            %The soil below the point of rotation is not taken into account
        output.Bhc=output.Bhc+0;    %integration of mobilized resistance
        output.PU=output.PU+0;        %integration of ultimate resistance. 
        zero_defl=1;    
    end
end
%% New addition
Defl_top=[1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.4 1.5 1.6 1.7 1.85 2.0 2.15 2.3 2.45 2.6 2.8 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.5 6.0 6.5 7.0 7.5 8.0 9.0 10 11 12 13 14 15.5 17 18.5 20 35 50 75 100 300 500 5000];
upy_top_ite=abs(ustep(1*3-2,settings.n_max));
output.PU2=zeros(1,length(Defl_top));
for c=1:nelem
	upy_top_hor(c)=ustep(c*3-2,settings.n_max)/ustep(1*3-2,settings.n_max);
end
[Y,I]=min(abs(upy_top_hor(1:ceil(length(upy_top_hor)*0.9))));
p_fit=polyfit(upy_top_hor(I-1:I+1)',node.level(I-1:I+1),1);
xrot=-p_fit(2);
for i=1:length(Defl_top)
	utop_pile(i)=upy_top_ite*Defl_top(i);		
	for c=1:nelem
		xtop    = node.level(c);    % [m] equivalent height of element top node
        xbot    = node.level(c+1);    % [m] equivalent height of element bottom node
		utop=utop_pile(i)-abs(utop_pile(i)/(xrot)*(xtop));
		ubot=utop_pile(i)-abs(utop_pile(i)/(xrot)*(xbot));
		[ksppynode_top  ksppynode_bot] = secspringstiff(element,pile,loads,[utop ubot],c); % secant stiffness values
        p_UR(1) = ksppynode_top*utop; % force corresponding to utop
        p_UR(2) = ksppynode_bot*ubot; % force corresponding to ubot
		if abs(xtop)<abs(xrot)
            if abs(xbot)<abs(xrot)
                output.PU2(1,i)=output.PU2(1,i)+(sum(p_UR(:)*0.5*(node.level(c)-node.level(c+1))));
            else
                output.PU2(1,i)=output.PU2(1,i)+(p_UR(1)*0.5*(node.level(c)+xrot));
            end
		end
    end
end
[maxPU2,indexPU2]=max(output.PU2);
if settings.ULS == 1
    figure
    plot(upy_top_ite,output.Bhc,'bo',utop_pile,output.PU2,'k.-',utop_pile(indexPU2),maxPU2,'ro')   
    xlabel('Pile deflection at seabed [m]')
    ylabel({'Integration of mobilised soil resistance';'from seabed to pivot point [kN]'})
    legend('Winkler model','Linear deflection profiles','R_{ult,d}','Location','Southeast')
    axis([0 ceil(max(utop_pile)/10)*10 0 ceil(max(output.PU2)/10)*10])
end
% figure(99)
% plot(Defl_top,output.PU2,'k',1,output.Bhc,'ro',Defl_top(1),output.PU2(1),'ko')
% xlabel('Multiplication factor on pile deflection at seabed [-]')
% ylabel('Integrated resistance to point of rotation [kN]')
% legend('Rigid pile','Non-rigid pile','Location','Southeast')
%% End of new addition

%% Calculation of global failure acc DIN 1054

B_hd                = output.Bhc * SF.B_hd; % design shear force incl psf (mobilized horizontal earth pressure integrated from mudline to rotational point)
R_ult_d        = max(output.PU2) / SF.R_ult; % design value of horizontal earth resistance
output.UR_global    = B_hd / R_ult_d; % utilisation ratio of global earth resistance
% disp(['mobilizied horizontal shear force: ',num2str(output.Bhc),' kN / shear force from turbine supplier: ',num2str(B_hd),' kN']);

    function [UR_top, UR_bot, p_value_top, p_value_bot, pult_value_top, pult_value_bot] = Calc_UR(element,pile,loads,c,disp_top,disp_bot, disp_ult_top,disp_ult_bot)	% calculate UR value for the curves

        [Ks_top,  Ks_bot] = secspringstiff(element,pile,loads,[disp_top disp_bot],c); % calculate secant stiffness (to be representative of current displacements)
        [Ks_u_top,  Ks_u_bot] = secspringstiff(element,pile,loads,[disp_ult_top disp_ult_bot],c); % calculate dummy secant stiffness (to be representative of failure)

        p_value_top = Ks_top*disp_top; % force corresponding to upy_top
        p_value_bot = Ks_bot*disp_bot; % force corresponding to upy_bot
        pult_value_top = Ks_u_top*disp_ult_top; % force corresponding to upy_top
        pult_value_bot = Ks_u_bot*disp_ult_bot; % force corresponding to upy_bot

        if Ks_u_top == 0
            UR_top = 1;
        else
            UR_top = p_value_top/pult_value_top;
        end
        if Ks_u_bot == 0
            UR_bot = 1;
        else 
            UR_bot = p_value_bot/pult_value_bot;
        end

    end
	
	    function [UR_top, UR_bot, URm_top, URm_bot, p_value_top, p_value_bot, pult_value_top, pult_value_bot, rot_value_top, rot_value_bot, rult_value_top, rult_value_bot]...
									= Calc_UR_PISA(element,pile,loads,c,disp_top,disp_bot,rotation_top,rotation_bot, disp_ult_top,disp_ult_bot)	% calculate UR value for the curves
									
		%secant stiffness of py spring
        [Ks_top,  Ks_bot] = secspringstiff(element,pile,loads,[disp_top disp_bot],c); % calculate secant stiffness (to be representative of current displacements)
        [Ks_u_top,  Ks_u_bot] = secspringstiff(element,pile,loads,[disp_ult_top disp_ult_bot],c); % calculate dummy secant stiffness (to be representative of failure)

		%secant stiffness of mom-rot springs
		ult_rot = 1e6; % high number to reach plateau.
		[Ksm_top Ksm_bot]  = secmomstiff(element,pile,[rotation_top rotation_bot],[disp_top disp_bot],Ks_top,Ks_bot,c);% mobilized secant stiffness values
		[Ksm_u_top Ksm_u_bot]  = secmomstiff(element,pile,[ult_rot ult_rot],[disp_top disp_bot],Ks_top,Ks_bot,c);% mobilized secant stiffness values
		
		
        p_value_top = Ks_top*disp_top; % force corresponding to upy_top
        p_value_bot = Ks_bot*disp_bot; % force corresponding to upy_bot
        pult_value_top = Ks_u_top*disp_ult_top; % force corresponding to upy_top
        pult_value_bot = Ks_u_bot*disp_ult_bot; % force corresponding to upy_bot
		
		rot_value_top = Ksm_top*rotation_top;
		rot_value_bot = Ksm_bot*rotation_bot;
		rult_value_top = Ksm_u_top*ult_rot; % force corresponding to ult_rot
		rult_value_bot = Ksm_u_bot*ult_rot; % force corresponding to ult_rot
		
        if Ks_u_top == 0
            UR_top = 1;
        else
            UR_top = p_value_top/pult_value_top;
        end
        if Ks_u_bot == 0
            UR_bot = 1;
        else 
            UR_bot = p_value_bot/pult_value_bot;
        end
		if Ksm_u_top == 0
            URm_top = 1;
        else
            URm_top = rot_value_top/rult_value_top;
        end
        if Ksm_u_bot == 0
            URm_bot = 1;
        else 
            URm_bot = rot_value_bot/rult_value_bot;
        end

    end

end