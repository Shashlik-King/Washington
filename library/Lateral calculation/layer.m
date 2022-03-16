function [element] = layer(element,node,pile,settings,loads)
%--------------------------------------------------------------------------
% PURPOSE
% Compute parameters employed in order to take effects of layering into
% account. The method of Georgiadis (1983) is used.
%
% INPUT:  springinput   : cf. mainfile
%         PUtot         : Total horz. capacity in the former element
%         Coord         : Global coordinate matrix
%         i             : Counter referring to element number
%		  D   			: pile outer diameter
%         L  			: pile length
%         nelem  		: number of elements
%		  settings		: structure that contains the settings for calculations
%
% OUTPUT: pu      : Horisontal capacity at the top and the bottom of the
%                   segment [kN/m]  [pu_top pu_bot]. For PISA models this
%					is normalized [-] Size=(nelem,2)
%         mu      : For PISA models normalized ultimate moment capacity at the top and the bottom of the
%                   segment [-]  [mu_top mu_bot]. Size=(nelem,2)
%         HBu     : For PISA models normalized ultimate base shear capacity at the top and the bottom of the
%                   segment [-]  [MBu_top MBu_bot]. Size=(nelem,2). In the
%                   final calculations only the last node value is used.
%         MBu     : For PISA models normalized ultimate base moment capacity at the top and the bottom of the
%                   segment [-]  [MBu_top MBu_bot]. Size=(nelem,2).  In the
%                   final calculations only the last node value is used.
%         zr      : Transition depth [m] determined based on homogeneous
%                   soil.
%         heqv    : Equivalent height for top and bottom of the segment
%                   [m]. [heqv_top heqv_bot]
%
% Log:
% 20.05.2008    AHAU    Programming - Sand and soft clay
% 01.08.2008    AHAU    Free standing length
% 2010.11.04    MMOL    Integrating into calculation routine + stiff clay
% 2011.04.28    JALY    Streamlining code, new interpretation of Georgiadis
% 21.09.2017    EVVA    Added PISA models
% 2018.04.04	DATY	On/Off switch for Georgiadis approach
%--------------------------------------------------------------------------

nelem=length(element.model_py);

%Preallocation
element.pu=NaN(nelem,2);
element.mu=NaN(nelem,2);
element.HBu=NaN(nelem,2);
element.MBu=NaN(nelem,2);
element.heqv=NaN(nelem,2);
element.hr=NaN(nelem,2);

%General values
D = pile.diameter;          % Outer diameter [m]
L = pile.L;      % Pile length [m]
PU = 0;                     % accumulated PU, initial value

display=1; 

for i=1:nelem %through all elements
    %General values
    gamma = element.gamma_eff(i);   % Effective unit weight [kN/m^3]
    su_top   = element.cu(i,1);      % Undrained shear strength in the top node[kPa]
    su_bot   = element.cu(i,2);   % Undrained shear strength in the bottom node[kPa]
    J       = element.J(i);              % Dimensionless empirical constant
    phi     = element.phi(i);         % Triaxial design value of friction angle [degrees]
    h_segment = node.level(i)-node.level(i+1); %length of pile segment [m]

    %% ------------------------------------------------------------------------
    % Sand - API(1993)
    % -------------------------------------------------------------------------
    if strcmp(element.model_py(i),'API sand') || strcmp(element.model_py(i),'Kirsch sand') || strcmp(element.model_py(i),'Kallehave sand')  || strcmp(element.model_py(i),'Kirsch sand 2015')
        %General values
        K=0.4;
        Beta=45+phi/2;
        rad=pi/180;
        
        %Factors according to Mosher and Dawkins 2008.
        C1 = (K*tan(phi*rad)*sin(Beta*rad))/(tan((Beta-phi)*rad)*cos((phi/2)*rad))+((tan(Beta*rad))^2*tan((phi/2)*rad))/(tan((Beta-phi)*rad))+K*tan(Beta*rad)*(tan(phi*rad)*sin(Beta*rad)-tan((phi/2)*rad));
        C2 = tan(Beta*rad)/tan((Beta-phi)*rad)-(tan((45-phi/2)*rad))^2;
        C3 = K*tan(phi*rad)*(tan(Beta*rad))^4+(tan((45-phi/2)*rad))^2*((tan(Beta*rad))^8-1);

        for top_bottom=1:2 %looping through top and bottom node of each layer

            
            % -------- Determination of equivalent height, according to Georgiadis (1983) -------
            % Determine the imaginary height, heqv, the amount of identical
            % material as the actual layer on top of the actual node that
            % would provide the same accumulated resistance as the actual PU_top
            
            %Transition depth, if the actual layer was uniform soil
            hr = D*(C3-C2)/C1;  
            % Integrated resistance for shallow failure from h=0 to h=hr;
            R_shallow = gamma*C1/3*hr^3 + gamma*C2*D/2*hr^2;
            
            if R_shallow > PU % shallow failure only down to actual node

                p  = [gamma*C1/3 gamma*C2*D/2 0 -PU]; %Factors in polynomium that
                % is solved to determine heqv based on expreesions for moderate depth.
                r = roots(p);       % Roots for the polynomium described by p.
                
                rr=[];
                for k=1:3
                    if isreal(r(k))
                        rr=[rr r(k)]; %extract all real results
                    end
                end
                heqv = min(rr(rr>=0)); %finding the minimum positive root
                
            else % both shallow and deep failure above actual node (ie at actual node is deep failure)
                % heqv is determined as the solution to PU = R_shallow + C3 * D * gamma * (heqv^2-hr^2) / 2
                 heqv = sqrt(hr^2 + 2*(PU-R_shallow)/(C3*D*gamma));
            end
            
            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 0     % Georgiadis is turned off
                heqv    = -element.level(i,top_bottom);
            end
            
            %"Old" interpretation of georgardis, used on LA
%             if top_bottom==1 %top node
%                 sigma = gamma*heqv;
%             elseif top_bottom==2 %bottom node
%                 sigma = gamma*heqv+h_segment;
%             end
%             if display; disp('Old Georgiadis interpretation'); display=0; end
            
%             test(i,:)=[gamma*heqv node.sigma_v_eff(i)];
            %"New" interpretation of georgardis
            if top_bottom==1 %top node
                sigma = node.sigma_v_eff(i);
            elseif top_bottom==2 %bottom node
                sigma = node.sigma_v_eff(i+1);
            end
            heqv = sigma/gamma;
%             if display; disp('New Georgiadis interpretation'); display=0; end


            % -------- Determination of hr, taking Georgiadis' principle into consideration
            % This means that the layer is assumed to extent infinitely
            % below the actual node, but the actual stress at the node
            % level is taken into consideration together with the
            % determined equivalent depth
            % An auxiliary height, haux, is introduced. It is the distance
            % from the considered node (with equivalent depth heqv) to hr depth.
            % Thus, hr = heqv + haux.
            % Thus, the equation to solve may be expressed as:
            % (C1*hr+C2*D) = C3*D
            % If haux is negative or not real.
            
            haux = (C3*D-C1*heqv-C2*D)/C1;
            hr = heqv + haux;
                        
            % -------- Determination of pu
            if hr>heqv %Shallow failure
                pu = (C1*heqv+C2*D)*sigma;
            else %Deep failure
                pu = C3*D*sigma;
            end
            
           
            % -------- Saving the computed values
            element.pu(i,top_bottom)=pu;
            element.heqv(i,top_bottom)=heqv;
            % hr is not saved as it is not needed in
            % the further calculations
            
            % -------- Integrating pu down each element
            if top_bottom==1 %if this is an evaluation of the top spring
                if hr>heqv+h_segment %if only shallow failure exists down to the bottom spring
                   delta_PU = 1/3*C1*gamma*((heqv+h_segment)^3-heqv^3) + 0.5*(C1*(sigma-gamma*heqv)+C2*D*gamma)*((heqv+h_segment)^2-heqv^2)+C2*D*(sigma-gamma*heqv)*h_segment; 
                elseif hr<heqv %if deep failure both at top and bottom spring
                    delta_PU = 0.5*C3*D*gamma*((heqv+h_segment)^2-heqv^2)+C3*D*(sigma-gamma*heqv)*h_segment;
                else %combined failure
                    h_upper   = hr - heqv;
                    delta_PU1 = 1/3*C1*gamma*((heqv+h_upper)^3-heqv^3) + 0.5*(C1*(sigma-gamma*heqv)+C2*D*gamma)*((heqv+h_upper)^2-heqv^2) + C2*D*(sigma-gamma*heqv)*h_upper; %shallow part
                    delta_PU2 = 0.5*C3*D*gamma*((heqv+h_segment)^2-hr^2)+C3*D*(sigma-gamma*heqv)*(heqv+h_segment-hr); %deep part
                    delta_PU  = delta_PU1 + delta_PU2;
                end
                PU=PU+delta_PU;
            end

        end
        
        %% -------------------------------------------------------------------------
        % Soft clay - API(1993)
        % -------------------------------------------------------------------------
    elseif strcmp(element.model_py(i),'API clay') || strcmp(element.model_py(i),'Stiff clay w/o free water') || strcmp(element.model_py(i),'Kirsch soft clay') || strcmp(element.model_py(i),'Kirsch stiff clay')
        for top_bottom=1:2 %looping through top and bottom node of each layer
            
            % Correcting errors causing numerical instability
            if su_top == 0
                su_top = 1;
            elseif su_bot == 0
                su_bot = 1;
            end
            
            % Setting parameters
            if top_bottom == 1 %top node
                su = su_top;
            elseif top_bottom == 2 %bottom node
                su = su_bot;
            end
            
            % -------- Determination of equivalent height, according to Georgiadis (1983) -------
            % Determine the imaginary height, heqv, the amount of identical
            % material as the actual layer on top of the actual node that
            % would provide the same accumulated resistance as the actual PU_top
            
            hr =(6*su*D)/(gamma*D+J*su); %Transition depth, if the actual layer was uniform soil
            
            % Integrated resistance for shallow failure from h=0 to h=hr;
            R_shallow = 0.5*(gamma*D+J*su)*hr^2 + 3*su*D*hr;
            
            if R_shallow > PU % shallow failure only down to actual node
                
                p  = [0.5*(gamma*D+J*su_top) 3*su_top*D -PU]; %Factors in polynomium that
                % is solved to determine heqv based on expreesions for moderate depth.
                
                r = roots(p);       % Roots for the polynomium described by p.
                
                if isreal(r)==0
                    error('Equivalent length determination is wrong: Clay')
                end
                heqv = min(r(r>=0)); %finding the minimum positive root
                
            else % both shallow and deep failure above actual node (ie at actual node is deep failure)
                % heqv is determined as the solution to PU = R_shallow + 9 * su * D *(heqv-hr)
                heqv = (PU - R_shallow)/(9 * su * D) + hr;
            end

            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 0     % Georgiadis is turned off
                heqv    = -element.level(i,top_bottom);
            end
 
             %"Old" interpretation of georgardis, used on LA
%             sigma = gamma*heqv;
%             if display; disp('Old Georgiadis interpretation'); display=0; end
            
%             
            %"New" interpretation of georgardis
            if top_bottom==1 %top node
                sigma = node.sigma_v_eff(i);
            elseif top_bottom==2 %bottom node
                sigma = node.sigma_v_eff(i+1);
            end
            heqv = sigma/gamma; %GMME F2
%             if display; disp('New Georgiadis interpretation'); display=0; end

            % -------- Determination of hr, taking Georgiadis' principle into consideration
            % This means that the layer is assumed to extent infinitely
            % below the actual node, but the actual stress at the node
            % level is taken into consideration together with the
            % determined equivalent depth
            % An auxiliary height, haux, is introduced. It is the distance
            % from the considered node (with equivalent depth heqv) to hr depth.
            % Thus, hr = heqv + haux.
            % Thus, the equation to solve may be expressed as:
            % 3 + (sigma_node + gamma_eff * haux)/su + J * hr / D = 9
            % If haux is negative, the failure is deep (since hr < heqv).
            
            haux = (6-sigma/su-J*heqv/D)/(gamma/su+J/D);
            hr = heqv + haux;
            
            % -------- Determination of pu
            if hr>heqv %Shallow failure
                pu = (3*su + sigma)*D + J*su*heqv;
            else %Deep failure
                pu = 9 * su * D;
            end
            
            % -------- Saving the computed values
            element.pu(i,top_bottom)=pu;
            element.heqv(i,top_bottom)=heqv;
            element.hr(i,top_bottom)=hr;
            
            % -------- Integrating pu down each element
            if top_bottom==1 %if this is an evaluation of the top spring
                if hr>heqv+h_segment %if only shallow failure exists down to the bottom spring
                    delta_PU=0.5*(gamma*D+J*su)*((heqv+h_segment)^2-heqv^2) + (3*su+sigma-gamma*heqv)*D*h_segment;
                elseif hr<heqv %if deep failure both at top and bottom spring
                    delta_PU = 9 * su * D * h_segment;
                else %combined failure
                    h_upper = hr - heqv; 
                    h_lower = heqv + h_segment - hr;
                    delta_PU1 = 0.5*(gamma*D+J*su)*((heqv+h_upper)^2-heqv^2) + (3*su+sigma-gamma*heqv)*D*h_upper; %shallow part
                    delta_PU2 = 9 * su * D * h_lower; %deep part
                    delta_PU  = delta_PU1 + delta_PU2;
                end
                PU=PU+delta_PU;
            end
        end

        %% -------------------------------------------------------------------------
        % Jeanjean (2009) - Clay py curve
        % -------------------------------------------------------------------------
    elseif strcmp(element.model_py(i),'Jeanjean clay')
        for top_bottom=1:2 %looping through top and bottom node of each layer
            
            % Correcting errors causing numerical instability
            if su_top == 0
                su_top = 1;
            elseif su_bot == 0
                su_bot = 1;
            end
            
            % Setting parameters
            if top_bottom == 1 %top node
                su = su_top;
            elseif top_bottom == 2 %bottom node
                su = su_bot;
            end
            
            % -------- Determination of equivalent height, according to Georgiadis (1983) -------
            % Determine the imaginary height, heqv, the amount of identical
            % material as the actual layer on top of the actual node that
            % would provide the same accumulated resistance as the actual PU_top
            
            hr =(6*su*D)/(gamma*D+J*su); %Transition depth, if the actual layer was uniform soil
            
            % Integrated resistance for shallow failure from h=0 to h=hr;
            R_shallow = 0.5*(gamma*D+J*su)*hr^2 + 3*su*D*hr;
            
            if R_shallow > PU % shallow failure only down to actual node
                
                p  = [0.5*(gamma*D+J*su_top) 3*su_top*D -PU]; %Factors in polynomium that
                % is solved to determine heqv based on expreesions for moderate depth.
                
                r = roots(p);       % Roots for the polynomium described by p.
                
                if isreal(r)==0
                    error('Equivalent length determination is wrong: Clay')
                end
                heqv = min(r(r>=0)); %finding the minimum positive root
                
            else % both shallow and deep failure above actual node (ie at actual node is deep failure)
                % heqv is determined as the solution to PU = R_shallow + 9 * su * D *(heqv-hr)
                heqv = (PU - R_shallow)/(9 * su * D) + hr;
            end

            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 0     % Georgiadis is turned off
                heqv    = -element.level(i,top_bottom);
            end
 
             %"Old" interpretation of georgardis, used on LA
%             sigma = gamma*heqv;
%             if display; disp('Old Georgiadis interpretation'); display=0; end
            
%             
            %"New" interpretation of georgardis
            if top_bottom==1 %top node
                sigma = node.sigma_v_eff(i);
            elseif top_bottom==2 %bottom node
                sigma = node.sigma_v_eff(i+1);
            end
            heqv = sigma/gamma; %GMME F2
%             if display; disp('New Georgiadis interpretation'); display=0; end

            % -------- Determination of hr, taking Georgiadis' principle into consideration
            % This means that the layer is assumed to extent infinitely
            % below the actual node, but the actual stress at the node
            % level is taken into consideration together with the
            % determined equivalent depth
            % An auxiliary height, haux, is introduced. It is the distance
            % from the considered node (with equivalent depth heqv) to hr depth.
            % Thus, hr = heqv + haux.
            % Thus, the equation to solve may be expressed as:
            % 3 + (sigma_node + gamma_eff * haux)/su + J * hr / D = 9
            % If haux is negative, the failure is deep (since hr < heqv).
            
            haux = (6-sigma/su-J*heqv/D)/(gamma/su+J/D);
            hr = heqv + haux;
            
            % -------- Determination of pu according to Jeanjean (2009)

            pu = element.Np(i)*su*D;

            
            % -------- Saving the computed values
            element.pu(i,top_bottom)=pu;
            element.heqv(i,top_bottom)=heqv;
            element.hr(i,top_bottom)=hr;
            
            % -------- Integrating pu down each element
            if top_bottom==1 %if this is an evaluation of the top spring
                if hr>heqv+h_segment %if only shallow failure exists down to the bottom spring
                    delta_PU=0.5*(gamma*D+J*su)*((heqv+h_segment)^2-heqv^2) + (3*su+sigma-gamma*heqv)*D*h_segment;
                elseif hr<heqv %if deep failure both at top and bottom spring
                    delta_PU = 9 * su * D * h_segment;
                else %combined failure
                    h_upper = hr - heqv; 
                    h_lower = heqv + h_segment - hr;
                    delta_PU1 = 0.5*(gamma*D+J*su)*((heqv+h_upper)^2-heqv^2) + (3*su+sigma-gamma*heqv)*D*h_upper; %shallow part
                    delta_PU2 = 9 * su * D * h_lower; %deep part
                    delta_PU  = delta_PU1 + delta_PU2;
                end
                PU=PU+delta_PU;
            end
        end
                
        %% -------------------------------------------------------------------------
        % Stiff clay - Reese et al. (1975)
        % -------------------------------------------------------------------------
    elseif strcmp(element.model_py(i),'Reese stiff clay')
        for top_bottom=1:2 %looping through top and bottom node of each layer
            %Setting parameters
            if top_bottom==1 %top node
                su = su_top;
            elseif top_bottom==2 %bottom node
                su = su_bot;
            end
            
            % -------- Determination of equivalent height, according to Georgiadis (1983) -------
            % Determine the imaginary height, heqv, the amount of identical
            % material as the actual layer on top of the actual node that
            % would provide the same accumulated resistance as the actual PU_top
            
            ca = su; %(su_top*heqvref+su*h)/(heqvref+h); %Average undrained soil shear strength - to be checked!!!!
            hr = ((11*su-2*ca)*D)/(gamma*D+2.83*ca); % Transition depth [m] determined
            %by equating expressions for deep and shallow depth, if the actual layer was uniform soil
            
            % Integrated resistance for shallow failure from h=0 to h=hr;
            R_shallow = 0.5*(gamma*D+2.83*ca)*hr^2 + 2*ca*D*hr;
            
            if R_shallow > PU % shallow failure only down to actual node
                p  = [0.5*(gamma*D+2.83*ca) 2*ca*D -PU]; %Factors in polynomium that
                % is solved to determine heqv based on expreesions for moderate depth.
                
                r = roots(p);       % Roots for the polynomium described by p.
                
                if isreal(r)==0
                    error('Equivalent length determination is wrong: Stiff Clay')
                end
                heqv = min(r(r>=0)); %finding the minimum positive root
                
            else % both shallow and deep failure above actual node (ie at actual node is deep failure)
                % heqv is determined as the solution to PU = R_shallow + 11 * su * D *(heqv-hr)
                heqv = (PU - R_shallow)/(11 * su * D) + hr;
            end

            
            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 0     % Georgiadis is turned off
                heqv    = -element.level(i,top_bottom);
            end
 
             %"Old" interpretation of georgardis, used on LA
            sigma = gamma*heqv;
%             if display; disp('Old Georgiadis interpretation'); display=0; end
            
%             
%             %"New" interpretation of georgardis
%             if top_bottom==1 %top node
%                 sigma = node.sigma_v_eff(i);
%             elseif top_bottom==2 %bottom node
%                 sigma = node.sigma_v_eff(i+1);
%             end
%             if display; disp('New Georgiadis interpretation'); display=0; end

            % -------- Determination of hr, taking Georgiadis' principle into consideration
            % This means that the layer is assumed to extent infinitely
            % below the actual node, but the actual stress at the node
            % level is taken into consideration together with the
            % determined equivalent depth
            % An auxiliary height, haux, is introduced. It is the distance
            % from the considered node (with equivalent depth heqv) to hr depth.
            % Thus, hr = heqv + haux.
            % Thus, the equation to solve may be expressed as:
            % 3*ca + sigma_node + gamma_eff * haux + 2.83 * ca * hr / D = 11*cu
            % If haux is negative, the failure is deep (since hr < heqv).
            
            haux = (11*su-2*ca-sigma-2.83*ca*heqv/D)/(gamma+2.83*ca/D);
            hr = heqv + haux;
            
            % -------- Determination of pu
            if hr>heqv %Shallow failure
                pu = (2*ca + sigma)*D + 2.83*ca*heqv;
            else %Deep failure
                pu = 11 * su * D;
            end
            
            % -------- Saving the computed values
            element.pu(i,top_bottom)=pu;
            element.heqv(i,top_bottom)=heqv;
            % hr is not saved as it is not needed in
            % the further calculations
            
            % -------- Integrating pu down each element
            if top_bottom==1 %if this is an evaluation of the top spring
                if hr>heqv+h_segment %if only shallow failure exists down to the bottom spring
                    delta_PU = 0.5*(gamma*D+2.83*ca)*((heqv+h_segment)^2-heqv^2) + (2*ca+sigma-gamma*heqv)*D*h_segment;
                elseif hr<heqv %if deep failure both at top and bottom spring
                    delta_PU = 11 * su * D * h_segment;
                else %combined failure
                    h_upper   = hr - heqv; h_lower = heqv + h_segment - hr;
                    delta_PU1 = 0.5*(gamma*D+2.83*ca)*((heqv+h_upper)^2-heqv^2) + (2*ca+sigma-gamma*heqv)*D*h_upper;%shallow part
                    delta_PU2 = 11 * su * D * h_lower; %deep part
                    delta_PU  = delta_PU1 + delta_PU2;
                end
                PU=PU+delta_PU;
            end
        end
        
        %% -------------------------------------------------------------------------
        % Weak Rock - Reese
        % -------------------------------------------------------------------------
    elseif strcmp(element.model_py(i),'Modified Weak rock')
        for top_bottom = 1:2
            % -------- Determination of pu
            q_ur = element.q_ur(i,top_bottom); % [kPa]
            alpha_r = 1-0.67*element.RQD(i);  
            first_weak_rock = find(strcmp(element.model_py,'Modified Weak rock'),1); % find first time weak rokc py-model is applied
            pu = min(alpha_r*q_ur*D*(1+1.4*(abs(element.level(i,top_bottom)-element.level(first_weak_rock,top_bottom)))/D),5.2*alpha_r*q_ur*D);
            % always choosing lowest value corresponds to choosing shallow
            % capacity for shallow depths and deep capacity for deep depths
            
            % -------- Saving the computed values
            element.pu(i,top_bottom) = pu;
            element.heqv(i,top_bottom) = -element.level(i,top_bottom);
            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 1     % Georgiadis is turned on
                error(['Georgiadis approach is not implemented for ',element.model_py{i},' model. Please turn off Georgiadis approach.']);
            end
        end
        %% -------------------------------------------------------------------------
        % PISA soils
        % -------------------------------------------------------------------------
    elseif settings.PISA_database == 1 
        for top_bottom = 1:2
            % -------- Determination of normalized ultimate values.
            pu = element.PISA_param(i,1+top_bottom); 
            mu =  element.PISA_param(i,8+top_bottom);
			HBu = element.PISA_param(i,16);  
            MBu = element.PISA_param(i,20);
            % -------- Saving the computed values
            
            element.heqv(i,top_bottom) = -element.level(i,top_bottom);
            
            if strcmp(loads.static_cyclic,'cyclic')
                if strcmp(loads.A,'API')
                    cyclic    =   0.9; % acc. API
                    Cyclic_Degradation=cyclic/max(0.9 , 3.0 - 0.8*element.heqv(i,top_bottom)/pile.diameter); % Ratio between Cyclic and static in API
                    Cyclic_Degradation_clay=0.72;
                elseif strcmp(loads.A,'TUHH')
                    cyclic =   min(0.343*element.heqv(i,top_bottom)/pile.diameter,0.9);  % for more than 100 load cycles, acc. EA-Piles (p.443, Equ. D3.6)
                    Cyclic_Degradation=cyclic/max(0.9 , 3.0 - 0.8*element.heqv(i,top_bottom)/pile.diameter); % Ratio between Cyclic and static in API
                    Cyclic_Degradation_clay=0.72;
                elseif strcmp(loads.A,'P_NGI')
                    Cyclic_Degradation = element.cyclic_ult(i);
                    Cyclic_Degradation_clay = element.cyclic_ult(i);
                end
			elseif strcmp(loads.static_cyclic,'static')
            Cyclic_Degradation=1;
			Cyclic_Degradation_clay=1;
            end 
                
            element = A_factor(element,pile,loads,i);
            
            if settings.psf_switch
                if strcmp(element.type{i},'Sand')
                    element.pu(i,top_bottom) = (pu*Cyclic_Degradation)/settings.psf_drained;
                    %warning('PSF to pu in Sand is applied')
					element.HBu(i,top_bottom) = HBu*Cyclic_Degradation/settings.psf_drained;
					element.MBu(i,top_bottom) = MBu*Cyclic_Degradation/settings.psf_drained;
                    element.mu(i,top_bottom) = mu;
                elseif strcmp(element.type{i},'Clay')  
                    element.pu(i,top_bottom) = pu*Cyclic_Degradation_clay/settings.psf_undrained;
                    element.mu(i,top_bottom) = mu*Cyclic_Degradation_clay/settings.psf_undrained;
					element.HBu(i,top_bottom) = HBu*Cyclic_Degradation_clay/settings.psf_undrained;
					element.MBu(i,top_bottom) = MBu*Cyclic_Degradation_clay/settings.psf_undrained;
                end
            else
                if strcmp(element.type{i},'Sand')
                    element.pu(i,top_bottom) = (pu*Cyclic_Degradation);
					element.HBu(i,top_bottom) = HBu*Cyclic_Degradation;
					element.MBu(i,top_bottom) = MBu*Cyclic_Degradation;
					
%                     warning('PSF to pu in Sand is applied')
                    element.mu(i,top_bottom) = mu;
                elseif strcmp(element.type{i},'Clay')
                    element.pu(i,top_bottom) = pu*Cyclic_Degradation_clay;
                    element.mu(i,top_bottom) = mu*Cyclic_Degradation_clay;
					element.HBu(i,top_bottom) = HBu*Cyclic_Degradation_clay;
					element.MBu(i,top_bottom) = MBu*Cyclic_Degradation_clay;
                end
            end
  
            
            % switch in order to turn on/off Georgiadis approach
            if settings.Georgiadis == 1     % Georgiadis is turn on
                error(['Georgiadis approach is not implemented for ',element.model_py{i},' model. Please turn off Georgiadis approach.']);
            end
        end 

        %% -------------------------------------------------------------------------
        % Free standing length, zero soil
        % -------------------------------------------------------------------------
    elseif strcmp(element.model_py(i),'Zero soil')
               element.pu(i,:)= [0 0];
        %% -------------------------------------------------------------------------
        % Not recognised soil model
        % -------------------------------------------------------------------------
    else
        error('The specified soil model is not supported')
    end     
end