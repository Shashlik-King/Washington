function [fs qp A G] = skin_tip(element,reduction,pile,settings,loads,plots)
%% CALCULATION OF SKIN FRICTION
% API MODELS FOR PILES IN SAND AND CLAY 
% PROJECT SPECIFIC MODELS BASED ON PROJECT DESIGN BASIS
%--------------------------------------------------------------------------
% CHANGE LOG
% 2013.08.09    MUOE - PROGRAMMING, SEPARATING CALCULATION OF SKIN FRICTION
%                      FROM AXIAL ULS CALCULAITION
% 2014.07.15    MMOL - cleaning up simplifying
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% fs:           skin friction - inner, outer, total     [kPa]
% qp:           tip resistance                          [kN]
% A:            area (inside and outside)               [m^2]
% G:            weight of soil inside pile              [kN]
% plug_unplug:  used to determine state                 [-]
%--------------------------------------------------------------------------
%% Pre-allocating
%--------------------------------------------------------------------------
alpha   = zeros(length(element.model_axial)-1,1); % shaft friction factor [-]
qp      = zeros(length(element.model_axial)-1,1); % tip resistance [kN]

%--------------------------------------------------------------------------
%% Axial design parameters for cohesionless siliceous soil (API 1993)
%--------------------------------------------------------------------------
        
% % Axial_parameters = [15   47.8	8	 1.9 
% %                    20	67      12	 2.9
% %                    25	81.3    20  4.8
% %                    30	95.7    40  9.6
% %                    35   114.8   50  12]; % Column 1: delta_eff, 2:limiting skin resistance, 3: Nq factor, 4: Limiting tip resistance
% % for i = 1:element.nelem
% %     if element.phi(i) == 0
% %         element.delta_eff(i,:) = NaN;
% %     end
% %     if isnan(element.delta_eff(i)) == 0
% %         for j = 1:length(Axial_parameters)
% %             if element.delta_eff(i)     <= Axial_parameters(1,1)
% %                 element.limit_skin(i)   = Axial_parameters(1,2);
% %                 element.Nq(i)           = Axial_parameters(1,3);
% %                 element.limit_tip(i)    = Axial_parameters(1,4);
% %             elseif element.delta_eff(i) >= Axial_parameters(5,1)
% %                 element.limit_skin(i)   = Axial_parameters(5,2);
% %                 element.Nq(i)           = Axial_parameters(5,3);
% %                 element.limit_tip(i)    = Axial_parameters(5,4);
% %             else 
% %                 element.limit_skin(i)   = interp1(Axial_parameters(:,1),Axial_parameters(:,2),element.delta_eff(i));
% %                 element.Nq(i)           = interp1(Axial_parameters(:,1),Axial_parameters(:,3),element.delta_eff(i));
% %                 element.limit_tip(i)    = interp1(Axial_parameters(:,1),Axial_parameters(:,4),element.delta_eff(i));
% %             end
% %         end
% %     end
% % end

%--------------------------------------------------------------------------
%% Determination of skin (shaft) friction and tip resistance
%--------------------------------------------------------------------------
for i = 1:(element.nelem-1)
    G.soil(i) = pi*((pile.diameter/2-element.thickness(i))^2)*element.height(i)*element.gamma_eff(i); % weight of soil inside pile
    Nc = 9; % factor accounting for tip resistance
    for top_bottom = 1:2
        A.so(i,top_bottom) = pi*pile.diameter*element.height(i)/2; % outer skin friction area
        A.si(i,top_bottom) = pi*(pile.diameter-2*element.thickness(i))*element.height(i)/2; % inner skin friction area
        
        %% API clay
        if strcmp(element.model_axial(i),'API clay')
            % ALPHA-METHOD
            if strcmp(settings.clay_type,'alpha')
                if element.cu(i,top_bottom)/element.sigma_v_eff(i,top_bottom) <= 1.0
                    alpha(i) = min(element.limit_alpha(i),1/(2*(element.cu(i,top_bottom)/element.sigma_v_eff(i,top_bottom))^(1/2)));
                else
                    alpha(i) = min(element.limit_alpha(i),1/(2*(element.cu(i,top_bottom)/element.sigma_v_eff(i,top_bottom))^(1/4)));
                end
                if isnan(alpha(i));
                    alpha(i) = 0;
                end
                fs.o(i,top_bottom) = alpha(i)*element.cu(i,top_bottom); % outer skin friction on the element in the clay case
                % BETA-METHOD
            elseif strcmp(settings.clay_type,'beta')
                fs.o(i,top_bottom) = min(element.c_eff(i)+element.K0(i)*element.sigma_v_eff(i,top_bottom)*tand(element.delta_eff(i),element.limit_skin(i)));
            else
                error('Choose ''alpha'' or ''beta'' for settings.clay_type') % Error message
            end
            % COMMON FOR ALPHA AND BETA-METHODS FOR CLAY
            qp(i,top_bottom) = min(Nc*element.cu(i,top_bottom),element.limit_tip(i)); % tip resistance on the element in the clay case (only makes sense for last (bottom) element)
            
            %% API SAND
        elseif strcmp(element.model_axial(i),'API sand')
            fs.o(i,top_bottom) = min(element.K0(i)*element.sigma_v_eff(i,top_bottom)*tand(element.delta_eff(i)),element.limit_skin(i)); % outer skin friction on the element (in the static sand case)
            
            qp(i,top_bottom) = min(element.Nq(i)*element.sigma_v_eff(i,top_bottom),element.limit_tip(i)); % tip resistance on the element in the sand case (only makes sense for last (bottom) element)
            %% ZERO SOIL
        elseif strcmp(element.model_axial(i),'Zero soil')
            fs.o(i,top_bottom) = 0;
            fs.i(i,top_bottom) = 0;
            qp(i,top_bottom) = 0;
        else
            error('Please select an axial model that has been implemented')
        end
%% Common for sand and clay        
        
    fs.o(i,top_bottom) = element.degradation_tz_t(i)*fs.o(i,top_bottom); % outer skin friction on reduced due to degradation    
    fs.i(i,top_bottom) = reduction.skin_inner*fs.o(i,top_bottom); % inner skin friction on the element             
    
    end
end

%% Account for torsional loading when determining axial capacity


if settings.torsion == 1 && loads.Mz > 0 % then account for torsional loading
    
    if plots.res_vs_pilelength         == 1
    % calculate outer skin surface
    
    
     r_outer = pile.diameter/2; % outer pile radius
   fs.T=zeros(length(fs.o(:,1)),1);
    
    
     for i = 1:element.nelem-1
        for top_bottom = 1:2
            if fs.o(i,top_bottom) > 0 % take only outer skin area into account, where torsional load can be transferred to subsoil
                fs.A_outer(i,top_bottom) = pi*pile.diameter*element.height(i)/2; % outer skin friction area
                
            else
                fs.A_outer(i,top_bottom) = 0; % outer skin friction area
            end
        end
     end
     for i = 1:element.nelem-1
            if i==1
            fs.A_outer_total(i,1) = (fs.A_outer(i,1))+(fs.A_outer(i,2));
            else
                fs.A_outer_total(i,1) =fs.A_outer_total(i-1,1)+ (fs.A_outer(i,1))+(fs.A_outer(i,2));
            end
        
     end
    %fs.A_outer_total = sum(sum(fs.A_outer)); % outer skin surface area of the pile, where torsional load can be transferred to subsoil
   
   for i = 1:element.nelem-1
        if fs.A_outer_total(i,1)==0
            fs.T(i)=0;  
        else
         fs.T(i) = abs(loads.Mz) / fs.A_outer_total(i) / r_outer; % uniform distributed skin friction due to torsional loading
           
        end
   end
   
    
    for i = 1:element.nelem-1
  
        for top_bottom = 1:2  
            if fs.o(i,top_bottom) > 0
            if fs.o(i,top_bottom) >= fs.T(pile.axial_elem) % resulting skin friction has to be > 0
                fs.o(i,top_bottom)  = fs.o(i,top_bottom) * sqrt(1 - (fs.T(pile.axial_elem)/ fs.o(i,top_bottom))^2);
            end
            end
                
        
        end
      % f_sT
    end
    
    else
    for i = 1:element.nelem-1
        for top_bottom = 1:2
            if fs.o(i,top_bottom) > 0 % take only outer skin area into account, where torsional load can be transferred to subsoil
                fs.A_outer(i,top_bottom) = pi*pile.diameter*element.height(i)/2; % outer skin friction area
            else
                fs.A_outer(i,top_bottom) = 0; % outer skin friction area
            end
        end
    end
    fs.A_outer_total = sum(sum(fs.A_outer)); % outer skin surface area of the pile, where torsional load can be transferred to subsoil
    r_outer = pile.diameter/2; % outer pile radius
    fs.T = abs(loads.Mz) / fs.A_outer_total / r_outer; % uniform distributed skin friction due to torsional loading
    for i = 1:element.nelem-1
        for top_bottom = 1:2  
            if fs.o(i,top_bottom) >= fs.T % resulting skin friction has to be > 0
                fs.o(i,top_bottom)  = fs.o(i,top_bottom) * sqrt(1 - (fs.T / fs.o(i,top_bottom))^2);
            end
        end
    end
    end
elseif settings.torsion == 0 % do not account for torsional loading
else
    disp('Influence of torsional load on axial capacity not defined');
end