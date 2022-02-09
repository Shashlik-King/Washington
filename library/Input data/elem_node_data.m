% Function that establishes input data to the winkler-py function
% From the basis of the soil profile and the pile information, a set of
% data vectors are made for each beam-element in the pile and for each
% node.
%
% Units: kN, m, s, kPa
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Log of changes------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Date            Initials        Change
%2011.04.27      JALY            Programming
%2016.04.13      EVVA            ep [E A I D_eq Gm ksf]
%2017.09.21		 EVVA			 Added extra element for base shear/moment calculation
%2019.12.02		 FKMV			 PISA soil springs calibration values for
%cu and gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element,node, pile]=elem_node_data(pile,soil,scour,settings,data,plots,c,Database)

%Overburden reduction and scour level
ORD_level    = soil.toplevel(1)-scour.ORD;
scour_level  = soil.toplevel(1)-scour.local;

if pile.head<soil.toplevel(1)
    error('Pile head is specified below the either seabed or cross section data.')
end

if plots.res_vs_pilelength == 1 && settings.axial_loading == 1
    if settings.interimgeometry        == 0
        pile.length_start_or           = pile.length;
        pile.length_start_or           = pile.length_start_or+pile.extra_L; % [m] original embedment length (pile_tip is defined in mLAT) %%%%%%%%%%%%%%%%%%%%%%%%%
        pile.length_start = pile.length_start_or ;   
        pile.length_end             = pile.length_start; % [m] for general calculations only one pile embedment length is considered
        pile.length_inc             = 1; % [m] is not used for general calculations

        pile_tip=-scour.water_depth-pile.length_start_or;
        pile.length                     = (pile.length_start:pile.length_inc:pile.length_end);
        pile.tip = pile_tip; 
        pile.L = pile.length;
    end
end

pile.L = pile.length(c); % Selecting pile length
disp(['Pile length: ',num2str(pile.L),' m'])
pile.toe = pile.head - pile.L; % pile toe level


 % Selecting pile length
% if sum(pile.toe==soil.toplevel)
%     warning('Routine:PileToeAtBoundary','Pile toe level is located at a soil layer boundary. The layer directly below the pile toe is disregarded.')
% end

% Dividing the pile into parts, where either a change in pile cross section
% or soil profile is present.
toplevels=sort(unique([scour_level ; ORD_level ; soil.toplevel ; pile.cross_section.toplevel ]),1,'descend'); %all layer boundaries


% sometimes the sort+unique does not work if two levels have a minimum difference (e.g. 45.00000001 and 44.999999999)so the below fix is added%
toplevels(abs(diff(toplevels))<0.0001)=[];



%The actual needed levels for the pile
j=find(toplevels<pile.toe,1);
if isempty(j) 
    j=length(toplevels)+1; 
end %if the pile toe is the lower specified level
k=find(toplevels>pile.head,1,'last');
if isempty(k) 
    k = 0; 
end %if the pile head is the upper specified level
section_levels=sort(unique([pile.head ; toplevels(k+1:j-1) ; pile.toe]),'descend');
% section_levels2=sort(unique([pile.head ; toplevels(k+1:j-1) ; pile.toe; (pile.toe-0.1)]),'descend'); %(pile.toe-0.1) is an extra section for the toe curves
nsections=length(section_levels)-1;
% nsections2=length(section_levels2)-1;

%Preallocation
node.level                      = zeros(1000,1);
node.level(1)                   = pile.head;
node.sigma_v_eff                = zeros(1000,1);
node.sigma_v_eff_noscour        = zeros(1000,1);
element.model_py                = cell(1000,1);
element.model_axial             = cell(1000,1);
element.epsilon50               = zeros(1000,2);
element.J                       = zeros(1000,1);
element.gamma_eff               = zeros(1000,1);
element.phi                     = zeros(1000,1);
element.cu                      = zeros(1000,2);
element.Es                      = zeros(1000,2);
element.G0                      = zeros(1000,2);
element.delta_eff               = zeros(1000,1);
element.c_eff                   = zeros(1000,1);
element.K0                      = zeros(1000,1);
element.thickness               = zeros(1000,1);
element.limit_skin              = zeros(1000,1);
element.limit_alpha             = zeros(1000,1);
element.limit_tip               = zeros(1000,1);
element.poisson                 = zeros(1000,1);
element.k_rm                    = zeros(1000,1);    
element.RQD                     = zeros(1000,1);
element.q_ur                    = zeros(1000,2);
element.level                   = zeros(1000,3);
element.shaft_friction_factor   = zeros(1000,1);
element.Nq                      = zeros(1000,1);
element.degradation_tz_t        = zeros(1000,1);
element.degradation_tz_z        = zeros(1000,1);
element.degradation_py_p        = zeros(1000,1);
element.degradation_py_y        = zeros(1000,1);
element.degradation_mt_m        = zeros(1000,1);
element.degradation_mt_t        = zeros(1000,1);

element.cyclic_ult               = zeros(1000,1);
element.zpeak                   = zeros(1000,1);
element.zres                    = zeros(1000,1);
element.zpeak_Qz                = zeros(1000,1);
element.tres                    = zeros(1000,1);
nel_tot                         = 0;
nel                             = zeros(1,nsections);
element.PISA_param              = zeros(1000,22);
element.soil_layer              = zeros(1000,1);
%Generating element and node data
for i=1:nsections %Through all sections
    
    dsec          = section_levels(i)-section_levels(i+1);    % Length of section
    nel(i)        = ceil(dsec*settings.nelem_factor);         % Number of elements in section
    
    %%%%% Node coordinates
    node.level(nel_tot+1:nel_tot+1+nel(i)) = node.level(nel_tot+1) + dsec * linspace(0,-1,nel(i)+1);
    
    %%%%% Fetching material parameters for current section
    soil_index = find(section_levels(i) <= soil.toplevel,1,'last');
    pile_index = find(section_levels(i) <= pile.cross_section.toplevel,1,'last');
    
    if isempty(pile_index)
        error('No pile cross section data is given for the top of the pile.')
    end
    
    if isempty(soil_index) % if no soil is specified, above seabed
        cu                                              = 0;   
        delta_cu                                        = 0; 
        dist_cu                                         = 0;
        G0                                              = 0;   
        delta_G0                                        = 0;
        epsilon50                                       = 0;
        delta_epsilon50                                 = 0;
		Es 												= 0;
		delta_Es 										= 0;
        q_ur                                            = 0;
        delta_q_ur                                      = 0;
        phi                                             = 0;
        delta_phi                                       = 0;
        element.model_py(nel_tot+1:nel_tot+nel(i))      = {'Zero soil'};
        element.model_axial(nel_tot+1:nel_tot+nel(i))   = {'Zero soil'};
        element.model_tz(nel_tot+1:nel_tot+nel(i))      = {'Zero soil'};
        element.gamma_eff(nel_tot+1:nel_tot+nel(i))     = 0;
        element.epsilon50(nel_tot+1:nel_tot+nel(i))     = 0;
        element.J(nel_tot+1:nel_tot+nel(i))             = 0;
%         element.phi(nel_tot+1:nel_tot+nel(i))           = 0;
%         element.delta_phi(nel_tot+1:nel_tot+nel(i))     = 0;
        element.delta_eff(nel_tot+1:nel_tot+nel(i))     = 0;
        element.c_eff(nel_tot+1:nel_tot+nel(i))         = 0;
        element.K0(nel_tot+1:nel_tot+nel(i))            = 0;
        element.limit_skin(nel_tot+1:nel_tot+nel(i))    = 0;
        element.limit_alpha(nel_tot+1:nel_tot+nel(i))   = 0;
        element.limit_tip(nel_tot+1:nel_tot+nel(i))     = 0;
        element.poisson(nel_tot+1:nel_tot+nel(i))       = 0;
        element.k_rm(nel_tot+1:nel_tot+nel(i))          = 0;
        element.RQD(nel_tot+1:nel_tot+nel(i))           = 0;
        element.Nq(nel_tot+1:nel_tot+nel(i))            = 0;
        element.Nkt(nel_tot+1:nel_tot+nel(i))           = 0;
        element.degradation_tz_t(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.degradation_tz_z(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.degradation_py_p(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.degradation_py_y(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.degradation_mt_m(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.degradation_mt_t(nel_tot+1:nel_tot+nel(i))= ones(nel(i),1);
        element.zpeak(nel_tot+1:nel_tot+nel(i))         = 0; % FKMV
        element.zres(nel_tot+1:nel_tot+nel(i))          = 0; % FKMV
        element.zpeak_Qz(nel_tot+1:nel_tot+nel(i))      = 0; % FKMV
        element.tres(nel_tot+1:nel_tot+nel(i))          = 0; % FKMV
        element.soil_layer(nel_tot+1:nel_tot+nel(i))    = 0; % FKMV
    else % below seabed (scour hole is taken into account below where zero soil is applied
        cu                                              = soil.cu(soil_index);
        delta_cu                                        = soil.delta_cu(soil_index);
        dist_cu                                         = soil.toplevel(soil_index)-node.level(nel_tot+1); % The distance from top of the section to top of the soil layer (where cu/G0 was defined)
        G0                                              = soil.G0(soil_index);
        delta_G0                                        = soil.delta_G0(soil_index);
        epsilon50                                       = soil.epsilon50(soil_index);
        delta_epsilon50                                 = soil.delta_epsilon50(soil_index);
		Es												= soil.Es(soil_index);
		delta_Es 										= soil.delta_Es(soil_index);
        q_ur                                            = soil.q_ur(soil_index);
        delta_q_ur                                      = soil.delta_q_ur(soil_index);
        phi                                             = soil.phi(soil_index);
        delta_phi                                       = soil.delta_phi(soil_index);
        delta_eff                                       = soil.delta_eff(soil_index);
        element.model_py(nel_tot+1:nel_tot+nel(i))      = soil.model_py(soil_index);
        element.model_axial(nel_tot+1:nel_tot+nel(i))   = soil.model_axial(soil_index);
        element.gamma_eff(nel_tot+1:nel_tot+nel(i))     = soil.gamma_eff(soil_index);
        element.J(nel_tot+1:nel_tot+nel(i))             = soil.J(soil_index);
%         element.phi(nel_tot+1:nel_tot+nel(i))           = soil.phi(soil_index);
%         element.delta_phi(nel_tot+1:nel_tot+nel(i))     = soil.delta_phi(soil_index);
%         element.delta_eff(nel_tot+1:nel_tot+nel(i))     = soil.delta_eff(soil_index);
        element.c_eff(nel_tot+1:nel_tot+nel(i))         = soil.c_eff(soil_index);
        element.K0(nel_tot+1:nel_tot+nel(i))            = soil.K0(soil_index);
        element.limit_skin(nel_tot+1:nel_tot+nel(i))    = soil.limit_skin(soil_index);
        element.limit_alpha(nel_tot+1:nel_tot+nel(i))   = soil.limit_alpha(soil_index);
        element.limit_tip(nel_tot+1:nel_tot+nel(i))     = soil.limit_tip(soil_index);
        element.poisson(nel_tot+1:nel_tot+nel(i))       = soil.poisson(soil_index);
        element.k_rm(nel_tot+1:nel_tot+nel(i))          = soil.k_rm(soil_index);
        element.RQD(nel_tot+1:nel_tot+nel(i))           = soil.RQD(soil_index);
        element.level(nel_tot+1:nel_tot+nel(i),3)       = soil.toplevel(soil_index);
        element.Nq(nel_tot+1:nel_tot+nel(i))            = soil.Nq(soil_index);
        element.degradation_tz_t(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_tz_t(soil_index);
        element.degradation_tz_z(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_tz_z(soil_index);
        element.degradation_py_p(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_py_p(soil_index);
        element.degradation_py_y(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_py_y(soil_index);
        element.degradation_mt_m(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_mt_m(soil_index);
        element.degradation_mt_t(nel_tot+1:nel_tot+nel(i))= soil.degradation.value_mt_t(soil_index);
        element.cyclic_ult(nel_tot+1:nel_tot+nel(i))=soil.degradation.Cyclic_Ult(soil_index);
        
        element.zpeak(nel_tot+1:nel_tot+nel(i))         = soil.zpeak; 
        element.zres(nel_tot+1:nel_tot+nel(i))          = soil.zres; 
        element.zpeak_Qz(nel_tot+1:nel_tot+nel(i))      = soil.zpeak_Qz; 
        element.tres(nel_tot+1:nel_tot+nel(i))          = soil.tres; 
        element.soil_layer(nel_tot+1:nel_tot+nel(i))    = soil.layer(soil_index);
        
    end
    
    element.thickness(nel_tot+1:nel_tot+nel(i),1) = pile.cross_section.thickness(pile_index);
    
    %%%%% Correcting for local scour and consequently ORD
    if node.level(nel_tot+1)>scour_level %if the section is above scour level
        gamma_eff_corr = 0;
        element.gamma_eff(nel_tot+1:nel_tot+nel(i)) = gamma_eff_corr; % elements above scour need to have 0 weight in order to calculate plug weight correctly (this is done based on elements)
        element.model_py(nel_tot+1:nel_tot+nel(i)) = {'Zero soil'};
        element.model_axial(nel_tot+1:nel_tot+nel(i)) = {'Zero soil'};
        element.model_tz(nel_tot+1:nel_tot+nel(i)) = {'Zero soil'};
    elseif node.level(nel_tot+1)>ORD_level %if the section is above ORD
        gamma_eff_corr = (element.gamma_eff(nel_tot+1)*(ORD_level-node.level(nel_tot+1))-(node.sigma_v_eff_noscour(nel_tot+1)-node.sigma_v_eff(nel_tot+1)))/(ORD_level-node.level(nel_tot+1));
%         element.gamma_eff(nel_tot+1:nel_tot+nel(i)) = gamma_eff_corr;
    else % if below ORD
        gamma_eff_corr = element.gamma_eff(nel_tot+1);
    end
    
    %%%%% Effective vertical stresses (not taking ORD and local scour into account)
    node.sigma_v_eff_noscour(nel_tot+1:nel_tot+1+nel(i)) = node.sigma_v_eff_noscour(nel_tot+1) + element.gamma_eff(nel_tot+1) * dsec * linspace(0,1,nel(i)+1);
    %%%%% Effective vertical stresses (taking ORD and local scour into account)
    node.sigma_v_eff(nel_tot+1:nel_tot+1+nel(i)) = node.sigma_v_eff(nel_tot+1) + gamma_eff_corr * dsec * linspace(0,1,nel(i)+1);
    
    %%%%% cu in top and bottom of elements
    element.cu(nel_tot+1:nel_tot+nel(i),1) = cu + delta_cu * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.cu(nel_tot+1:nel_tot+nel(i),2) = cu + delta_cu * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);

    %%%%% phi in top and bottom of elements
    element.phi(nel_tot+1:nel_tot+nel(i),1) = phi + delta_phi * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.phi(nel_tot+1:nel_tot+nel(i),2) = phi + delta_phi * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);

    %%%%% delta_eff in top and bottom of elements
    element.delta_eff(nel_tot+1:nel_tot+nel(i),1) = delta_eff + delta_phi * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.delta_eff(nel_tot+1:nel_tot+nel(i),2) = delta_eff + delta_phi * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);
    
    %%%%% G0 in top and bottom of elements
    element.G0(nel_tot+1:nel_tot+nel(i),1) = G0 + delta_G0 * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.G0(nel_tot+1:nel_tot+nel(i),2) = G0 + delta_G0 * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);
    
    %%%%% epsilon50 in top and bottom of elements
    element.epsilon50(nel_tot+1:nel_tot+nel(i),1) = epsilon50 + delta_epsilon50 * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.epsilon50(nel_tot+1:nel_tot+nel(i),2) = epsilon50 + delta_epsilon50 * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);
	
	%%%%% Es in top and bottom of elements
    element.Es(nel_tot+1:nel_tot+nel(i),1) = Es + delta_Es * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.Es(nel_tot+1:nel_tot+nel(i),2) = Es + delta_Es * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);
    
    %%%%% q_ur in top and bottom of elements
    element.q_ur(nel_tot+1:nel_tot+nel(i),1) = q_ur + delta_q_ur * (dsec * linspace(0,(nel(i)-1)/nel(i),nel(i)) + dist_cu);
    element.q_ur(nel_tot+1:nel_tot+nel(i),2) = q_ur + delta_q_ur * (dsec * linspace(1/nel(i),1,nel(i)) + dist_cu);
    

    
    nel_tot       = nel_tot+nel(i);      % Total number of elements;
    
end
    
%Cutting away the preallocated, empty part of vectors;
node.level=node.level(1:nel_tot+1);
node.sigma_v_eff=node.sigma_v_eff(1:nel_tot+1);
node.sigma_v_eff_noscour=node.sigma_v_eff_noscour(1:nel_tot+1);
for k = 1:length(node.level)-1
    element.height(k,1) = abs(node.level(k+1)-node.level(k)); 
end

element.model_py                = element.model_py(1:nel_tot);
element.model_axial             = element.model_axial(1:nel_tot);
element.gamma_eff               = element.gamma_eff(1:nel_tot);
element.epsilon50               = [element.epsilon50(1:nel_tot,1) element.epsilon50(1:nel_tot,2)];
element.J                       = element.J(1:nel_tot);
% element.phi                     = element.phi(1:nel_tot);
element.phi                     = [element.phi(1:nel_tot,1) element.phi(1:nel_tot,2)];
element.cu                      = [element.cu(1:nel_tot,1) element.cu(1:nel_tot,2)];
element.G0                      = [element.G0(1:nel_tot,1) element.G0(1:nel_tot,2)];
element.Es                      = [element.Es(1:nel_tot,1) element.Es(1:nel_tot,2)];
element.q_ur                    = [element.q_ur(1:nel_tot,1) element.q_ur(1:nel_tot,2)];
element.delta_eff               = [element.delta_eff(1:nel_tot,1) element.delta_eff(1:nel_tot,2)];
% element.delta_eff               = element.delta_eff(1:nel_tot);
element.c_eff                   = element.c_eff(1:nel_tot);
element.K0                      = element.K0(1:nel_tot);
element.thickness               = element.thickness(1:nel_tot);
element.limit_skin              = element.limit_skin(1:nel_tot);
element.limit_alpha             = element.limit_alpha(1:nel_tot);
element.limit_tip               = element.limit_tip(1:nel_tot);
element.sigma_v_eff             = [node.sigma_v_eff(1:nel_tot) node.sigma_v_eff(2:nel_tot+1)];
element.level                   = [node.level(1:nel_tot) node.level(2:nel_tot+1) element.level(1:nel_tot,3)];
element.poisson                 = element.poisson(1:nel_tot);
element.k_rm                    = element.k_rm(1:nel_tot);
element.RQD                     = element.RQD(1:nel_tot);
element.shaft_friction_factor   = element.shaft_friction_factor(1:nel_tot);
element.Nq                      = element.Nq(1:nel_tot);
element.degradation_tz_t        = element.degradation_tz_t(1:nel_tot);
element.degradation_tz_z        = element.degradation_tz_z(1:nel_tot);
element.degradation_py_p        = element.degradation_py_p(1:nel_tot);
element.degradation_py_y        = element.degradation_py_y(1:nel_tot);
element.degradation_mt_m        = element.degradation_mt_m(1:nel_tot);
element.degradation_mt_t        = element.degradation_mt_t(1:nel_tot);
element.cyclic_ult =element.cyclic_ult(1:nel_tot);
element.zpeak                   = element.zpeak(1:nel_tot);
element.zres                    = element.zres(1:nel_tot);
element.zpeak_Qz                = element.zpeak_Qz(1:nel_tot);
element.tres                    = element.tres(1:nel_tot);
element.soil_layer              = element.soil_layer(1:nel_tot);

% %adding an extra element below pile tip but with the same properties as the
%last pile element
extra_el_height=element.height(end)/5; %extra element height
node.level=[node.level(1:nel_tot+1); node.level(end)-extra_el_height];
node.sigma_v_eff=[node.sigma_v_eff(1:nel_tot+1); node.sigma_v_eff(end)];
node.sigma_v_eff_noscour=[node.sigma_v_eff_noscour(1:nel_tot+1); node.sigma_v_eff_noscour(end)];
element.model_py                = [element.model_py(1:nel_tot); element.model_py(end)];
element.model_axial             = [element.model_axial(1:nel_tot); element.model_axial(end)];
element.gamma_eff               = [element.gamma_eff(1:nel_tot); element.gamma_eff(end)];
element.epsilon50               = [element.epsilon50(1:nel_tot,1)  element.epsilon50(1:nel_tot,2); element.epsilon50(end,1) element.epsilon50(end,2)];
element.J                       = [element.J(1:nel_tot); element.J(end)];
% element.phi                     = [element.phi(1:nel_tot); element.phi(end)];
element.phi                     = [element.phi(1:nel_tot,1) element.phi(1:nel_tot,2); element.phi(end,1) element.phi(end,2)];
element.cu                      = [element.cu(1:nel_tot,1) element.cu(1:nel_tot,2); element.cu(end,1) element.cu(end,2)];
element.G0                      = [element.G0(1:nel_tot,1) element.G0(1:nel_tot,2);element.G0(end,1) element.G0(end,2)];
element.Es                      = [element.Es(1:nel_tot,1) element.Es(1:nel_tot,2); element.Es(end,1) element.Es(end,2)];
element.q_ur                    = [element.q_ur(1:nel_tot,1) element.q_ur(1:nel_tot,2); element.q_ur(end,1) element.q_ur(end,2)];
element.delta_eff               = [element.delta_eff(1:nel_tot,1) element.delta_eff(1:nel_tot,2); element.delta_eff(end,1) element.delta_eff(end,2)];
% element.delta_eff               = [element.delta_eff(1:nel_tot); element.delta_eff(end)];
element.c_eff                   = [element.c_eff(1:nel_tot); element.c_eff(end)];
element.K0                      = [element.K0(1:nel_tot); element.K0(end)];
element.thickness               = [element.thickness(1:nel_tot); element.thickness(end)];
element.limit_skin              = [element.limit_skin(1:nel_tot); element.limit_skin(end)];
element.limit_alpha             = [element.limit_alpha(1:nel_tot); element.limit_alpha(end)];
element.limit_tip               = [element.limit_tip(1:nel_tot); element.limit_tip(end)];
element.sigma_v_eff             = [node.sigma_v_eff(1:nel_tot) node.sigma_v_eff(2:nel_tot+1); node.sigma_v_eff(end) node.sigma_v_eff(end)];
element.level                   = [node.level(1:nel_tot)  node.level(2:nel_tot+1)  element.level(1:nel_tot,3); node.level(end-1)  node.level(end) element.level(end,3)];%SPSO changed
element.poisson                 = [element.poisson(1:nel_tot); element.poisson(end)];
element.k_rm                    = [element.k_rm(1:nel_tot); element.k_rm(end)];
element.RQD                     = [element.RQD(1:nel_tot); element.RQD(end)];
element.shaft_friction_factor   = [element.shaft_friction_factor(1:nel_tot); element.shaft_friction_factor(end)];
element.Nq                      = [element.Nq(1:nel_tot); element.Nq(end)];
element.degradation_tz_t        = [element.degradation_tz_t(1:nel_tot); element.degradation_tz_t(end)];
element.degradation_tz_z        = [element.degradation_tz_z(1:nel_tot); element.degradation_tz_z(end)];
element.degradation_py_p        = [element.degradation_py_p(1:nel_tot); element.degradation_py_p(end)];
element.degradation_py_y        = [element.degradation_py_y(1:nel_tot); element.degradation_py_y(end)];
element.degradation_mt_m        = [element.degradation_mt_m(1:nel_tot); element.degradation_mt_m(end)];
element.degradation_mt_t        = [element.degradation_mt_t(1:nel_tot); element.degradation_mt_t(end)];

element.cyclic_ult        = [element.cyclic_ult(1:nel_tot); element.cyclic_ult(end)];


element.zpeak                   = [element.zpeak(1:nel_tot); element.zpeak(end)];
element.zres                    = [element.zres(1:nel_tot); element.zres(end)];
element.zpeak_Qz                = [element.zpeak_Qz(1:nel_tot); element.zpeak_Qz(end)];
element.tres                    = [element.tres(1:nel_tot); element.tres(end)];
element.soil_layer              = [element.soil_layer(1:nel_tot); element.soil_layer(end)];
%Calculating ep, pile element parameters
% ep [E A I D_eq Gm ksf]
E = ones(size(element.thickness))*pile.E;
A = ((pile.diameter-0*element.thickness).^2-(pile.diameter-2*element.thickness).^2)*pi/4;
I = ((pile.diameter-0*element.thickness).^4-(pile.diameter-2*element.thickness).^4)*pi/64;
D_eq = sqrt(pile.diameter^2-(pile.diameter-2*element.thickness).^2);
Gm = ones(size(element.thickness))*pile.G;
ksf = ones(size(element.thickness))*pile.ksf;
element.ep = [E A I D_eq Gm ksf];
element.nelem = length(element.model_py); %number of elements in pile

if plots.res_vs_pilelength         == 1
    Min=(abs(abs(element.level(:,1))-(pile.L-pile.extra_L)));
    pile.axial_elem=find(Min==min(Min))-1;
    if size(pile.axial_elem,1) > 1
        pile.axial_elem = pile.axial_elem(1);
    end
end


%% PISA parameters
element.PISA_switch = settings.PISA_database;
if settings.PISA_database 
for iii = 1:element.nelem % FKMV
    [element] = PISA_formulations_DB(pile,soil,element,scour,settings,data,iii,Database); % FKMV - New PISA organisation giving PISA parameters as a vector  
    [element] = PISA_formulations_calculation(pile,soil,element,scour,settings,data,iii); % calculates the final parameters
end
element.PISA_param                      = element.PISA_param(1:nel_tot+1,:); % cutting away non needed pre allocated end
if settings.axial_loading == 0
    [element] = check_PISA_param(element); % Checks for unacceptable PISA parameters
end
end

element.height = [element.height; element.height(end)/5];