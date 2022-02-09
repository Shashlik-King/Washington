function [element Ndof Coord Mmax Es u ustep Ex Ey output] = winkler_axial(element,node,pile,loads,reduction,plug_unplug,G,settings,output,plots,SF,A,fs,qp,ii)
% -------------------------------------------------------------------------
% This program calculates deformations, moments etc. of a monopile based
% on a Winkler approach where the soil resistance is represented by t-z and
% q-z curves.
% CALFEM has been applied.
%
% VARYING MATERIAL PROPERTIES OF PILE AND SOIL WITH DEPTH
%
% -------------------------------------------------------------------------
% CODE            : EBGX
% APPROVED        : 
%
% Log
%
% 24.10.2016    EBGX     Programming based on Winkler for horizontal loads.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
% kstop;                   % Secant stiffness of soil spring in top and bot-
% ksbot;                   % tom of pile segment [KN/m/m].
% qespring;                % Internal forces in beam due to spil springs.
% qebeam;                  % Internal forces in beam.
% qenormal;                % Internal forces in beam due to normal force.
% q;                       % Global internal force vector.
% kttop;                   % Tangent stiffness of soil spring in top and bot-
% ktbot;                   % tom of pile segment [KN/m/m].
% Ktspring1;               % Global tangent stifness matrix based on soil.
% Ktbeam;                  % Global tangent stifness matrix based on the beam.
% Ktnormal;                % Global tangent stifness - contribution from 
%                          % normal force.
% Kt;                      % Global tangent stiffnes matrix for the beam 
%                          % and soil.
% r;                       % Residual.
% Q;                       % Reaction forces.
% du;                      % Increment in displacement.
% u;                       % Displacement vector

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
if strcmp(settings.analysis_loading(ii),'Tens') == 1
    P = [1 0 
         2 -loads.Vt+sum(G.pile)*SF.weight_advan 
         3 0]; %Loads given at pile head in global coordinates (dof=1, 2 and 3), horizontal force neglected in winkler-analysis
elseif strcmp(settings.analysis_loading(ii),'Comp') == 1
    P = [1 0 
         2 loads.Vc+sum(G.pile)*SF.weight_advan 
         3 0];
end 

% nelem = element.nelem; %number of elements in pile
nelem = element.nelem-1; %number of elements in pile FKMV

% BC = [(3*nelem+2) 0]; % Free pile head 
%BC = [ 3 0 ; (3*nelem+2) 0 ]; % Fixed pile head    
BC  = [ 1 0 ; 3 0 ]; 
%BC = [ (3*nelem+2) 0];
% Global prescribed nodal degrees of freedom matrix, C = [dof value]
% -------------------------------------------------------------------------
% Topology for beam
% -------------------------------------------------------------------------
% Coord = [zeros(nelem+1,1) node.level];
% Dof = 1:3:(3*nelem+1);
% Dof = [Dof'  Dof'+1  Dof'+2];
% Edof = [(1:nelem)' Dof(1:end-1,:) Dof(2:end,:)];

Coord = [zeros(nelem+1,1) node.level(1:nelem+1)]; %FKMV:
Dof = 1:3:(3*nelem+1);
Dof = [Dof'  Dof'+1  Dof'+2];
Edof = [(1:nelem)' Dof(1:end-1,:) Dof(2:end,:)];

%Edof;                   % The element topology matrix
                         % Edof = [el dof1 .... dof6 ]
                         % el = element number
                         % dof = global degree of freedom for el

%Coord;                  % Global node coordinate matrix
                         % Coord = [x1 y1; x2 y2 etc.], 1,2 = Nodes

%Dof;                    % The global topology matrix
                         % DOF = [k1 l1 m1; k2 l2 m2]
                         % k,l,m = dof for given node, 1,2 = Node

nen = 2;                 % Number of element nodes in each element.


[Ex Ey] = coordxtr(Edof,Coord,Dof,nen);
                         % Element nodal coordinates
                         % Ex = [x11 x21; x12 x22; ... x1nelem x2nelem]
                         % x-cordinates, first number 1,2 = element node
                         % second number = element
                         
% -------------------------------------------------------------------------
% Stiffness matrix, load vector and solving equations
% -------------------------------------------------------------------------
% -------- Initializing stiffness matrix and load vector ------------------

nnodes  = length(Dof(:,1));     % Number of nodes
dofnode = length(Dof(1,:));     % Degrees of freedom per node
Ndof    = nnodes*dofnode;       % Total number of degrees of freedom

F       = zeros(Ndof,1);        % Initializing global load vector
f       = zeros(Ndof,1);        % Initializing load vector for iteration
u       = zeros(Ndof,1);        % Initializing global displacement vector
Q       = zeros(Ndof,1);        % Initializing global reaction vector

for j   = 1:length(P(:,1))
    F(P(j,1)) = P(j,2); %Inserting the prescribed loads in the gload load vector
end

% -------- Global element stiffness matrix and assembling -----------------
[nelem,n]       = size(Edof);
t               = Edof(:,2:n);
if plots.load_deflection_axial == 0
	df              = F/(settings.n_max);         % Load step size
elseif plots.load_deflection_axial == 1
	df              = F/settings.n_max*settings.max_load_ratio;
	disp({'Attempts to apply' settings.max_load_ratio 'times input load'});
end

ustep = NaN(Ndof,settings.n_max); %Preallocation
fstep = NaN(Ndof,settings.n_max); 
stop = 0; % Preallocation
n_possible = settings.n_max; % Preallocation
output.n_possible = n_possible; 
no_convergence = 0; 

for n = 1:settings.n_max                    % Looping over load increments

        f = f + df;
   
    j = 1; 
    while j < settings.j_max+1                 % N-R Iteration procedure

        q  = zeros(Ndof,1);        % Initializing internal force vector
        q_spring = zeros(Ndof,1); 
        Kt = zeros(Ndof);          % Initializing global stiffness matrix
        Kts = zeros(Ndof);

        %Preallocation
        ks_tot=NaN(nelem,2);
       
        
        for i = 1:nelem % FKMV
          
            % Internal forces calculated and assembling
            % Description of variables after the N-R procedure 

            z_topbottom = [u(3*i-1) u(3*i+2)]; % vertical disp. at the top and bottom of the pile segment [m]
            
            [kstop ksbot]       = secspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii);
            
            qespring            = qspring_axial(Ex(i,:),Ey(i,:),kstop,ksbot,u,t,i);
			if settings.beam_theory
                qebeam              = qbeam2t(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
            else
                qebeam              = qbeam2e(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
			end
			
            qelem               = qebeam + qespring;
            q_spring(t(i,:))    = q_spring(t(i,:)) + qespring; 
            q(t(i,:))           = q(t(i,:)) + qelem;
                
            ks_tot(i,1)=kstop;
            ks_tot(i,2)=ksbot;
            
            if i == nelem && strcmp(settings.analysis_loading(ii),'Comp') == 1
                zbot = u(3*i+2);         % Horizontal disp. at the bottom of the last pile element [m]
                Q_end    = Qz_spring(element,pile,plug_unplug,zbot,qp,i);
                q(i*3+2) = q(i*3+2)+Q_end;
            end
            % Global tangent stiffness matrix, Kt
            % Description of variables after the N-R procedure
                   
            z_topbottom = [abs(u(3*i-1)) abs(u(3*i+2))];
            
            [kttop ktbot]       = tanspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii);
            
            Ktspring1           = Ktspring_axial(Ex(i,:),Ey(i,:),kttop,ktbot);
			if settings.beam_theory
                Ktbeam              = beam2t(Ex(i,:),Ey(i,:),element.ep(i,:));
            else
                Ktbeam              = beam2e(Ex(i,:),Ey(i,:),element.ep(i,:));
			end
			
            Ktelem              = Ktbeam + Ktspring1;
            Kts(t(i,:),t(i,:))  = Kts(t(i,:),t(i,:)) + Ktspring1; 
            Kt(t(i,:),t(i,:))   = Kt(t(i,:),t(i,:)) + Ktelem;
            
        end
        % Residual force vector where the reaction forces Q are taken into
        % consideration. This is done since the internal forces q are
        % calculated based on boundary conditions and f are not. Q balances
        % thereby the internal forces due to the boundary conditions.
        
        r = f - q - Q;

        % Solving the equations (determination of disp. increments) by
        % taking boundary conditions into account. Determination of new
        % displacement vector (u + du) and reaction force vector.

        [du,Q]                  = solveq(Kt,r,BC);

        u                       = u + du;
        
        if isnan(u(2))
            fprintf(['Displacement vector u contain entrances with NaN in load step ' num2str(n) ' - increment ' num2str(j)])
            j = settings.j_max; 
            u = zeros(Ndof,1); 
            Q = zeros(Ndof,1); 
        end 

        % Termination criterion for the iteration procedure
        if rednorm(r,Ndof,BC) < settings.TOL*rednorm(df,Ndof,BC)
        %if norm(r)<settings.TOL*norm(df)
             fprintf(' Convergence in step %2i - increment %2i\n',n,j)
             if no_convergence == 1 
                stop = 1; 
             end    
             break
		elseif j == settings.j_max && plots.load_deflection_axial == 1
            stop = 1;
            break
        elseif j == settings.j_max
            error([' No convergence in step ' num2str(n) ' - increment ' num2str(j) '\n'])
%             df = 0.5*df; 
%             f = f-df; 
%             no_convergence = 1; 
%             j=0; 
        end            
        j = j+1; 
    end
    
    ustep(:,n) = u;
    fstep(:,n) = f; % Save load for plotting
	
    if stop == 1
        if no_convergence == 1
            n_possible = n; 
        else 
            n_possible = n-1;
        end 
        output.n_possible = n_possible;
        break
    end

end
  
output.deflections = ustep;
output.ks=ks_tot;
output.loads = fstep; % Move fstep to output

if round(output.loads(2,output.n_possible)) < round(F(2)) 
    fprintf(2,' Final load step is smaller than defined load - disregard plot \n') 
end 


% -------------------------------------------------------------------------
% Determination of section forces and controls
% -------------------------------------------------------------------------

% -------- Deformed and undeformed pile -----------------------------------
Ed=cell(settings.n_max,1); %preallocation
for m =1:settings.n_max
Ed{m}  = extract(Edof,ustep(:,m));          % Element displacements from 
                                            % global solution vector, u.
                                            % Ed = [u1', u2' .... ], e.g.
                                            % every row contains the disp.
                                            % from each element 
end

eq     = [0 0];                             % Distributed load along element

% Computation of section forces and displacements in local directions for
% all load steps
for l = 1:min([settings.n_max n_possible])
    for k = 1:nelem      
		if settings.beam_theory
                [Es{l}(k,:,:) Edi{l}(k,:,:) Eci{l}(k,:)] = beam2tsec(Ex(k,:),Ey(k,:),element.ep(k,:),Ed{l}(k,:),eq,2);  
        else
                [Es{l}(k,:,:) Edi{l}(k,:,:) Eci{l}(k,:)] = beam2sec(Ex(k,:),Ey(k,:),element.ep(k,:),Ed{l}(k,:),eq,2);        
        end
	end
end

% -------- Control of yielding of pile ------------------------------------
for m = 1:nelem
    I               = element.ep(m,3);
    A               = element.ep(m,2);
    W               = I/(pile.diameter);        % Moment of resistance [m^3]
    M(m,:)          = abs(Es{min([settings.n_max n_possible])}(m,:,3));    % Moment in pile segment
    sigmasteel(m,:) = abs(Es{min([settings.n_max n_possible])}(m,:,1))/A;  % Normal stress
end

sigmasteelmax       = max(max(sigmasteel));     % Max. stress in steel

Mmax                = max(max(abs(Es{min([settings.n_max n_possible])}(:,:,3)))); % Max. moment

Efficiency          = sigmasteelmax/pile.sigma_y*100;     % Rate of utilization [%]