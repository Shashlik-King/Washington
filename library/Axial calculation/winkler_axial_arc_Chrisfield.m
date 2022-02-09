function [element Ndof Coord Mmax Es u ustep Ex Ey output] = winkler_axial_arc_Chrisfield(element,node,pile,loads,reduction,plug_unplug,G,settings,output,plots,SF,A,fs,qp,ii)
% -------------------------------------------------------------------------
% This program calculates deformations, moments etc. of a monopile based
% on a Winkler approach where the soil resistance is represented by t-z and
% q-z curves.
% CALFEM has been applied. The iteration scheme applied is arc-length with
% the Chrisfield approach. 
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
warning('Chosen Arc-length solver is not fully QA-ed. Its use is not recomended until completely checked and approved')
% z_peak = 0.01*pile.diameter;
% z_res = 0.02*pile.diameter; % FKMV
%        strcmp(setting.analysis_loading(ii),'Tens')
if strcmp(settings.analysis_loading(ii),'Tens') == 1
    P = [1 0 
         2 -loads.Vt+sum(G.pile)*SF.weight_advan 
         3 0]; %Loads given at pile head in global coordinates (dof=1, 2 and 3), horizontal force neglected in winkler-analysis
elseif strcmp(settings.analysis_loading(ii),'Comp') ==1
    P = [1 0 
         2 loads.Vc+sum(G.pile)*SF.weight_advan 
         3 0];
end 

% nelem = element.nelem; %number of elements in pile
nelem = element.nelem-1; %number of elements in pile

BC  = [ 1 0 ; 3 0 ]; 

% Global prescribed nodal degrees of freedom matrix, C = [dof value]


% -------------------------------------------------------------------------
% Topology for beam
% -------------------------------------------------------------------------
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
r       = zeros(Ndof,1); 
u       = zeros(Ndof,1); 
Delta_lambda_bar = settings.Delta_lambda_bar;        % Initial load increment. Used for first load step and to determine arc length  
Q       = zeros(Ndof,1);        % Initializing global reaction vector

for j   = 1:length(P(:,1))
    F(P(j,1)) = P(j,2); %Inserting the prescribed loads in the gload load vector
end

% -------- Preallocation and assempling -----------------
[nelem,n]       = size(Edof);
t               = Edof(:,2:n);

ustep = NaN(Ndof,settings.n_max); %Preallocation
fstep = NaN(Ndof,settings.n_max); % Preallocation 
stop = 0; % Preallocation
n_possible = settings.n_max; % Preallocation
output.n_possible = n_possible; 
n=0; 
lambda0 = 0; 

% Find z_stop  
for i = 1:nelem 
    z_stop(i) = element.zpeak(i)*element.zres(i); 
end 
if strcmp(settings.analysis_loading(ii),'Tens') == 0
    z_stop = min(z_stop); 
elseif strcmp(settings.analysis_loading(ii),'Tens') == 1
    z_stop = -min(z_stop);
end 
 
while n<settings.n_max , n=n+1;              % Looping over load increments
    
    if n==1
        lambdap = lambda0;%+Delta_lambda_bar; % Initial lambdap
        lambdai = lambda0+Delta_lambda_bar; % Intial load increment 
        up = zeros(Ndof,1); % Initializing global displacement vector to contain displacement from last converged point 
        ui = zeros(Ndof,1); % Initializing global displacement vector 
    end 
    j = 1; 
    while j<settings.j_max+1 % N-R Iteration procedure                 
        % Preallocation
        ks_tot=NaN(nelem,2);
        q  = zeros(Ndof,1);        % Initializing internal force vector
        Kt = zeros(Ndof);          % Initializing global stiffness matrix
        Kts = zeros(Ndof);
        DeltauRi = zeros(Ndof,1); 
        DeltauFi = zeros(Ndof,1);

        
        % Check conditions
        TF = isnan(ui); 
        if any(TF==1) 
            stop = 1; 
            break; 
        end 
        if strcmp(settings.analysis_loading(ii),'Tens') == 0
            if max(ui) > z_stop
                stop = 1; 
                break;
            end 
            if plots.load_deflection_axial == 0
            if max(f) > max(F) 
                stop = 1; 
                break; 
            end 
            end 
        elseif strcmp(settings.analysis_loading(ii),'Tens') == 1
            if min(ui) < z_stop 
                stop = 1; 
                break; 
            end
            if plots.load_deflection_axial == 0
            if min(f) < min(F)
                stop = 1; 
                break; 
            end 
            end 
        end        
             
        [Kt Kts] = Kt_axial(ui,element,pile,reduction,plug_unplug,loads,Ex,Ey,nelem,settings,t,Kts,Kt,A,fs,ii); % Determine tangent stiffness matrix 
        [q,ks_tot] = internal_forces(nelem,ui,element,pile,reduction,plug_unplug,loads,Ex,Ey,t,settings,q,ks_tot,A,fs,qp,ii); % Determine internal forces 
        
        Df = lambdai*F; % Calculate applied load 
        r = Df-q-Q; % Residual load 
        
        [DeltauRi,Q] = solveq(Kt,r,BC); % Displacements and reaction forces from residual force 
        
        DeltauFi = solveq(Kt,F,BC); % Displacements from defined load 
        
        if j==1 
            if n==1
                DeltaS = Delta_lambda_bar*sqrt(DeltauFi'*DeltauFi); % Arc-length 
                Deltaui = Delta_lambda_bar*DeltauFi; % Displacement increment for j = 1 and n = 1
                ui = ui + Deltaui; % Total global displacement vector 
            else 
                AA = (up-up_1)'*DeltauFi;
                Delta_lambdai = sign(AA)*DeltaS/sqrt(DeltauFi'*DeltauFi); % Change in load increment 
                lambdai = lambdap + Delta_lambdai; % Load increment for iteration step j = 1
                Deltaui = DeltauRi+Delta_lambdai*DeltauFi; % Displacement increment for iteration step j=1 
                ui = ui+Deltaui;  % Total global displacement vector 
            end 
        else % Determine load increment and displacement increment for iteration step j>1
            vi = (ui-up)+DeltauRi;
            a = DeltauFi'*DeltauFi; 
            b = DeltauFi'*vi; 
            c = vi'*vi-DeltaS^2; 
            upi = ui-up; 
            if (b^2-a*c) > 0 
                Delta_lambdai=[(-b-sqrt(b^2-a*c))/a,(-b+sqrt(b^2-a*c))/a]; % Change in load increment 
                Deltaui=[DeltauRi+Delta_lambdai(1)*DeltauFi,DeltauRi+Delta_lambdai(2)*DeltauFi]; % Displacement increment calculated with new change in load increment
                ui1=[ui+Deltaui(:,1),ui+Deltaui(:,2)]; % Total displacements 
                upi1=[ui1(:,1)-up,ui1(:,2)-up]; % Difference between the new total displacement and the displacements from last converged point 
                cosines=[upi'*upi1(:,1)/abs(upi'*upi1(:,1)),upi'*upi1(:,2)/abs(upi'*upi1(:,2))]; % Cosines values 
                 % root selection (least posive value of the two cosines)
                [~,IX]=sort(cosines);
                Delta_lambdai=Delta_lambdai(IX(2));
                Deltaui=Deltaui(:,IX(2));
                ui=ui1(:,IX(2));
            else % the roots contain imaginary numbers 
                Delta_lambdai=-b/a;
                Deltaui=DeltauRi+Delta_lambdai*DeltauFi; % eqs.(28)&(29) in [3]
                ui=ui+Deltaui; % eq.(3) in [3]
                stop = 1; 
                break
            end
            lambdai = lambdai+Delta_lambdai; % The change in the load increment is added to the total load increment 
 
   %% Termination criterion for the iteration procedure
            if rednorm(r,Ndof,BC) < settings.TOL*rednorm((lambdai*F - lambdap*F),Ndof,BC)   
                fprintf(' Convergence in step %2i - increment %2i\n',n,j)
                break
            elseif j==settings.j_max
                fprintf(' No convergence in step %2i - increment %2i\n',n,j) 
                stop = 1; 
                break 
            end  
        end 



        j=j+1;

    end 

  f = F*lambdai;  
     ustep(:,n) = ui; % Save displacement for plotting 
     fstep(:,n) = f; % Save load for plotting 
lambdap = lambdai; 
up_1 = up; 
up = ui; 
u = ui; 

	if stop == 1
        n_possible = n-1;
        output.n_possible = n_possible;
        break; 
    end


end 
   

output.loads = fstep; % Move fstep to output 
output.deflections = ustep; % Move ustep to output 
output.ks=ks_tot; % Move ks_tot to output 

if plots.load_deflection_axial == 0 
    if max(output.loads(2,:)) < F(2) 
        fprintf(2,'Final load step is smaller than defined load \n') 
    end 
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
end