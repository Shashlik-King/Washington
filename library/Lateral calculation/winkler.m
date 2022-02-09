function [element Ndof Coord Mmax Es u ustep output] = winkler(element,node,pile,loads,settings,output,plots,data)
% -------------------------------------------------------------------------
% This program calculates deformations, moments etc. of a monopile based
% on a Winkler approach where the soil resistance is represented by p-y
% curves.
% CALFEM has been applied.
%
% VARYING MATERIAL PROPERTIES OF PILE AND SOIL WITH DEPTH
% NORMAL FORCE INFLUENCES DEFLECTION AND SECTION FORCES
% EFFECTS OF LAYERED SOIL
%
% -------------------------------------------------------------------------
% Log
%
% 13.08.2007    AHA     Programming.
% 14.08.2007    AHA/LA  Correcting - including reaction forces when
%                       calculating the residual.
%                       Correcting - qbeam2e, qspring and secspringstiff
% 21.08.2007    AHA     Ex and Ey is applied in the calculations instead of
%                       ex and ey so the calculations are independent on
%                       way the elements are numbered.
% 26.10.2007    AHA     Including the effects of a normal force.
% 21.11.2007    AHA     Introduction of varying cross section properties
%                       and soil models.
% 22.05.2008    AHA     Effects of layered soil included.
% 26.05.2008    AHA     Cyclic behaviour associated with soft clay
%                       criterion
% 19.08.2008    AHA     Free standing length and modification of
%                       calculation of shear forces
% 16.01.2009    AHA     Plotting p-y curves
%                       Plotting results corresponding to more than one
%                       calculation step
%                       Saving relevant parameters
%                       Plotting soil capacity on screen
% 2010.11.04    MMOL    Integrating into calculation routine
% 2011.05.05    JALY    Streamlining code
% 2011.07.07    MMOL    Re-programming the calculation of itilization ration
% 2013.01.08    MUOE    Taking utilization ratio and calculation of p-y
%                       curves of out the file.
% 2016.04.13    EVVA    Added Timoshenko beam theory part.
% 2017.09.20    EVVA    Added PISA model
% 2018.04.04	DATY	On/Off switch for Georgiadis approach
% 2019.11.21    ASSV    Averaging of calculated internal base moment forces
%                       to enhance convergence of the solver
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
% kstop;                   % Secant stiffness of soil spring in top and bot-
% ksbot;                   % tom of pile segment [KN/m/m].
% qespring;                % Internal forces in beam due to spil springs.
% qebeam;                  % Internal forces in beam.
% qenormal;                % Internal forces in beam due to normal force.
% qemoment;                % Internal forces in beam due to distrib. moment.
% q;                       % Global internal force vector.
% kttop;                   % Tangent stiffness of soil spring in top and bot-
% ktbot;                   % tom of pile segment [KN/m/m].
% Ktspring1;               % Global tangent stifness matrix based on soil.
% Ktbeam;                  % Global tangent stifness matrix based on the beam.
% Ktnormal;                % Global tangent stifness - contribution from
%                          % normal force.
% Ktmoment1;               % Global tangent stifness - contribution from
%                          % distributed moment.
% Kt;                      % Global tangent stiffnes matrix for the beam
%                          % and soil.
% r;                       % Residual.
% Q;                       % Reaction forces.
% du;                      % Increment in displacement.
% u;                       % Displacement vector

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
P = [1 loads.H
    2 0
    3 loads.M]; %Loads given at pile head in global coordinates (dof=1, 2 and 3), vertical force neglected in winkler-analysis

nelem = element.nelem; %number of elements in pile

BC = [  (3*nelem+2) 0 ]; % Free pile head
% BC = [ 3 0 ; (3*nelem+2) 0 ]; % Fixed pile head
% Global prescribed nodal degrees of freedom matrix, C = [dof value]
% -------------------------------------------------------------------------
% Topology for beam
% -------------------------------------------------------------------------
Coord = [zeros(nelem+1,1) node.level];
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
if plots.load_deflection == 0
    df              = F/settings.n_max;         % Load step size
elseif plots.load_deflection == 1
    df              = F/settings.n_max*settings.max_load_ratio;
    disp({'Attempts to apply' settings.max_load_ratio 'times input load'});
end

ustep = NaN(Ndof,settings.n_max); %Preallocation
stop = 0; % Preallocation
n_possible = settings.n_max; % Preallocation
for n = 1:settings.n_max                    % Looping over load increments
    
    f = f + df;
    qm_stored = zeros(Ndof,settings.j_max); % Initialize array to store
    % internal base moment forces for every iterations
    for j = 1:settings.j_max                % N-R Iteration procedure
        
        q  = zeros(Ndof,1);        % Initializing internal force vector
        Kt = zeros(Ndof);          % Initializing global stiffness matrix
        
        % Including effects due to layered soil - The method of
        % Georgiadis (1983) is employed.
        % Description of variables after the N-R procedure
        [element] = layer(element,node,pile,settings,loads);
        
        %Preallocation
        ks_tot=NaN(nelem,2);
        ksm_tot=NaN(nelem,2);
        for i = 1:nelem
            if i == 36
               ptest = 1;  
            end
            % Internal forces calculated and assembling
            % Description of variables after the N-R procedure
            
            y_topbottom         = [abs(u(3*i-2)) abs(u(3*i+1))]; % Horizontal disp. at the top and bottom of the pile segment [m]
            teta_topbottom      = [abs(u(3*i)) abs(u(3*i+3))]; % Nodal section rotation at the top and bottom of the pile segment [rad]
            if i<=(nelem-1) %only pile nodes
                [kstop ksbot]       = secspringstiff(element,pile,loads,y_topbottom,i); % secant stiffness for distributed lateral load
                if settings.mteta
                    [ksmomtop ksmombot] = secmomstiff(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);
                else
                    ksmomtop = 0;
                    ksmombot = 0;
                end
            elseif settings.toe_shear %only the extra node for base element
                [kstop ksbot]       = secHBstiff(element,pile,y_topbottom,i); % secant stiffness for base shear load
                [ksmomtop ksmombot] = secMBstiff(element,pile,teta_topbottom,i);  % secant stiffness for base moment
            else %only the extra node for base element
                kstop = 0;
                ksbot = 0;
                ksmomtop = 0;
                ksmombot = 0;
            end
            qespring            = qspring(Ex(i,:),Ey(i,:),kstop,ksbot,u,t,i);
            qemoment            = ((ksmomtop+ksmombot)/2)*qmoment(Ex(i,:),Ey(i,:),u,t,i);
            if settings.toe_shear && j>2 && i>(nelem-1)
                qemoment = (qemoment+qm_stored(t(i,:),j-1))/2;%averaging with previous iteration
            end
            if settings.beam_theory
                qebeam              = qbeam2t(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
            else
                qebeam              = qbeam2e(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
            end
            
            qenormal            = P(2,2)*qnormal(Ex(i,:),Ey(i,:),u,t,i);
            if settings.mteta
                qelem               = qebeam + qespring + qenormal + qemoment;
            else
                qelem               = qebeam + qespring + qenormal;
            end
            q(t(i,:))           = q(t(i,:)) + qelem;
            
            ks_tot(i,1)=kstop;
            ks_tot(i,2)=ksbot;
            ksm_tot(i,1)=ksmomtop;
            ksm_tot(i,2)=ksmombot;
            qm_stored(t(i,:),j)=qemoment;
            
            % if i == nelem
            % ybot = u(3*i+1);         % Horizontal disp. at the bottom of the pile segment [m]
            % if settings.toe_shear == 1
            % ts = toe_shear(element,node,pile,ybot,i);
            % elseif settings.toe_shear == 0
            % ts = 0;
            % end
            % q(i*3+1) = q(i*3+1)+ts; % Addition of toe shear force to load vector
            % end
            % Global tangent stiffness matrix, Kt
            % Description of variables after the N-R procedure
            
            if i<=(nelem-1) %only pile nodes
                [kttop ktbot]       = tanspringstiff(element,pile,loads,y_topbottom,i);
                if settings.mteta
                    [ktmomtop ktmombot] = tanmomstiff(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,i);  % tangent stiffness for distributed moment
                else
                    ktmomtop = 0;
                    ktmombot = 0;
                end
                
            elseif settings.toe_shear %only the extra node for base element
                [kttop ktbot]       = tanHBstiff(element,pile,y_topbottom,i); % tangent stiffness for base shear load
                [ktmomtop ktmombot] = tanMBstiff(element,pile,teta_topbottom,i);  % tangent stiffness for base moment
            else
                kttop = 0;
                ktbot = 0;
                ktmomtop = 0;
                ktmombot = 0;
            end
            Ktspring1           = Ktspring(Ex(i,:),Ey(i,:),kttop,ktbot);
            Ktmoment1           = Ktmoment(Ex(i,:),Ey(i,:),ktmomtop,ktmombot);
            if settings.beam_theory
                Ktbeam              = beam2t(Ex(i,:),Ey(i,:),element.ep(i,:));
            else
                Ktbeam              = beam2e(Ex(i,:),Ey(i,:),element.ep(i,:));
            end
            
            Ktnormal            = P(2,2)*Ktnormalforce(Ex(i,:),Ey(i,:));
            if settings.mteta
                Ktelem              = Ktbeam + Ktspring1 + Ktnormal + Ktmoment1;
            else
                Ktelem              = Ktbeam + Ktspring1 + Ktnormal;
            end
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
        
        % Termination criterion for the iteration procedure
        if rednorm(r,Ndof,BC) < settings.TOL*rednorm(df,Ndof,BC)
            fprintf(' Convergence in step %2i - increment %2i\n',n,j)
            break
        elseif j == settings.j_max && plots.load_deflection == 1
            stop = 1;
            break
        elseif j == settings.j_max
            error([' No convergence in step ' num2str(n) ' - increment ' num2str(j)])
        end
    end
    if stop == 1
        n_possible = n-1;
        output.n_possible = n_possible;
        break
    end
    ustep(:,n) = u;
end

output.deflections = ustep;
output.ks=ks_tot;
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
    sigmasteel(m,:) = M(m,:)/W + abs(P(2,2))/A; % Naviers equation
end

sigmasteelmax       = max(max(sigmasteel));     % Max. stress in steel

Mmax                = max(max(abs(Es{min([settings.n_max n_possible])}(:,:,3)))); % Max. moment

Efficiency          = sigmasteelmax/pile.sigma_y*100;     % Rate of utilization [%]
end