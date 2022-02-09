function qebeam = qbeam2t(ex,ey,ep,u,t,i)
% 
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the internal forces in a beam element.
% 
% INPUT:  ex      : Local nodal ccordinates of the beam elements
%         ey      : ex = [x1 x2], ey = [y1 y2]
%         ep      : Global material property matrix ep = [E A I]
%         u       : Global displacment vector
%         t       : Each row in t describes the global dof associated with
%                   Edof
%         i       : Counter referring to element number
%
% OUTPUT: qebeam  : Internal forces in a beam element (6 x 1)
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA  

% LAST MODIFIED   : AHA    13.07.2007   Programming
%                   AHA/LA 14.08.2007   Correcting - Transformation matrix
%                                       G and G' is multiplied G instead of
%                                       only G
%                   EVVA   15.02.2016   Changed Kle to Timoshenko beam  
%--------------------------------------------------------------------------


% -------- Initializing pile and soil parameters --------------------------

b=[ ex(2)-ex(1); ey(2)-ey(1) ];
L=sqrt(b'*b);  n=b/L;

E=ep(1);  A=ep(2);  I=ep(3); Gm=ep(5); ksf=ep(6); 
 
m=(12/L^2)*(E*I/(Gm*A*ksf));

ue = u(t(i,:)');


% -------- Internal forces in beam element --------------------------------

Kle  = E/(1+m)*[A*(1+m)/L    0         0       -A*(1+m)/L   0            0;
                  0      12*I/L^3  6*I/L^2        0      -12*I/L^3   6*I/L^2;
                  0      6*I/L^2  4*I*(1+m/4)/L   0      -6*I/L^2  2*I*(1-m/2)/L;
              -A*(1+m)/L    0         0        A*(1+m)/L    0           0;
                  0     -12*I/L^3  -6*I/L^2       0      12*I/L^3   -6*I/L^2;
                  0      6*I/L^2  2*I*(1-m/2)/L   0      -6*I/L^2   4*I*(1+m/4)/L];

     
G    = [ n(1) n(2)  0    0    0   0;
        -n(2) n(1)  0    0    0   0;
          0    0    1    0    0   0;
          0    0    0   n(1) n(2) 0;
          0    0    0  -n(2) n(1) 0;
          0    0    0    0    0   1 ];


qebeam = G'*Kle*G*ue;
        
