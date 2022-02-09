function Ktmoment1 = Ktmoment(ex,ey,ktmomtop,ktmombot)
% 
% 
%---------------------------------------------------------------------
% PURPOSE
% Compute the internal forces in a pile segment due to the contributions 
% from the uniformly distributed moment.
% 
% INPUT:  ex      = [x1 x2]
%         ey      = [y1 y2]       Element node coordinates
%         u       : Global displacement vector
%         t       : Each row in t describes the global dof associated with
%                   Edof
%         i       : Counter referring to element number
%
% OUTPUT: qemoment: Internal force vector due to the contribution from
%                   uniformly distributed moment (6 x 1)
%
% Log:
% EVVA    18.08.2016   Programming
%                   
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1) ];
  L     =   sqrt(b'*b);  n=b/L;
  
  
  
  Ktmom = [ 0   0          0       0   0         0;
             0   6/(5*L)    1/10    0  -6/(5*L)   1/10;
             0   1/10       2*L/15  0  -1/10     -L/30; 
             0   0          0       0   0         0;
             0   -6/(5*L)  -1/10    0   6/(5*L)  -1/10;
             0   1/10      -L/30    0  -1/10      2*L/15 ];

         
            
  G      = [ n(1) n(2)  0    0    0   0;
            -n(2) n(1)  0    0    0   0;
              0    0    1    0    0   0;
              0    0    0   n(1) n(2) 0;
              0    0    0  -n(2) n(1) 0;
              0    0    0    0    0   1 ];
  
  Ktmoment1 = G'*((ktmomtop+ktmombot)/2)*Ktmom*G;
  