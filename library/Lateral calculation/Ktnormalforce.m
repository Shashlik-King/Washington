function Ktnormal = Ktnormalforce(ex,ey)
% 
% 
%---------------------------------------------------------------------
% PURPOSE
% Compute the global element tangent stiffness matrix due to the contribu-
% tion from the normal force to the deflections. It is assumed that normal
% force is constant along the considered pile segment (element).
% 
% INPUT:  ex      = [x1 x2]
%         ey      = [y1 y2]       element node coordinates
%
% OUTPUT: Ktnormal: Element tangent stiffness matrix (6 x 6)
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    26.10.2007
% 
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1) ];
  L     =   sqrt(b'*b);  n=b/L;

  
  Ktnorm = [ 0   0          0       0   0         0;
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
  
  Ktnormal = G'*Ktnorm*G;
