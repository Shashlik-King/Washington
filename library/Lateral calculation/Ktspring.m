function Ktspring1 = Ktspring(ex,ey,kttop,ktbot)
% 
% 
%---------------------------------------------------------------------
% PURPOSE
% Compute the global element tangent stiffness matrix for the soil 
% supporting the pile. It is taking into account that the stiffness of the
% soil springs varies linearly along the considered pile segment (element).
% 
% INPUT:  ex      = [x1 x2]
%         ey      = [y1 y2]       element node coordinates
%         kttop   : Soil tangent stiffness at the top of the element
%         ktbot   : Soil tangent stiffness at the bottom of the element
%
% OUTPUT: Ktspring: element tangent stiffness matrix (6 x 6)
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    11.07.2007
% 
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1) ];
  L     =   sqrt(b'*b);  n=b/L;

  
  Ktetop = kttop * [ 0   0          0          0   0           0;
                     0   2/7*L      1/28*L^2   0   9/140*L    -1/60*L^2;
                     0   1/28*L^2   1/168*L^3  0   1/70*L^2   -1/280*L^3; 
                     0   0          0          0   0           0;
                     0   9/140*L    1/70*L^2   0   3/35*L     -1/60*L^2;
                     0  -1/60*L^2  -1/280*L^3  0  -1/60*L^2    1/280*L^3 ];

               
  Ktebot = ktbot * [ 0   0          0          0    0          0;
                     0   3/35*L     1/60*L^2   0    9/140*L   -1/70*L^2;
                     0   1/60*L^2   1/280*L^3  0    1/60*L^2  -1/280*L^3; 
                     0   0          0          0    0          0;
                     0   9/140*L    1/60*L^2   0    2/7*L     -1/28*L^2;
                     0  -1/70*L^2  -1/280*L^3  0   -1/28*L^2   1/168*L^3 ];

  Kte    = Ktetop + Ktebot;      
            
            
  G      = [ n(1) n(2)  0    0    0   0;
            -n(2) n(1)  0    0    0   0;
              0    0    1    0    0   0;
              0    0    0   n(1) n(2) 0;
              0    0    0  -n(2) n(1) 0;
              0    0    0    0    0   1 ];
  
  Ktspring1 = G'*Kte*G;
