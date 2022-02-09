function Ktspring1 = Ktspring_axial(ex,ey,kttop,ktbot)
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
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    18.01.2017
% 
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1) ];
  L     =   sqrt(b'*b);  
  n=b/L;
  
  Ktetop = 0.5*kttop * [ L/3  0  0   L/6   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0; 
                     L/6  0  0   L/3   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0 ];

               
  Ktebot = 0.5*ktbot * [ L/3  0  0   L/6   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0; 
                     L/6  0  0   L/3   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0 ];

%   Ktetop = 0.5*kttop * [ 0    0    0   0   0    0;
%                          0    L/3  0   0   L/6  0;
%                          0    0    0   0   0    0; 
%                          0    0    0   0   0    0;
%                          0    L/6  0   0   L/3  0;
%                          0    0    0   0   0    0 ];
% 
%                
%   Ktebot = 0.5*ktbot * [ 0    0    0   0   0    0;
%                          0    L/3  0   0   L/6  0;
%                          0    0    0   0   0    0; 
%                          0    0    0   0   0    0;
%                          0    L/6  0   0   L/3  0;
%                          0    0    0   0   0    0 ];

  Kte    = Ktetop + Ktebot;      
            
            
  G      = [ n(1) n(2)  0    0    0   0;
            -n(2) n(1)  0    0    0   0;
              0    0    1    0    0   0;
              0    0    0   n(1) n(2) 0;
              0    0    0  -n(2) n(1) 0;
              0    0    0    0    0   1 ];
  
  Ktspring1 = G'*Kte*G;