function qespring = qspring(ex,ey,kstop,ksbot,u,t,i)
% 
% 
%---------------------------------------------------------------------
% PURPOSE
% Compute the internal forces in a pile segment due to the contributions 
% from the soil springs. It is taking into account that the stiffness of
% the soil springs varies linearly along the considered pile segment 
% (element).
% 
% INPUT:  ex      = [x1 x2]
%         ey      = [y1 y2]       Element node coordinates
%         kstop   : Soil secant stiffness at the top of the element
%         ksbot   : Soil secant stiffness at the bottom of the element
%         u       : Global displacement vector
%         t       : Each row in t describes the global dof associated with
%                   Edof
%         i       : Counter referring to element number
%
% OUTPUT: qespring: Internal force vector due to springs (6 x 1)
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    13.07.2007   Programming
%                   AHA/LA 14.08.2007   Correcting - Transformation matrix
%                                       G and G' are multiplied Kse instead
%                                       of only G
% 
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1) ];
  L     =   sqrt(b'*b);  n=b/L;

  ue = u(t(i,:)');
  
  Ksetop = kstop * [ 0   0          0          0   0           0;
                     0   2/7*L      1/28*L^2   0   9/140*L    -1/60*L^2;
                     0   1/28*L^2   1/168*L^3  0   1/70*L^2   -1/280*L^3; 
                     0   0          0          0   0           0;
                     0   9/140*L    1/70*L^2   0   3/35*L     -1/60*L^2;
                     0  -1/60*L^2  -1/280*L^3  0  -1/60*L^2    1/280*L^3 ];

               
  Ksebot = ksbot * [ 0   0          0          0    0          0;
                     0   3/35*L     1/60*L^2   0    9/140*L   -1/70*L^2;
                     0   1/60*L^2   1/280*L^3  0    1/60*L^2  -1/280*L^3; 
                     0   0          0          0    0          0;
                     0   9/140*L    1/60*L^2   0    2/7*L     -1/28*L^2;
                     0  -1/70*L^2  -1/280*L^3  0   -1/28*L^2   1/168*L^3 ];

  Kse    = Ksetop + Ksebot;      
            
  G      = [ n(1) n(2)  0    0    0   0;
            -n(2) n(1)  0    0    0   0;
              0    0    1    0    0   0;
              0    0    0   n(1) n(2) 0;
              0    0    0  -n(2) n(1) 0;
              0    0    0    0    0   1 ];
  
  qespring = G'*Kse*G*ue;
  
% Dette skal med i beregningen af internal force vector from spring, i.e. in qespring.
% if abs(q(i)) > pu
%         q(i) = pu;
%       end

