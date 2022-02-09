function qespring = qspring_axial(ex,ey,kstop,ksbot,u,t,i)
% 
% 
%---------------------------------------------------------------------
% PURPOSE
% Compute the internal forces in a pile segment due to the contributions 
% from the soil springs. It is taken into account that the spring stiffness 
% is different in the top and bottom of the element. Meaning that the total 
% stiffness divided 50/50 between top and bottom.  
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
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    24.10.2016   Programming
%                  
%-------------------------------------------------------------

  b     =   [ ex(2)-ex(1); ey(2)-ey(1); ];
  L     =   sqrt(b'*b);  n=b/L;

  ue = u(t(i,:)');
  
  Ksetop = 0.5*kstop * [ L/3  0  0   L/6   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0; 
                     L/6  0  0   L/3   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0 ];

               
  Ksebot = 0.5*ksbot * [ L/3  0  0   L/6   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0; 
                     L/6  0  0   L/3   0   0;
                     0    0  0   0     0   0;
                     0    0  0   0     0   0 ];

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

