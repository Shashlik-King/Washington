function [Kt Kts] = Kt_axial(u,element,pile,reduction,plug_unplug,loads,Ex,Ey,nelem,settings,t,Kts,Kt,A,fs,ii)  
%--------------------------------------------------------------------------
% PURPOSE
% Compute the tangent stiffness matrix for the pile by applying t-z curves for the soil.
% 
% INPUT:  element       : Element information about soil 
%         element.degradation: Degradation of t-z curve to take cyclic into
%                              account 
%         pile          : Pile information for each element
%         z_topbottom   : Vertical displacement in top and bottom of
%                         element
%         i             : Counter referring to element number
%
% OUTPUT: Kt   : Global stiffness matrix  [kN/m]
%         Kt   : Global stiffness matrix, only for the soil part [kN/m]
%
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    24.11.2016   Programming
%--------------------------------------------------------------------------
for i=1:(element.nelem-1)
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
end 