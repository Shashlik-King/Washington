function [q,ks_tot] = internal_forces(nelem,u,element,pile,reduction,plug_unplug,loads,Ex,Ey,t,settings,q,ks_tot,A,fs,qp,ii)
%--------------------------------------------------------------------------
% PURPOSE
% Compute the internal forces [kN] in the pile segment by applying t-z curves for the soil.
% 
% INPUT:  element       : Element information about soil 
%         element.degradation: Degradation of t-z curve to take cyclic into
%                              account 
%         pile          : Pile information for each element
%         z_topbottom   : Vertical displacement in top and bottom of
%                         element
%         i             : Counter referring to element number
%
% OUTPUT: q       : Total internatl forces [kN]
%         ks_tot  : Secant soil stiffness at the top and bottom of all elements [kN/m/m]
%
% CODE            : EBGX
% APPROVED        : 

% LAST MODIFIED   : EBGX    24.11.2016   Programming
%--------------------------------------------------------------------------
        for i = 1:(element.nelem-1)
          
            % Internal forces calculated and assembling
            % Description of variables after the N-R procedure

            z_topbottom         = [u(3*i-1) u(3*i+2)]; % vertical disp. at the top and bottom of the pile segment [m]
            
            [kstop ksbot]       = secspringstiff_axial(settings,element,pile,reduction,plug_unplug,loads,z_topbottom,A,fs,i,ii);
            
            qespring            = qspring_axial(Ex(i,:),Ey(i,:),kstop,ksbot,u,t,i);
            
			if settings.beam_theory
                qebeam              = qbeam2t(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
            else
                qebeam              = qbeam2e(Ex(i,:),Ey(i,:),element.ep(i,:),u,t,i);
			end
			
            qelem               = qebeam + qespring;
            q(t(i,:))           = q(t(i,:)) + qelem;
                
            ks_tot(i,1)=kstop;
            ks_tot(i,2)=ksbot;
            
            if i == nelem && strcmp(settings.analysis_loading(ii),'Comp') == 1
                zbot = u(3*i+2);         % Horizontal disp. at the bottom of the last pile element [m]
                Q_end    = Qz_spring(element,pile,plug_unplug,zbot,qp,i);
                q(i*3+2) = q(i*3+2)+Q_end;
            end
        end 
        return 