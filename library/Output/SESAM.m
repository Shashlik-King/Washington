% creation of Excel file with linear soil springs for SESAM
function [Stiff_node_py,Stiff_node_mtheta] = SESAM(node,scour,pile,data,element,loads,u,ustep,settings,output)
%--------------------------------------------------------------------------
% PURPOSE
% 
%
% INPUT:  
%
% OUTPUT:
%
% CODE            : 
% APPROVED        : 

% LAST MODIFIED   : 
%--------------------------------------------------------------------------

%% create subfolder for documentation, if not existent
dir     = 'SESAM';
if exist(dir,'file') ~= 7   % if subfolder 'Data' does not exist -> create it
    [check, ~, ~]  =   mkdir(dir);
    if check ~= 1, error('subfolder for documentation ',dir1,' couldn''t be created'),return, end
end

nelem=length(element.level);%number of elements
Stiff_node_py=zeros(nelem,1);%Assignment of vector for distributed load stiffness at nodes
Stiff_node_mtheta=zeros(nelem,1);%Assignment of vector for distributed moment stiffness at nodes
for c = 1:nelem-1
%     upy_top=abs(ustep(c*3-2,settings.n_max));%lateral displacement at top of element
%     upy_bot=abs(ustep((c+1)*3-2,settings.n_max)); %lateral displacement at bottom of element
    teta_top=abs(ustep(c*3,settings.n_max));%rotation at top of element
    teta_bot=abs(ustep(c*3+3,settings.n_max));%rotation at bottom of element
	
    
	Stiff_py(c,1)=abs(output.p_UR(c,1)/ustep(c*3-2,settings.n_max));%Secant stiffness of distributed load springs at top of element. Unit kN/m/m.
	Stiff_py(c,2)=abs(output.p_UR(c,2)/ustep(c*3+1,settings.n_max));%Secant stiffness of distributed load springs at bottom of element. Unit kN/m/m.
    
    Stiff_node_py(c)=Stiff_node_py(c)+(2/6*Stiff_py(c,1)+1/6*Stiff_py(c,2))*(element.level(c,1)-element.level(c,2))*1000;%distributed load contribution to top node in element. Conversion to unit (N/m)
    Stiff_node_py(c+1)=Stiff_node_py(c+1)+(1/6*Stiff_py(c,1)+2/6*Stiff_py(c,2))*(element.level(c,1)-element.level(c,2))*1000;%distributed load contribution to bottom node in element. Conversion to unit (N/m)
    
    if settings.mteta
% 	Stiff_mtheta(c,1)=abs(output.m_UR(c,1)/teta_bot(c*3,settings.n_max));%secant stiffness of distributed moment spring at top of element. Unit kNm/m/rad
% 	Stiff_mtheta(c,2)=abs(output.m_UR(c,2)/teta_bot(c*3+3,settings.n_max));%secant stiffness of distributed moment spring at bottom of element. Unit kNm/m/rad
    
    Stiff_mtheta(c,1)=abs(output.m_UR(c,1)/teta_top);%secant stiffness of distributed moment spring at top of element. Unit kNm/m/rad FKMV proposal
	Stiff_mtheta(c,2)=abs(output.m_UR(c,2)/teta_bot);%secant stiffness of distributed moment spring at bottom of element. Unit kNm/m/rad FKMV proposal
    
    Stiff_node_mtheta(c)=Stiff_node_mtheta(c)+(2/6*Stiff_mtheta(c,1)+1/6*Stiff_mtheta(c,2))*(element.level(c,1)-element.level(c,2))*1000;%distributed moment contribution to top node in element. Conversion to unit (Nm/rad)
    Stiff_node_mtheta(c+1)=Stiff_node_mtheta(c+1)+(1/6*Stiff_mtheta(c,1)+2/6*Stiff_mtheta(c,2))*(element.level(c,1)-element.level(c,2))*1000;%distributed moment contribution to bottom node in element. Conversion to unit (Nm/rad)
    end
        
end
if settings.toe_shear
   y_topbottom         = [abs(u(end-5,settings.n_max)) abs(ustep(end-2,settings.n_max))]; % Horizontal disp. at the top and bottom of the pile segment [m]
   teta_topbottom      = [abs(u(end-3,settings.n_max)) abs(u(end,settings.n_max))]; % Nodal section rotation at the top and bottom of the pile segment [rad]
   [kstop, ksbot]       = secspringstiff(element,pile,loads,y_topbottom,nelem); % secant stiffness for distributed lateral load 
   Stiff_node_py(end)=Stiff_node_py(end)+(kstop+ksbot)/2*(element.level(end,1)-element.level(end,2))*1000;%adding contribution from base shear to bottom node
    if settings.mteta
		[ksmomtop, ksmombot] = secmomstiff(element,pile,teta_topbottom,y_topbottom,kstop,ksbot,nelem);
        Stiff_node_mtheta(end)=Stiff_node_mtheta(end)+(ksmomtop+ksmombot)/2*(element.level(end,1)-element.level(end,2))*1000;%adding contribution from base moment spring to bottom node
    end
    
end

input=[node.level(1:end-1) Stiff_node_py Stiff_node_mtheta];

excelname = 'saaaaaaggggggggg.xlsx';

 xlswrite(excelname ,input(:,1), data.analysis, 'A3')
 xlswrite(excelname ,input(:,2), data.analysis, 'B3')
 xlswrite(excelname ,input(:,3), data.analysis, 'C3')
 xlswrite(excelname ,cellstr('Depth below mudline'), data.analysis, 'A1')
 xlswrite(excelname ,cellstr('Lateral spring stiffness'), data.analysis, 'B1')
 xlswrite(excelname ,cellstr('Rotational spring stiffness'), data.analysis, 'C1')
 xlswrite(excelname ,cellstr('[m]'), data.analysis, 'A2')
 xlswrite(excelname ,cellstr('[N/m]'), data.analysis, 'B2')
 xlswrite(excelname ,cellstr('[Nm/rad]'), data.analysis, 'C2')
end