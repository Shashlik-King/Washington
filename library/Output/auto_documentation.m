%% automatic documentation of calculation settings and results/variables
function auto_documentation(data,settings,soil,plots,loads);
%--------------------------------------------------------------------------
% PURPOSE
% Automatic generation of documentation. Generation of ASCI file of Run_COSPIN and .mat file of CO-SPIN variables
%
% INPUT:  data		    : structure containing project input data
%         settings      : structure containing settings for calculations
%         soil          : structure containing soil input data
%         plots         : structure containing settings for plots
%         loads         : structure containing load input
%
% OUTPUT:
%
% CODE            : THHM
% APPROVED        : DATY

% LAST MODIFIED   : 
%--------------------------------------------------------------------------

%% create subfolder for documentation, if not existent
dir     = 'documentation';
if exist(dir,'file') ~= 7   % if subfolder 'Data' does not exist -> create it
    [check, tmp1, tmp2]  =   mkdir(dir);
    if check ~= 1, error('subfolder for documentation ',dir1,' couldn''t be created'),return, end
end

%% create filename and save data
if settings.PSI
    filename_tmp    = [data.location,'_PSI'];
    copyfile('run_COSPIN.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
    save([dir,'\',filename_tmp,'_workspace']);
else
	if settings.axial_loading
		filename_tmp    = [data.location,'_Axial_Capacity'];
		copyfile('run_COSPIN_axial_capacity.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end   
	if plots.pilehead_vs_length
		filename_tmp    = [data.location,'_Critical_Pile_Length'];
		copyfile('run_COSPIN_critical_length.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end 
	if plots.deflection_plot
		filename_tmp    = [data.location,'_Deflection'];
		copyfile('run_COSPIN_lat_deflection.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
	if plots.utilization_ratio
		filename_tmp    = [data.location,'_Lat_Soil_Utilisation'];
		copyfile('run_COSPIN_lat_UR.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
	if plots.deflection_bundle
		filename_tmp    = [data.location,'_Deflection_Bundle'];
		copyfile('run_COSPIN.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
	if plots.toe_shear_graph
		filename_tmp    = [data.location,'_Toe_Shear_Graph'];
		copyfile('run_COSPIN.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
	if plots.permanent_rot_def
		filename_tmp    = [data.location,'_Permanent_Rotation_',num2str(loads.n_cycles),'_cycles'];
		copyfile('run_COSPIN_perm_rotation.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
	if plots.load_deflection && soil.psf==1
		filename_tmp    = [data.location,'_Load_deflection'];
		copyfile('run_COSPIN_load_deflection_UR_DNVGL.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
    
    elseif plots.load_deflection && soil.psf==0
		filename_tmp    = [data.location,'_Load_deflection'];
		copyfile('run_COSPIN_load_deflection_UR_BSH.m' , [dir,'\',filename_tmp,'_run_COPSIN.dat']);
		save([dir,'\',filename_tmp,'_workspace']);
	end
end

return