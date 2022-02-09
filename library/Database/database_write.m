function database_write(output,settings,loads,data,plots,save_name,pile,i,soil)
%% WIKINGER DATABASE MODULE
%--------------------------------------------------------------------------
% CHANGE LOG
% 2010.11.04    MMOL    Programming
% 2011.0XXXX    JALY    Streamlining code
% 2011.06.23    MMOL    Re-programming to match Dantysk OWF
% 2013.09.30    MUOE    Re-programming to match Wikinger OWF
% Access MySQL-database
% mysql('open','212.202.161.166','wkdb_viewer','ibdvvdd'); % ('server','username','password')
%mysql('open','sql.ims-ing.de','wkdb_viewer','ibdvvdd'); % ('server','username','password')
% mysql('open','DKLYCOPILOD1','wkdb_viewer','ibdvvdd'); % ('server','username','password')
mysql('open',settings.db_server,settings.db_user,settings.db_pass); % ('open','server','username','password')
mysql(['use ',settings.db_name]); % name of database

%%
% ----- Assumptions to be checked in db before use ------------------------------------
% -------------------------------------------------------------------------

% scour.ORD = 0;   % set to zero, not used for Merkur
% scour.local = 0; % set to zero, scour protection is applied for all locations at Merkur

% origin of local coordinate system (in CO-SPIN) is defined at initial mudline (defined water depth)

% pile head is defined at initial mudline for positions at Merkur:
%    pile.head                   = 0; % [m VREF] pile head is set to be at mudline (internal COSPIN definition)

% pile diameter and steel properties of that can, which enters the seabed
% are used for the whole pile

% Definiton of dimensions for determination of critical pile length, now:
%    pile.length_start           = 25; % [m] 
%    pile.length_end             = 50; % [m]
%    pile.length_inc             = 0.5; % [m]

% Definition of psf applied on soil parameters in ULS case, if applicable
% (now: drained: SF.drained = 1.15 ; undrained: SF.undrained = 1.25)

% loads are interpolated to pile.head elevation

% shear forve Fxy is defined positiv and moment Mxy is defined negativ,
% since it is an monopile design at Merkur

% sign of axial loading for compression and tension 
% (now: compression = negativ, tension = positiv) 

%--------------------------------------------------------------------------
%% Database unique id's
%--------------------------------------------------------------------------
location = ['''',data.id,'''']; % name of location
rev_global = data.revision.global; % global revision no. for specified location to be used

if rev_global == -1 % saved results only if a global revision (configuration) is used
    rev_global = data.revision.output;
% else % a global revision (configuration) is used -> results will be saved
%--------------------------------------------------------------------------
% Check revision no. (table: soil_results_geo)
%--------------------------------------------------------------------------
end
% check, if specified global revision is available for specified location
table = 'soil_results_geo';
[rev]       = mysql(['select rev from ',table,' where id=',location]);
if ismember(rev_global,str2double(rev)) == 0 % if specified global revision doesn't exist
    if rev_global<10
        rev_global1          = ['''0',num2str(rev_global),'''']; % revision no. to be used
        rev_global2          = ['0',num2str(rev_global)]; % revision no. to be used by perl
    else
        rev_global1          = ['''',num2str(rev_global),'''']; % revision no. to be used
        rev_global2          = num2str(rev_global); % revision no. to be used by perl
    end 
    % axial capacity
    if (plots.res_vs_pilelength == 1  && i == length(pile.length)) % save results for axial capacity
        axial_capacity_comp  = round(output.Rd(end,2)); % axial capacity in compression [kN], incl. psf
        axial_capacity_tens  = round(output.Rd(end,1)); % axial capacity in tnesion [kN], incl. psf
%       axial_util_ratio= round(output.axial_util_ratio(end,1)*100) / 100;   % axial utilisation ratio
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,axial_capa_comp,axial_capa_ten) VALUES (',location,',',rev_global1...        %axial_util_ratio
            ,',',num2str(axial_capacity_comp),',',num2str(axial_capacity_tens),');']; %,',',num2str(axial_util_ratio)
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo axial_plot ',save_name{1}]);
        cd(oldFolder);
        delete(save_name{1});
    end
    % critical pile length
    if (plots.pilehead_vs_length == 1 && i == length(pile.length)) % save results for critical pile length
        crit_pile_length  = ceil(output.pile_length_deter.crit_length(1,1)*10)/10; % critical pile length (10% criterium)
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,crit_pile_length'...
            ') VALUES (',location,',',rev_global1...
            ,',',num2str(crit_pile_length),');'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo crit_pile_length_plot ',save_name{2}]);
        (plots.pilehead_vs_length == 1 && i == length(pile.length)) 
        cd(oldFolder);
        delete(save_name{2});
    end
    % horizontal deflection
    if plots.deflection_plot == 1 % save results for deflection
        defl_pile_head  = round(output.hor_defl(1,1) * 10000) / 10000; % pile head deflection
        defl_pile_tip   = round(output.hor_defl(end,1) * 10000) / 10000; % pile tip deflection
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,hori_defl_mudline,'...
            'hori_defl_pile_tip) VALUES (',location,',',rev_global1...
            ,',',num2str(defl_pile_head),',',num2str(defl_pile_tip),');'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo hori_defl_plot ',save_name{3}]);
        cd(oldFolder);
        delete(save_name{3});
    end
    % permanent rotation ULS loads
    if (plots.permanent_rot_def == 1 && loads.n_cycles == 100) % save results for permanent rotation
        rot_pile_head  = round(output.perm_rot*(180/pi)*1000) / 1000; % pile head permanent rotation
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,perm_rot_mudline'...
            ') VALUES (',location,',',rev_global1...
            ,',',num2str(rot_pile_head),');'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo perm_rot_mudline_plot ',save_name{5}]);
        cd(oldFolder);
        delete(save_name{5});
    end
    % permanent rotation operational loads
    if (plots.permanent_rot_def == 1 && loads.n_cycles == 10000) % save results for permanent rotation
        rot_pile_head  = round(output.perm_rot*1000) / 1000; % pile head permanent rotation
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,perm_rot_mudline_oper'...
            ') VALUES (',location,',',rev_global1...
            ,',',num2str(rot_pile_head),');'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo perm_rot_mudline_oper_plot ',save_name{5}]);
        cd(oldFolder);
        delete(save_name{5});
    end
    % horizontal utilisation ratio (global failure)
    if plots.utilization_ratio == 1 % save results for permanent rotation
        [SF reduction]  = factors(loads); % include partial safety factors
        lat_capacity  = round(output.PU / SF.R_ult); % lateral capacity [kN], incl. psf
        lat_util_ratio= round(output.UR_global*100) / 100;   % lateral utilisation ratio
        mysqlstr        = ['INSERT INTO ',table,'(id,rev,lat_capa,'...
            'lat_util_ratio) VALUES (',location,',',rev_global1...
            ,',',num2str(lat_capacity),',',num2str(lat_util_ratio),');'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo lat_plot ',save_name{4}]);
        cd(oldFolder);
        delete(save_name{4});
    end
else    % if specified global revision already exists
    if rev_global<10
        rev_global1          = ['''0',num2str(rev_global),'''']; % revision no. to be used
        rev_global2          = ['0',num2str(rev_global)]; % revision no. to be used by perl
    else
        rev_global1          = ['''',num2str(rev_global),'''']; % revision no. to be used
        rev_global2          = num2str(rev_global); % revision no. to be used by perl
    end 
    % axial capacity
    if settings.axial_loading==1 % save results for axial capacity
        axial_capacity_comp  = round(output.Rd(end,2)); % axial capacity in compression [kN], incl. psf
        axial_capacity_tens  = round(output.Rd(end,1)); % axial capacity in tension [kN], incl. psf
         axial_util_ratio= round(output.axial_util_ratio(end,1) * 100) / 100;   % axial utilisation ratio
        mysqlstr        = ['Update ',table,' set axial_capa_comp=',...
            num2str(axial_capacity_comp),', axial_capa_ten=',num2str(axial_capacity_tens),',axial_util_ratio=',num2str(axial_util_ratio,'%0.2f')...     %',axial_util_ratio=',num2str(axial_util_ratio),
            ' where id=',location,' and rev=',rev_global1,';'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo axial_plot ',save_name{1}]);
        cd(oldFolder);
        delete(save_name{1});
    end
    % critical pile length
    if (plots.pilehead_vs_length == 1 && i == length(pile.length))% save results for critical pile length
        crit_pile_length  = ceil(output.pile_length_deter.crit_length(1,1)*10)/10; % critical pile length (10% criterium)
        mysqlstr        = ['Update ',table,' set crit_pile_length=',...
            num2str(crit_pile_length),' where id=',location,' and rev=',rev_global1,';'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo crit_pile_length_plot ',save_name{2}]);
        cd(oldFolder);
        delete(save_name{2});
    end
    % horizontal deflection
    if plots.deflection_plot == 1 % save results for deflection
        defl_pile_head  = round(output.hor_defl(1,1) * 10000) / 10000; % pile head deflection
        defl_pile_tip   = round(output.hor_defl(end,1) * 10000) / 10000; % pile tip deflection
        mysqlstr = ['Update ',table,' set hori_defl_mudline=',num2str(defl_pile_head),...
            ', hori_defl_pile_tip=',num2str(defl_pile_tip),' where id=',location,' and rev=',rev_global1,';'];
        mysql([mysqlstr])
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo hori_defl_plot ',save_name{3}]);
        cd(oldFolder);
        delete(save_name{3});
    end
    % permanent rotation ULS loads
    if (plots.permanent_rot_def == 1 && loads.n_cycles == 100) % save results for permanent rotation
        rot_pile_head  = round(output.perm_rot*(180/pi)*1000) / 1000; % pile head permanent rotation
        mysqlstr        = ['Update ',table,' set perm_rot_mudline='...
            ,num2str(rot_pile_head),' where id=',location,' and rev=',rev_global1,';'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo perm_rot_mudline_plot ',save_name{5}]);
        cd(oldFolder);
        delete(save_name{5});
    end  
    % permanent rotation operational loads
    if (plots.permanent_rot_def == 1 && loads.n_cycles == 10000)% save results for permanent rotation
        rot_pile_head  = round(output.perm_rot*(180/pi)*1000) / 1000; % pile head permanent rotation
        mysqlstr        = ['Update ',table,' set perm_rot_mudline_oper='...
            ,num2str(rot_pile_head),' where id=',location,' and rev=',rev_global1,';'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo perm_rot_mudline_oper_plot ',save_name{5}]);
        cd(oldFolder);
        delete(save_name{5});
    end  
    % horizontal utilisation ratio (global failure)
    if plots.load_deflection == 1 && soil.psf==0% save results for permanent rotation
        %[SF reduction]  = factors(loads); % include partial safety factors
        lat_capacity  = output.H_load_10_bsh;%round(output.PU / SF.R_ult); % lateral capacity [kN], incl. psf
        lat_util_ratio = output.UR_lat_bsh;%round(output.UR_global*100) / 100;   % lateral utilisation ratio
        mysqlstr        = ['Update ',table,' set lat_capa=',...
            num2str(lat_capacity),',lat_util_ratio=',num2str(lat_util_ratio),...
             ' where id=',location,' and rev=',rev_global1,';'];
        mysql(mysqlstr);
        %change to local folder including perl script, run perl script, change back and delete tempfile:
        oldFolder = cd('library\Database\perl\prj\ao1db\soil');
        system(['insert-file-geo.pl ',data.id,' ',rev_global2,' soil_results_geo lat_plot ',save_name{4}]);
        cd(oldFolder);
        delete(save_name{4});
    end
end

% Close MySQL-database
mysql('close')

