function database_write_tz_curve(settings,element,loads,data,t,z,j,plots,table)
%% WIKINGER DATABASE MODULE
%--------------------------------------------------------------------------
% CHANGE LOG
% 2016.03.18    THHM    Programming
% 2019.06.20    GINI    Programming for SNA database

%% Dialog to ensure, if data shall really be saved in the database
if (plots.pilehead_vs_length == 1 && j ~= length(pile.length)) % no dialog for each incremental pile length within critical pile length determination
    SaveSwitch = 0; %Off
elseif (plots.res_vs_pilelength == 1  && j ~= length(pile.length)) % no dialog for each incremental pile length within axial pile length determination
    SaveSwitch = 0; %Off
else            
%     Reply = questdlg('Do you want to save p-y curves in the database?','Save p-y curves to database','Yes','No','Yes');
%     if(strcmp(Reply,'Yes')==1)
        SaveSwitch = 1; %On
%     else
%         SaveSwitch = 0; %Off
%     end
end

%% Access MySQL-database
if(SaveSwitch == 1)
    
    disp('Exporting tz-curves to database')

    % Access MySQL-database
    mysql('open',settings.db_server,settings.db_user,settings.db_pass); % ('open','server','username','password')
    mysql(['use ',settings.db_name]); % name of database



    %--------------------------------------------------------------------------
    %% Database unique id's
    %--------------------------------------------------------------------------
    id = ['''',data.id,'''']; % name of id
    rev_global = data.revision.global; % global revision no. for specified id to be used

    if rev_global == -1 % saved results only if a global revision (configuration) is used
    else
     % a global revision (configuration) is used -> results will be saved
    %--------------------------------------------------------------------------
    % Check revision no. (table: results_geo)
    %--------------------------------------------------------------------------

        % check, if specified global revision is available for specified id
        [rev]       = mysql(['select rev from ',table,' where id=',id]);
        if rev_global<10
            rev_global          = ['''0',num2str(rev_global),'''']; % revision no. to be used
        else
            rev_global          = ['''',num2str(rev_global),'''']; % revision no. to be used
        end

        % if specified global revision exists for this id -> delete all 
        % tz-entries for this id and revision (necessary, if update of
        % tz-curves shall be done using a different discretisation of pile):
        if ismember(data.revision.global,str2double(rev)) > 0 
            mysqlstr = ['DELETE FROM ',table,' where id=',id,' and rev=',rev_global,' and stat_cycl=''',loads.static_cyclic,''';'];
            mysql(mysqlstr);
        end
        %%%%%% t_av is calculated as a weighted average with element thickness between t.top and t.bottom when these two are different in the same node due to a boundary between two different soil layers

        t_av=zeros(length(t.top),length(t.top(1,:)));
        for i =2:element.nelem-1
        for j=1:length(t.top(1,:))
        if t.top(i,j)<t.bottom(i-1,j)||t.top(i,j)>t.bottom(i-1,j)
        t_av(i,j) = (t.top(i,j).*(abs(element.level(i,2))-abs(element.level(i,1)))+t.bottom(i-1,j).*(abs(element.level(i,1))-abs(element.level(i-1,1))))./((abs(element.level(i,2))-abs(element.level(i,1)))+(abs(element.level(i,1))-abs(element.level(i-1,1)))); % maximum of 20 p values is allowed
        else
        t_av(i,j) = t.top(i,j);
        end
        end
        end

        for i =1:element.nelem-1
            % element top nodes
            depth              = abs(element.level(i,1));
            %depth_top_influence     = abs(element.level(i,1));
            %depth_bot_influence     = abs(element.level(i,1) + (element.level(i,2)-element.level(i,1)) / 2);
            model_tz                = element.model_axial{i};
            stat_cycl               = loads.static_cyclic;
            value_t_top             = zeros(1,15); % maximum of 15 t values is allowed
            value_z_top             = zeros(1,15);
            value_t_top(1,1:size(t_av,2)) = t_av(i,:); % maximum of 15 t values is allowed
            value_z_top(1,1:size(z.top,2)) = z.top(i,:);

            mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   't_z_value,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,model_tz,stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth),','''];
                   
                

            % insert p values
        mysqlstr = [mysqlstr_ini,'t'''];
        for j = 1:15 % t_top_values
            mysqlstr = [mysqlstr,',',num2str(value_t_top(1,j))];
        end
        mysqlstr_ttop = [mysqlstr,',''',model_tz,''',''', stat_cycl,''');'];
        
       
        
        % z_top values
        mysqlstr = [mysqlstr_ini,'z'''];
        for j = 1:15 
            mysqlstr = [mysqlstr,',',num2str(value_z_top(1,j))];
        end
        mysqlstr_ztop = [mysqlstr,',''',model_tz,''',''', stat_cycl,''');'];
            
          
            
            mysql(mysqlstr_ttop);
           % mysql(mysqlstr_tbot);
            mysql(mysqlstr_ztop);
           % mysql(mysqlstr_zbot);
            
         
        end
        
         % element bot node (only bot element)
        depth_2              = abs(element.level(length(element.level)-1,2));
        %depth_top_influence     = abs(element.level(i,2) - (element.level(i,2)-element.level(i,1)) / 2);
        %depth_bot_influence     = abs(element.level(i,2));
        model_tz                = element.model_axial{length(element.level)-1};
        stat_cycl               = loads.static_cyclic;
        value_t_bot             = zeros(1,15);
        value_z_bot             = zeros(1,15);
        value_t_bot(1,1:size(t.bottom,2)) = t.bottom(length(t.bottom),:);
        value_z_bot(1,1:size(z.bottom,2)) = z.bottom(length(z.bottom),:);
           
           
           mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   't_z_value,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,model_tz,stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth_2),','''];
        

            % insert p values
        mysqlstr = [mysqlstr_ini,'t'''];
        for j = 1:15 % t_top_values
            mysqlstr = [mysqlstr,',',num2str(value_t_bot(1,j))];
        end
        mysqlstr_tbot = [mysqlstr,',''',model_tz,''',''', stat_cycl,''');'];
        
       
        
        % z_top values
        mysqlstr = [mysqlstr_ini,'z'''];
        for j = 1:15 
            mysqlstr = [mysqlstr,',',num2str(value_z_bot(1,j))];
        end
        mysqlstr_zbot = [mysqlstr,',''',model_tz,''',''', stat_cycl,''');'];
            
          
            
           %mysql(mysqlstr_ttop);
           mysql(mysqlstr_tbot);
           %mysql(mysqlstr_ztop);
           mysql(mysqlstr_zbot);
        
    end


    % Close MySQL-database
    mysql('close')
    disp('Finished exporting tz-curves')

end
