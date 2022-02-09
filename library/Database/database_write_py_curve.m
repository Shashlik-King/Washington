function database_write_py_curve(settings,element,loads,data,p,y,j,plots,table)
%% WIKINGER DATABASE MODULE
%--------------------------------------------------------------------------
% CHANGE LOG
% 2016.03.18    THHM    Programming
% 2019.06.20    GINI    Programming for SNA database
% 2019.10.21    ASSV    Programming for ARCADIS database

%% Dialog to ensure, if data shall really be saved in the database - Not used 
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
    
    disp('Exporting py-curves to database')

    % Access MySQL-database
    mysql('open',settings.db_server,settings.db_user,settings.db_pass); % ('open','server','username','password')
    mysql(['use ',settings.db_name]); % name of database



    %--------------------------------------------------------------------------
    %% Database unique id's - Verify revisions
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
        % py-entries for this id and revision (necessary, if update of
        % py-curves shall be done using a different discretisation of pile):
        if ismember(data.revision.global,str2double(rev)) > 0 
            mysqlstr = ['DELETE FROM ',table,' where id=',id,' and rev=',rev_global,' and stat_cycl=''',loads.static_cyclic,''';'];
            mysql(mysqlstr);
        end

        %% Adjust p data
        % p_av is calculated as a weighted average with element thickness 
        % between p.top and p.bottom when these two are different in the same node
        % due to a boundary between two different soil layers
        
        p_av=zeros(length(p.top),length(p.top(1,:)));
       
        for i =2:element.nelem-1
        for j=1:length(p.top(1,:))
        if p.top(i,j)<p.bottom(i-1,j)||p.top(i,j)>p.bottom(i-1,j)
        p_av(i,j) = (p.top(i,j).*(abs(element.level(i,2))-abs(element.level(i,1)))+p.bottom(i-1,j).*(abs(element.level(i,1))-abs(element.level(i-1,1))))./((abs(element.level(i,2))-abs(element.level(i,1)))+(abs(element.level(i,1))-abs(element.level(i-1,1)))); 
        else
        p_av(i,j) = p.top(i,j);
        end
        end
        end       
        %% Generate strings and save into database until nelem-1
        for i =1:element.nelem-1
            % element top nodes
            depth              = abs(element.level(i,1));
            %depth_top_influence     = abs(element.level(i,1));
            %depth_bot_influence     = abs(element.level(i,1) + (element.level(i,2)-element.level(i,1)) / 2);
            model_py                = element.model_py{i};
            stat_cycl               = loads.static_cyclic;
            value_p_top             = zeros(1,20); % maximum of 20 p values is allowed
            value_y_top             = zeros(1,20);
            value_p_top(1,1:size(p_av,2)) = p_av(i,:); % maximum of 20 p values is allowed
            value_y_top(1,1:size(y.top,2)) = y.top(i,:);

            mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   'p_y_value,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,model_py,stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth),','''];
                   
                

            % insert p values
        mysqlstr = [mysqlstr_ini,'p'''];
        for j = 1:20 % p_top_values
            mysqlstr = [mysqlstr,',',num2str(value_p_top(1,j))];
        end
        mysqlstr_ptop = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];
        
       
        
        % y_top values
        mysqlstr = [mysqlstr_ini,'y'''];
        for j = 1:20 
            mysqlstr = [mysqlstr,',',num2str(value_y_top(1,j))];
        end
        mysqlstr_ytop = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];
            
          
            
            mysql(mysqlstr_ptop);
           % mysql(mysqlstr_pbot);
            mysql(mysqlstr_ytop);
           % mysql(mysqlstr_ybot);
                           
        end       
        %% Generate strings and save into database only bottom elem
        depth_2              = abs(element.level(length(element.level)-1,2));
        %depth_top_influence     = abs(element.level(i,2) - (element.level(i,2)-element.level(i,1)) / 2);
        %depth_bot_influence     = abs(element.level(i,2));
        model_py                = element.model_py{length(element.level)-1};
        stat_cycl               = loads.static_cyclic;
        value_p_bot             = zeros(1,20);
        value_y_bot             = zeros(1,20);
        value_p_bot(1,1:size(p.bottom,2)) = p.bottom(length(p.bottom),:);
        value_y_bot(1,1:size(y.bottom,2)) = y.bottom(length(y.bottom),:);
           
           
           mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   'p_y_value,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,model_py,stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth_2),','''];
            % insert p values
        mysqlstr = [mysqlstr_ini,'p'''];
        for j = 1:20 % p_top_values
            mysqlstr = [mysqlstr,',',num2str(value_p_bot(1,j))];
        end
        mysqlstr_pbot = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];
        
       
        
        % y_top values
        mysqlstr = [mysqlstr_ini,'y'''];
        for j = 1:20 
            mysqlstr = [mysqlstr,',',num2str(value_y_bot(1,j))];
        end
        mysqlstr_ybot = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];
            
          
            
           % mysql(mysqlstr_ptop);
            mysql(mysqlstr_pbot);
            %mysql(mysqlstr_ytop);
           mysql(mysqlstr_ybot);
           
           
        
    end


    % Close MySQL-database
    mysql('close')
    disp('Finished exporting py-curves')

end
