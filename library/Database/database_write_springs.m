function database_write_springs(settings,element,loads,data,p,y,table,type)
%--------------------------------------------------------------------------
% CHANGE LOG
% 2019.11.21    ASSV    Programming

    

    %% Connect to database
    if strcmp(type,'py')
        disp('Exporting py-curves to database')
    elseif strcmp(type,'tz')
        disp('Exporting tz-curves to database')
    end
    

    % Access MySQL-database
    mysql('close')
    mysql('open',settings.db_server,settings.db_user,settings.db_pass); % ('open','server','username','password')
    mysql(['use ',settings.db_name]); % name of database
    %% Database unique id's - Verify revisions
    %--------------------------------------------------------------------------
    id = ['''',data.id,'''']; % name of id
    if data.revision.global == -1 % saved results only if a global revision (configuration) is used
        disp('Global revision not specified - > SSI springs are saved into the database with data.revision.output')
        rev_global = data.revision.output; % global revision no. for specified id to be used
        rev_check = rev_global;
        % check, if specified global revision is available for specified id
        [rev]       = mysql(['select rev from ',table,' where id=',id]);
        if rev_global<10
            rev_global          = ['''0',num2str(rev_global),'''']; % revision no. to be used
        else
            rev_global          = ['''',num2str(rev_global),'''']; % revision no. to be used
        end
    else    
        rev_global = data.revision.global; % global revision no. for specified id to be used
        % check, if specified global revision is available for specified id
        rev_check = rev_global;
        [rev]       = mysql(['select rev from ',table,' where id=',id]);
        if rev_global<10
            rev_global          = ['''0',num2str(rev_global),'''']; % revision no. to be used
        else
            rev_global          = ['''',num2str(rev_global),'''']; % revision no. to be used
        end
    end


    % if specified global revision exists for this id -> delete all 
    % py-entries for this id and revision (necessary, if update of
    % py-curves shall be done using a different discretisation of pile):
    if ismember(rev_check,str2double(rev)) > 0 
%         mysqlstr = ['DELETE FROM ',table,' where id=',id,' and rev=',rev_global,' and stat_cycl=''',loads.static_cyclic,''';'];
        mysqlstr = ['DELETE FROM ',table,' where id=',id,' and rev=',rev_global];
        mysql(mysqlstr);
    end
    %% Definining spring type values
    if strcmp(type,'py')
        spring.col_value = 'p_y_value';
        spring.col_p = 'p';
        spring.col_y = 'y';
        spring.col_model = 'model_py';
        spring.ncol = 20;
    elseif strcmp(type,'tz')
        spring.col_value = 't_z_value';
        spring.col_p = 't';
        spring.col_y = 'z';
        spring.col_model = 'model_tz';
        spring.ncol = 15;
    end

    %% Correction for nodes with 2 springs
    % p_av is calculated as a weighted average with element thickness 
    % between p.top and p.bottom when these two are different in the same node
    % due to a boundary between two different soil layers

    p_av=zeros(length(p.top),length(p.top(1,:)));

    for i =2:element.nelem-1
        for j=1:length(p.top(1,:))
            if p.top(i,j)<p.bottom(i-1,j)||p.top(i,j)>p.bottom(i-1,j)
                p_av(i,j) = (p.top(i,j).*(abs(element.level(i,2))-...
                    abs(element.level(i,1)))+p.bottom(i-1,j).*...
                    (abs(element.level(i,1))-abs(element.level(i-1,1))))...
                    ./((abs(element.level(i,2))-abs(element.level(i,1)))...
                    +(abs(element.level(i,1))-abs(element.level(i-1,1)))); 
            else
                p_av(i,j) = p.top(i,j);
            end
        end
    end   

    %% Create strings and save top srpings until nelem-1
    for i =1:element.nelem-1
        % element top nodes
        depth              = abs(element.level(i,1));
        model_py                = element.model_py{i};
        stat_cycl               = loads.static_cyclic;
        value_p_top             = zeros(1,spring.ncol); % maximum of 20 p values is allowed
        value_y_top             = zeros(1,spring.ncol);
        value_p_top(1,1:size(p_av,2)) = p_av(i,:); % maximum of 20 p values is allowed
        value_y_top(1,1:size(y.top,2)) = y.top(i,:);

        if strcmp(type,'py')
            mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   spring.col_value,',value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,',spring.col_model,',stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth),','''];
        elseif strcmp(type,'tz')
            mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
                   spring.col_value,',value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                    'value11,value12,value13,value14,value15,',spring.col_model,',stat_cycl) VALUES (',...
                    id,',',...
                    rev_global,',',...
                    num2str(depth),','''];            
        end

        % insert p values
        mysqlstr = [mysqlstr_ini,spring.col_p,''''];
        for j = 1:spring.ncol % p_top_values
            mysqlstr = [mysqlstr,',',num2str(value_p_top(1,j))]; %#ok<*AGROW>
        end
        mysqlstr_ptop = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];

        % y_top values
        mysqlstr = [mysqlstr_ini,spring.col_y,''''];
        for j = 1:spring.ncol 
            mysqlstr = [mysqlstr,',',num2str(value_y_top(1,j))];
        end
        mysqlstr_ytop = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];

        mysql(mysqlstr_ptop);
        mysql(mysqlstr_ytop);                           
    end 
    %% Create strings and save bot srpings for last element
    depth_2              = abs(element.level(length(element.level)-1,2));
    model_py                = element.model_py{length(element.level)-1};
    stat_cycl               = loads.static_cyclic;
    value_p_bot             = zeros(1,spring.ncol);
    value_y_bot             = zeros(1,spring.ncol);
    value_p_bot(1,1:size(p.bottom,2)) = p.bottom(length(p.bottom),:);
    value_y_bot(1,1:size(y.bottom,2)) = y.bottom(length(y.bottom),:);

    if strcmp(type,'py')
        mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
               spring.col_value,',value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                'value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,',spring.col_model,',stat_cycl) VALUES (',...
                id,',',...
                rev_global,',',...
                num2str(depth_2),','''];
    elseif strcmp(type,'tz')
        mysqlstr_ini        = ['INSERT INTO ',table,'(id,rev,depth,'...
               spring.col_value,',value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,'...
                'value11,value12,value13,value14,value15,',spring.col_model,',stat_cycl) VALUES (',...
                id,',',...
                rev_global,',',...
                num2str(depth_2),','''];        
    end
    % insert p values
    mysqlstr = [mysqlstr_ini,spring.col_p,''''];
    for j = 1:spring.ncol % p_top_values
        mysqlstr = [mysqlstr,',',num2str(value_p_bot(1,j))];
    end
    mysqlstr_pbot = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];

    % y_top values
    mysqlstr = [mysqlstr_ini,spring.col_y,''''];
    for j = 1:spring.ncol 
        mysqlstr = [mysqlstr,',',num2str(value_y_bot(1,j))];
    end
    mysqlstr_ybot = [mysqlstr,',''',model_py,''',''', stat_cycl,''');'];

    mysql(mysqlstr_pbot);
    mysql(mysqlstr_ybot);
        % Close MySQL-database
    mysql('close')
    if strcmp(type,'py')
        disp('Finished exporting py-curves')
    elseif strcmp(type,'tz')
        disp('Finished exporting tz-curves')
    end
    