function [database1, T1] = read_copcat_database()

%% Import COPCAT_Data_Base Table
db_server              ='DKLYCOPILOD1';     % Databse server
db_user                ='owdb_user';        % Database user
db_pass                ='ituotdowdb';       % Database pass
db_name                ='owdb';             % Database name  

            
tablename = 'COPCAT_Data_Base';
% column_name='name'; % query_string='Batch3_N1';

mysql('open',db_server,db_user,db_pass);
mysql(['use ', db_name]);

mysql_command_string = ['SELECT COUNT(*) AS Columns FROM INFORMATION_SCHEMA.COLUMNS WHERE table_schema = ', ['''',db_name,''''],' AND table_name = ',['''', tablename, '''']];
NumberColumns = mysql(mysql_command_string);


mysql_command_string = ['SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE table_schema = ', ['''',db_name,''''],' AND table_name = ',['''', tablename, '''']];
column_names = mysql(mysql_command_string);


%mysql_command_string=['SELECT * FROM ',tablename,' WHERE ', column_name ,' = ',['''',query_string,'''']];
mysql_command_string=['SELECT * FROM ',tablename];
% mysql(mysql_command_string)


T1=table(); % T1 will contain the desire tablename
vars_str = '[';
for ii=1:NumberColumns  
    vars_ii=['T1.',column_names{ii},','];
    vars_str=[vars_str,vars_ii]; 
end
vars_str(end)=']';


eval([vars_str,' = mysql(mysql_command_string);']);

mysql('close')
disp('COPCAT database was imported successfully!')

%% Create a unique row name for each row entry
% the values of the folowing columns will be concatenated:
% row_name = project___location___name___rev, there are three underscores as separater
%  e.g. row_name = 'Oceanwind___OSSC___U90S_9___0'

Project = T1.project;
Location = strcat('___',T1.location);
Name = strcat('___',T1.name);


Rev_string = cellstr(num2str([T1.rev])); % convert rev to string
Rev_string = replace(Rev_string,' ',''); % remove the spaces
Rev = strcat('___',Rev_string);

row_name = strcat(Project,Location,Name,Rev);

T1.row_name = row_name;
T1 = movevars(T1,'row_name','Before','name');

%% Check for duplicate rows
[unique_row_id,ia,~] = unique(row_name);
idx = 1:length(row_name); idx(ia)=[];

if length(row_name)-length(unique_row_id)>0
    disp('There are duplicate entries in the copcat database. Check:');
    disp(char(row_name(idx)));
    error('Remove duplicate entries in the copcat database!');    
end

%% 

params = {'y_u_F','y_u_1','y_u_2','y_u_3','P_u_F','P_u_1','P_u_2','P_u_3',...
    'k_p_F','k_p_1','k_p_2','k_p_3','n_p_F','n_p_1','n_p_2','n_p_3',...
    'tetam_u_F', 'tetam_u_1','tetam_u_2','tetam_u_3','m_u_F','m_u_1','m_u_2','m_u_3',...
    'k_m_F','k_m_1','k_m_2','k_m_3','n_m_F','n_m_1','n_m_2','n_m_3','yB_u_F',...
    'yB_u_1','yB_u_2','yB_u_3','HB_u_F','HB_u_1','HB_u_2','HB_u_3','k_H_F',...
    'k_H_1','k_H_2','k_H_3','n_H_F','n_H_1','n_H_2','n_H_3','tetaMb_u_F',...
    'tetaMb_u_1','tetaMb_u_2','tetaMb_u_3','MB_u_F','MB_u_1','MB_u_2',...
    'MB_u_3','k_Mb_F','k_Mb_1','k_Mb_2','k_Mb_3','n_Mb_F','n_Mb_1',...
    'n_Mb_2','n_Mb_3'};


database1.num = T1{:,params};
database1.txt = [T1.row_name, T1.project, T1.notes, T1.soil_type];

T2 = T1(:,['row_name','project','notes','soil_type',params]);

Headers = T2.Properties.VariableNames;
database1.raw = [Headers; Headers; Headers ; table2cell(T2)];


