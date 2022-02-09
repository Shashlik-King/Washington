function database_write_toesprings_multi(settings,element,loads,data,pile,soil,plots,output,i,p_toe,y_toe,m_toe,teta_toe,database,scour)
%--------------------------------------------------------------------------
% CHANGE LOG
% 2020.03.19    FKMV    Programming

L_increments = pile.toe_max+pile.toe_min;
if strcmp(settings.database,'Manual')
    pile.length_start_or = pile.length_start;
end
disp(['Exporting toe springs: from L = ',num2str(pile.length_start_or-pile.toe_min),...
    ' to L = ',num2str(pile.length_start_or+pile.toe_max)]),

settings.multi_toe_levels = settings.multi_toe_levels + 1; % change of switch to appy different embedded depth
%% Initial upload of first depth

for ii = 1:L_increments % 1m increments are assumed
    clear p_toe y_toe m_toe teta_toe
    pile.multi_length_start_or = pile.length_start_or + pile.toe_max - (ii-1); % reduce pile length by 1m per iteration
    % re-run of calculations for different pile depth
%     [scour,soil,pile,loads,data,settings] = database_input(pile,soil,data,...
%     loads,settings,plots); % load soil and pile information from database (SNA module)
    pile.L = pile.multi_length_start_or; % Selecting pile length
    disp(['Writing toe springs for pile length: ',num2str(pile.L),' m'])
    pile.toe = pile.head - pile.L; % pile toe level
    pile.length=pile.multi_length_start_or;
    [element,node,pile] = elem_node_data(pile,soil,scour,settings,data,plots,1,database); % assign properties to elements and nodes
    [element] = layer(element,node,pile,settings,loads);
    [p_toe,y_toe] = p_y_toe(node,pile,element,loads,output,settings);
    [m_toe,teta_toe] = m_teta_toe(node,pile,element,loads,output,settings);
    database_write_toesprings(settings,element,loads,data,p_toe,y_toe,...
    strcat('soil_py_toe_curves_',data.table_springs),'py'); % write py-curves into database
    database_write_toesprings(settings,element,loads,data,m_toe,teta_toe,...
    strcat('soil_Mt_toe_curves_',data.table_springs),'Mt'); % write Mt-curves into database

disp('Finished exporting toe springs')
end