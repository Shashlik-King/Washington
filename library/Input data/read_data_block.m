function [data_table]=read_data_block(data_block_header,raw0)

disp(['Reading datablock : ', data_block_header])

%data_block_header='Run settings';
cond1=strcmp(raw0,data_block_header);
if sum(cond1(:))>1
    error_string = ['Remove cell entries identical to the data block header : ' data_block_header];
    error(error_string);
elseif sum(cond1(:))==0 
    error_string = ['Data block header could not be found : ' data_block_header];
    error(error_string);
end


[row0,col0] = find(cond1);
Nrows=raw0{row0,col0+3};
Ncols=raw0{row0,col0+4};
if Nrows<2 || Ncols<2   
    error_string = ['Check Data block dimensions. Datablock should be at least 1 row and 1 column : ' data_block_header];
    error(error_string);    
end

rows=(row0:row0-1+Nrows)+1;
cols=(col0:col0-1+Ncols)+0;
data_block=raw0(rows,cols);

cond1=ismissing(string(data_block)); % rows containing nan
cond2=strcmp(string(data_block),''); % empty rows
cond3 = sum(cond1 | cond2,2) ==0;    % rows without nans or empty cells
data_block = data_block(cond3,:);    % filtered data block

data_table = cell2table(data_block(2:end,:),'VariableNames',data_block(1,:));




