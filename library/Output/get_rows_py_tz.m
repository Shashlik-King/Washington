function [row_py_tz,row_unique] = get_rows_py_tz(element,scour,n_curves)
%% FUNCTION TO GET THE CORRECT NUMBER OF ROWS TO PLOT P-Y AND T-Z CURVES
% MINIMUM 1 CURVE PR. LAYER, MINIMUM n_curves CURVES IN TOTAL
%--------------------------------------------------------------------------
% CHANGE LOG
% 07.11.2013    MUOE - PROGRAMMING
% 07.01.2014    MUOE - SCOUR IS TAKEN INTO ACCOUNT
%--------------------------------------------------------------------------
depths_unique = unique(element.level(:,3));
pos_under_scour = find(depths_unique < -scour.local);
depths_unique_scour = depths_unique(pos_under_scour);
depths_unique_scour(end+1) = -scour.local;
for i = 1:length(depths_unique_scour)
    if i < length(depths_unique_scour)
        row_unique(i) = find(element.level(:,3) == depths_unique_scour(i),1);
    else
        row_unique(i) = find(element.level(:,1) == depths_unique_scour(i));
    end
end
row_unique = sort(row_unique,'ascend');
row_choose = row_unique;
row_choose(end+1) = length(element.level);
while length(row_choose) <= n_curves
    for j = 1:length(row_choose)-1
        diff(j) = row_choose(j+1)-row_choose(j);
    end
    [val,pos] = max(diff);
    row_choose(end+1) = floor((row_choose(pos)+row_choose(pos+1))/2);
    row_choose = sort(row_choose,'ascend');
end
for k = 1:length(row_choose)-1
    row_py_tz(k) = floor((row_choose(k)+row_choose(k+1))/2);
end