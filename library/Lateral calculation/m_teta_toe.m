function [m teta] = m_teta_toe(node,pile,element,loads,output,settings)
%% FORMULATION FOR p-y CURVES
%--------------------------------------------------------------------------
% CHANGE LOG
% 01.08.2013    MUOE - SEPARATING FROM WINKLER_PY.M
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% pile.xx:      pile.diameter is used to determine yc
% element.xx:   model.py is used to determine which p-y curve to use
%               epsilon50 is used to determine yc
%--------------------------------------------------------------------------
%% Output parameters
%--------------------------------------------------------------------------
% p:        lateral resistance
% y:        lateral displacement
%--------------------------------------------------------------------------
%% Initial
%--------------------------------------------------------------------------
nelem = length(element.model_py);
y_tot = zeros(nelem,2);
%--------------------------------------------------------------------------
%% Calculation routine
%--------------------------------------------------------------------------
c = nelem;
if element.PISA_switch==1
    teta_top(c,:) = [0 0.0000003 0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.01 0.1];
    teta_bot(c,:) = teta_top(c,:);
end

npoints = size(teta_bot,2);

for j = 1:npoints
    [ksmtetanode_top ksmtetanode_bot] = secMBstiff(element,pile,[teta_top(c,j) teta_bot(c,j)],c);
    mmtetanode_top(c,j) = ksmtetanode_top*teta_top(c,j);  
    mmtetanode_bot(c,j) = ksmtetanode_bot*teta_bot(c,j);  
end


% Putting p and y data into same format as t and z
for i = 1:size(teta_top,1)
    m.top(i,:)      = mmtetanode_top(i,:);
    m.bottom(i,:)   = mmtetanode_bot(i,:);
    teta.top(i,:)      = teta_top(i,:);
    teta.bottom(i,:)   = teta_bot(i,:);
end