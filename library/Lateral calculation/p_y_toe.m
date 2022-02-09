function [p y] = p_y_toe(node,pile,element,loads,output,settings)
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
    upy_top(c,:) = [0 0.0001 0.0005 0.001 0.003 0.005 0.01 0.0125 0.015 0.02 0.03 0.04 0.06 0.080 0.12 0.2 0.4];
    upy_bot(c,:) = upy_top(c,:);
end

npoints = size(upy_top,2);

for j = 1:npoints
    [ksppynode_top ksppynode_bot] = secHBstiff(element,pile,[upy_bot(c,j) upy_bot(c,j)],c);
    ppynode_top(c,j) = ksppynode_top*upy_top(c,j);  
    ppynode_bot(c,j) = ksppynode_bot*upy_bot(c,j);  
end


% Putting p and y data into same format as t and z
for i = 1:size(upy_top,1)
    p.top(i,:)      = ppynode_top(i,:);
    p.bottom(i,:)   = ppynode_bot(i,:);
    y.top(i,:)      = upy_top(i,:);
    y.bottom(i,:)   = upy_bot(i,:);
end