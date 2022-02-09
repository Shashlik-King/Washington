function GeniE_shallow_pytz(pile,soil,data,element,scour,settings,t,z,p,y,Q,zQ,plug_unplug,iii,loads,fs,levels)
%% FILE FOR GENERATING INPUT FOR Sesam GeniE
%--------------------------------------------------------------------------
% CHANGE LOG
% 2019.06.25        EMSS        Programming
% 2019.07.22        EVVA        Programming
%--------------------------------------------------------------------------
%% Basic file information - message to Structural designer
%--------------------------------------------------------------------------
[~, ~, ~]=mkdir(settings.FolderSaveInputToSTR);

switch loads.static_cyclic
	case 'cyclic'
		scenario = 'BE';
	case 'static'
		scenario = 'UB';
end

name = ['OD',mat2str(pile.diameter),'m_', data.project(1:3),'_','id-REP','_',scenario, '.txt'];
F = fopen([settings.FolderSaveInputToSTR,'\',name],'wt');

%Loc = 'DEEP';
Loc =settings.GeniE_Location;

nSoils=length(soil.toplevel);
%F = fopen(name,'Wt');
fprintf(F, '//********************************************************************************\n');
fprintf(F, '// GEO input in the form of p-y, t-z and q-z curves\n');
fprintf(F, '// File created by:%s, %s \n',data.prepared_by, datestr(now));
fprintf(F, '//********************************************************************************\n');
fprintf(F, '// Location_%s - %s - %s soil profile\n', Loc, data.project, loads.static_cyclic);
fprintf(F, '//********************************************************************************\n');
fprintf(F, '// Scour is accounted for in this input, do not apply additionally in GeniE\n');
fprintf(F, '//********************************************************************************\n');
fprintf(F, '// Water depth: zMudline,  input in GeniE\n');
fprintf(F, '//********************************************************************************\n');
fprintf(F, '// Pile embedment length %.2f\n',pile.length(iii));
fprintf(F, '// The last soil curve set is additionally copied at 0.5 m below pile embedment\n');
fprintf(F, '// length to ensure that pile end bearing data is applied in SESAM');
fprintf(F, '\n//********************************************************************************\n\n\n');

%--------------------------------------------------------------------------
%% Soil curves
%--------------------------------------------------------------------------
% Generate string of y-values. Since they are identical for all p-y
% curves, we only need to do it once
yString=strjoin(compose('%11.8E m', [y.top(end,:), y.top(end,end)*2]), ', '); % Create string of all y-values with m unit, seperated by comma ','

% Generate string of z-values. Since they are identical for all t-z
% curves, we only need to do it once
zString=strjoin(compose('%11.8E m', [z.top(end,:), z.top(end,end)*2]), ', '); % Create string of all z-values with m unit, seperated by comma ','

% Calculate tip area
if strcmp(plug_unplug.comp_index,'Unplugged')
    Q_area = pi*(pile.diameter^2-(pile.diameter-2*element.thickness(end))^2)/4; %GMME
elseif strcmp(plug_unplug.comp_index,'Plugged')
    Q_area = pi*pile.diameter^2/4; %GMME
else
    error('Pile neither plugged nor unplugged')
end

soilId=0;
% Write custom curves to file
nCurves=size(p,1);
for i=1:nCurves % Loop over all curves 
    % Find relevant soil layer
    while (soilId == 0) || ((soilId < nSoils) && (soil.toplevel(soilId+1) > levels(i))) % If we have not reached last soil layer then check if toplevel of next soil layer is above current level, if so, switch soil layer.
        soilId=soilId+1;
        curveId=1;
        fprintf(F,'//********************************************************************************\n');
        fprintf(F,'// Soil curves for layer %d\n', soilId);
        fprintf(F,'//********************************************************************************\n');
    end
    
     % Generate string of p-values
%   pString=strjoin(compose('%11.8E kPa', [p(i,:), p(i, end)]/(pi*pile.diameter)), ', '); % Create string of all p-values with kPa unit, seperated by comma ','
   pString=strjoin(compose('%11.8E kPa', [p(i,:), p(i, end)]/pile.diameter), ', '); % Create string of all p-values with kPa unit, seperated by comma ','
   % Generate string of t-values
    tString=strjoin(compose('%11.8E kPa', [t(i,:), t(i,end)]), ', '); % Create string of all t-values with kPa unit, seperated by comma ','
    
    
   % Generate strings concerning tip resistance
% THIS IS Q IN [kN]
%     if i~=nCurves-1 % Zero stiffness for anything but the final layer with the pile tip %GMME 
%         QString='0 kN, 0 kN';
%         zQString='0 m, 1 m';
%         peakTipString='0 kN'; 
%     else
%         QString=strjoin(compose('%11.8E kN', [Q.bottom(end,:), Q.bottom(end,end)]*Q_area), ', ');
%         zQString=strjoin(compose('%11.8E m', [zQ.bottom(end,:), zQ.bottom(end,end)*2]), ', ');
%         peakTipString=sprintf('%11.8E kN;', Q.bottom(end,end)*Q_area);
%     end
%  
%     
%     fprintf(F, 'Layer%d_SoilCurves%d = SoilCurves(pyManual, tzManual, qzManual);\n', soilId, curveId);
%     fprintf(F, 'Layer%d_SoilCurves%d.addManualPY(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, yString, pString);
%     fprintf(F, 'Layer%d_SoilCurves%d.addManualTZ(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, zString, tString);
%     fprintf(F, 'Layer%d_SoilCurves%d.addManualQZ(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, zQString, QString);
% %   fprintf(F, 'Layer%d_SoilCurves%d.addAdjustPY(0, %11.8E, %11.8E, 0, 1, 1, 1);\n', soilId, curveId, soil.degradation.value_py_p(soilId), soil.degradation.value_py_y(soilId));
% %   fprintf(F, 'Layer%d_SoilCurves%d.addAdjustTZ(0, %11.8E, %11.8E, 1, 1, %11.8E);\n\n', soilId, curveId, soil.degradation.value_tz_t(soilId), soil.degradation.value_tz_t(soilId), soil.degradation.value_tz_z(soilId));
%     fprintf(F, 'Layer%d_SoilCurves%d.addAdjustPY(0, 1, 1, 0, 1, 1, 1);\n', soilId, curveId);
%     fprintf(F, 'Layer%d_SoilCurves%d.addAdjustTZ(0, 1, 1, 1, 1, 1);\n\n', soilId, curveId);
    
% THIS IS Q IN [kPa]
    if i<nCurves-2 % Zero stiffness for all nodes except the last three that are input in the model 
        QString='0 kPa, 0 kPa';
        zQString='0 m, 1 m';
        peakTipString='0 kPa'; 
    else  %The last three elements have the same q-z values that are in the pile tip vicinity. 
        QString=strjoin(compose('%11.8E kPa', [Q.bottom(end,:), Q.bottom(end,end)]), ', ');
        zQString=strjoin(compose('%11.8E m', [zQ.bottom(end,:), zQ.bottom(end,end)*2]), ', ');
        peakTipString=sprintf('%11.8E kPa;', Q.bottom(end,end));   
    end
 
    
    fprintf(F, 'Layer%d_SoilCurves%d = SoilCurves(pyManual, tzManual, qzManual);\n', soilId, curveId);
    fprintf(F, 'Layer%d_SoilCurves%d.addManualPY(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, yString, pString);
    fprintf(F, 'Layer%d_SoilCurves%d.addManualTZ(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, zString, tString);
    fprintf(F, 'Layer%d_SoilCurves%d.addManualQZ(%11.8E, Array(%s), Array(%s));\n', soilId, curveId, pile.diameter, zQString, QString);
% p, t multipliers are already included in the curves, so they are set to 1
% as input to GeniE. 
%   fprintf(F, 'Layer%d_SoilCurves%d.addAdjustPY(0, %11.8E, %11.8E, 0, 1, 1, 1);\n', soilId, curveId, soil.degradation.value_py_p(soilId), soil.degradation.value_py_y(soilId));
%   fprintf(F, 'Layer%d_SoilCurves%d.addAdjustTZ(0, %11.8E, %11.8E, 1, 1, %11.8E);\n\n', soilId, curveId, soil.degradation.value_tz_t(soilId), soil.degradation.value_tz_t(soilId), soil.degradation.value_tz_z(soilId));
    fprintf(F, 'Layer%d_SoilCurves%d.addAdjustPY(0, 1, 1, 0, 1, 1, 1);\n', soilId, curveId);
    fprintf(F, 'Layer%d_SoilCurves%d.addAdjustTZ(0, 1, 1, 1, 1, 1);\n\n', soilId, curveId);
    
    % Increment curveId
    curveId=curveId+1;

end

%--------------------------------------------------------------------------
%% Soil data and coefficients required as a template input to GeniE. These are not applied to soil calculations. 
%--------------------------------------------------------------------------
fprintf(F,'//********************************************************************************\n');
fprintf(F,'// Inserting Material Coefficients and Dummy with Soil Data\n');
fprintf(F,'//********************************************************************************\n');
fprintf(F,'\n');
fprintf(F,'Tan_Phi_Coeff = 1.0;\n'); 
fprintf(F,'Shear_Strength_Coeff = 1.0;\n'); 
fprintf(F,'Skin_Friction_Coeff = 1.0;\n'); 
fprintf(F,'Pile_Tip_Resistance_Coeff = 1.0;\n'); 
fprintf(F,'TZ_curve_Shape_Factor = 1.00;\n'); 
fprintf(F,'LowestShearStiff = 1.0 MPa;\n'); 
fprintf(F,'\n');

%--------------------------------------------------------------------------
%% Pile group effect
%--------------------------------------------------------------------------
% fprintf(F,'// For pile group effects \n'); 
% fprintf(F,'\n');
% fprintf(F,'E_MudL = 5 MPa;\n'); 
% fprintf(F,'E_MudL_1m = 14.0 MPa;\n'); 
% fprintf(F,'TZ_zone_Of_Influence_Ratio = 5;\n'); 
% fprintf(F,'\n');


%--------------------------------------------------------------------------
%% Dummy soil
%--------------------------------------------------------------------------
fprintf(F,'//********************************************************************************\n');
fprintf(F,'// Dummy soil \n'); 
fprintf(F,'//********************************************************************************\n');
fprintf(F,'\n');
fprintf(F,'Layer_Dummy = Sand(false, 1, 1.961 tonne/m^3, 27.5 deg, 5400 kN/m^3);\n'); 
fprintf(F,'SoilData_Dummy  = SoilData(0 kPa, 0.3, 72.6 kPa, 72.6 kPa, 0.01, 3660 kPa, 0.1);\n\n'); 

%--------------------------------------------------------------------------
%% Defining the Locations of the Soil Layers
%--------------------------------------------------------------------------
%% Layers for LWL, MSL and HWL (low water level, mean water level and high water level)
fprintf(F,'//********************************************************************************\n');
fprintf(F,'// Defining the Locations of the Soil Layers\n');
fprintf(F,'//********************************************************************************\n');
fprintf(F,'\n');

Locs = {'MSL'};

for j = 1:size(Locs)
    Loc = Locs{j};
    fprintf(F, 'Location_%s.relativeSoilLayers = false;\n', Loc);
    fprintf(F, 'Location_%s.clearSoilLayers();\n', Loc);
    soilId=0;
    
    for i=1:(nCurves)-1 % Loop over all curves and create a sublayer for each curve
        % Find relevant soil layer
        while (soilId == 0) || ((soilId < nSoils) && (soil.toplevel(soilId+1) > levels(i))) % If we have not reached last soil layer then check if toplevel of next soil layer is above current level, if so, switch soil layer.
            soilId=soilId+1;
            curveId=1;
        end
        fprintf(F, 'Location_%s.insertSoilBorder(%11.8E m +zMudline);\n', Loc, levels(i+1));
        fprintf(F, 'Location_%s.soil(%d).soilCurves = Layer%d_SoilCurves%d;\n', Loc, i, soilId, curveId);
        fprintf(F, 'Location_%s.soil(%d).numberOfSublayers = 1;\n', Loc, i);
        fprintf(F, 'Location_%s.soil(%d).soilType = Layer_Dummy;\n', Loc, i);
        fprintf(F, 'Location_%s.soil(%d).soilData = SoilData_Dummy;\n', Loc, i);  

        % Increment curveId
        curveId=curveId+1;
    end

    % Create a copy of the last soil layer 0.5 m below pile tip to ensure
    % that SESAM assigns pile end bearing higher than zero
    fprintf(F, 'Location_%s.insertSoilBorder(%11.8E m +zMudline);\n', Loc, levels(nCurves)-0.5);
    fprintf(F, 'Location_%s.soil(%d).soilCurves = Layer%d_SoilCurves%d;\n', Loc, nCurves, soilId, curveId-1);
    fprintf(F, 'Location_%s.soil(%d).numberOfSublayers = 1;\n', Loc, nCurves);
    fprintf(F, 'Location_%s.soil(%d).soilType = Layer_Dummy;\n', Loc, nCurves);
    fprintf(F, 'Location_%s.soil(%d).soilData = SoilData_Dummy;\n', Loc, nCurves);  
    
    fprintf(F,'//********************************************************************************\n\n');

end

fclose(F);