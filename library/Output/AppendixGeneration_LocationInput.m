function AppendixGeneration_LocationInput(Settings,AppendixData)
% Generate the .tex input for each location

%% Blow count per meter
% % [value1, value2] = max(SRD.(Position).SOD(:,5)); 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxBCm_BE.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1);
% % fclose(FID);
% % 
% % 
FID = fopen(char('AppendixGenerationFiles/ProjectLocation/ID_location.tex'), 'w+');
fprintf(FID, '%s', AppendixData.locationID);
fclose(FID);

% FID = fopen(char('AppendixGenerationFiles/ProjectLocation/WaterDepth.tex'), 'w+');
% fprintf(FID, '%s', results.ID{c}.def.pile.toe-results.ID{c}.def.pile.length);
% fclose(FID);

% FID = fopen(char('AppendixGenerationFiles/ProjectLocation/EmbDepth.tex'), 'w+');
% fprintf(FID, '%s', embed_length);
% fclose(FID);
% 
% FID = fopen(char('AppendixGenerationFiles/ProjectLocation/BSHcheck.tex'), 'w+');
% fprintf(FID, '%s', results.ID{c}.def.pile.length);
% fclose(FID);
% 
% FID = fopen(char('AppendixGenerationFiles/ProjectLocation/DNVGLcheck.tex'), 'w+');
% fprintf(FID, '%s', results.ID{c}.def.pile.length);
% fclose(FID);

if Settings.Appendix.Axial_pile_capacity
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/AxialCapacity.tex'), 'w+');
    fprintf(FID, '%s', AppendixData.results.axial_cap(1,6));
    fclose(FID);
end

if Settings.Appendix.GeneralInfo
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/pile_depth.tex'), 'w+');
    fprintf(FID, '%1.2f', AppendixData.pileLength);
    fclose(FID);
end

if Settings.Appendix.GeneralInfo
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/pile_length.tex'), 'w+');
    fprintf(FID, '%1.2f', abs(AppendixData.pileTip)+abs(AppendixData.pileTop));
    fclose(FID);
end

if Settings.Appendix.GeneralInfo
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/MP_top.tex'), 'w+');
    fprintf(FID, '%1.2f', AppendixData.pileTop);
    fclose(FID);
end

if Settings.Appendix.GeneralInfo
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/MP_bottom.tex'), 'w+');
    fprintf(FID, '%1.2f', AppendixData.pileTip);
    fclose(FID);
end

if Settings.Appendix.GeneralInfo
    FID = fopen(char('AppendixGenerationFiles/ProjectLocation/water_depth.tex'), 'w+');
    fprintf(FID, '%1.2f', abs(AppendixData.pileTip + AppendixData.pileLength));
    fclose(FID);
end

% % 
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxBCmdepth_BE.tex'), 'w+');
% % fprintf(FID, '%1.1f', SRD.(Position).SOD(value2,1));
% % fclose(FID);
% % 
% % [value1, value2] = max(SRD_UB.(Position).SOD(:,5));
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxBCm_UB.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1);
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxBCmdepth_UB.tex'), 'w+');
% % fprintf(FID, '%1.1f', SRD_UB.(Position).SOD(value2,1));
% % fclose(FID);

%% Accumulated blow     -   Check this 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/AccumulatedBC_BE.tex'), 'w+');
% % fprintf(FID, '%1.0f', max(cumsum(SRD.(Position).SOD(:,10))));
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/AccumulatedBC_UB.tex'), 'w+');
% % fprintf(FID, '%1.0f', max(cumsum(SRD_UB.(Position).SOD(:,10))));
% % fclose(FID);

%% SRD
% % [value1, value2] = max(SRD.(Position).SOD(:,2));
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxSRD_BE.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1/1000);
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxSRDdepth_BE.tex'), 'w+');
% % fprintf(FID, '%1.1f', SRD.(Position).SOD(value2,1));
% % fclose(FID);
% % 
% % 
% % [value1, value2] = max(SRD_UB.(Position).SOD(:,2));
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxSRD_UB.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1/1000);
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxSRDdepth_UB.tex'), 'w+');
% % fprintf(FID, '%1.1f', SRD_UB.(Position).SOD(value2,1));
% % fclose(FID);

%% Fatigue damage
% % [value1, value2] = max(SRD.(Position).D);
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxfatique_BE.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1*100);
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxfatiquedepth_BE.tex'), 'w+');
% % fprintf(FID, '%1.2f', Data.(Position).SCF{value2,1});
% % fclose(FID);
% % 
% % [value1, value2] = max(SRD_UB.(Position).D);
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxfatique_UB.tex'), 'w+');
% % fprintf(FID, '%1.1f', value1*100);
% % fclose(FID);
% % 
% % FID = fopen(char('AppendixGenerationFiles/ProjectLocation/maxfatiquedepth_UB.tex'), 'w+');
% % fprintf(FID, '%1.2f', Data.(Position).SCF{value2,1});
% % fclose(FID);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pile table for appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry(canID, top diameter, bottom diamter,  can height, thickness)

% geometry=[Data.(Position).PileGeometry(:,1) Data.(Position).PileGeometry(:,2) Data.(Position).PileGeometry(:,3) Data.(Position).PileGeometry(:,4)]; 

% FKMV
% for i = 1:size(PileGeometry,1)
%     if i == 1
%        pile_geometry(i,:) = [i, can.top, geometry(i,1), geometry(i,2), geometry(i,4), geometry(i,3), geometry(i,6), geometry(i,5)];
%     elseif i == size(geometry,1)
%         pile_geometry(i,:) = [i, cell2mat(pile_geometry(i-1,2))-cell2mat(pile_geometry(i-1,5)), geometry(i,1), geometry(i,2), geometry(i,4), geometry(i,3), NaN, NaN];
%     else
%        pile_geometry(i,:) = [i, cell2mat(pile_geometry(i-1,2))-cell2mat(pile_geometry(i-1,5)), geometry(i,1), geometry(i,2), geometry(i,4), geometry(i,3), geometry(i,6), geometry(i,5)];  
%     end
% end
if Settings.Appendix.Pile_Geometry
    F = fopen(char('AppendixGenerationFiles/ProjectLocation/pile_geometry.tex'),'wt');
    PileGeometry = AppendixData.geometryTable;
    for i = 1:size(PileGeometry,1)
       if i ==  size(PileGeometry,1)
           fprintf(F,'%.0f & %.3f & %.3f & %.3f & %.3f  \\\\\\hline \n',cell2mat(PileGeometry(i,1)),cell2mat(PileGeometry(i,2)),cell2mat(PileGeometry(i,3)),cell2mat(PileGeometry(i,4)),cell2mat(PileGeometry(i,5)));
       else
           fprintf(F,'%.0f & %.3f & %.3f & %.3f & %.3f  \\\\\\hline \n',cell2mat(PileGeometry(i,1)),cell2mat(PileGeometry(i,2)),cell2mat(PileGeometry(i,3)),cell2mat(PileGeometry(i,4)),cell2mat(PileGeometry(i,5))); 
       end
    end
    fclose(F);
end

%% SSI springs table
% p-y
if Settings.Appendix.SSI_springs
    F = fopen(char('AppendixGenerationFiles/ProjectLocation/SSI_p_y.tex'),'wt');
    SSI = AppendixData.SSI_py;
    for i = 1:size(SSI,1)/2
           fprintf(F,'%.2f & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\\hline \n', AppendixData.SSI_elev(i,1),'p', SSI(i,1), SSI(i,2), SSI(i,3), SSI(i,4), SSI(i,5), SSI(i,6), SSI(i,7), SSI(i,8), SSI(i,9), SSI(i,10), SSI(i,11), SSI(i,12), SSI(i,13), SSI(i,14), SSI(i,15), SSI(i,16), SSI(i,17) );
           fprintf(F,'%.2f & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\\hline \n', AppendixData.SSI_elev(i,1),'y', SSI(end/2+i,1), SSI(end/2+i,2), SSI(end/2+i,3), SSI(end/2+i,4), SSI(end/2+i,5), SSI(end/2+i,6), SSI(end/2+i,7), SSI(end/2+i,8), SSI(end/2+i,9), SSI(end/2+i,10), SSI(end/2+i,11), SSI(end/2+i,12), SSI(end/2+i,13), SSI(end/2+i,14), SSI(end/2+i,15), SSI(end/2+i,16), SSI(end/2+i,17) );
    end
    fclose(F);
end

% m-t
if Settings.Appendix.SSI_springs
    F = fopen(char('AppendixGenerationFiles/ProjectLocation/SSI_m_t.tex'),'wt');
    SSI = AppendixData.SSI_mt;
    for i = 1:size(SSI,1)/2
           fprintf(F,'%.2f & %s & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f  \\\\\\hline \n', AppendixData.SSI_elev(i,1),'m', SSI(i,1), SSI(i,2), SSI(i,3), SSI(i,4), SSI(i,5), SSI(i,6), SSI(i,7), SSI(i,8), SSI(i,9), SSI(i,10), SSI(i,11));
           fprintf(F,'%.2f & %s & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f  \\\\\\hline \n', AppendixData.SSI_elev(i,1),'theta', SSI(end/2+i,1), SSI(end/2+i,2), SSI(end/2+i,3), SSI(end/2+i,4), SSI(end/2+i,5), SSI(end/2+i,6), SSI(end/2+i,7), SSI(end/2+i,8), SSI(end/2+i,9), SSI(end/2+i,10), SSI(end/2+i,11));
    end
    fclose(F);
end

% p-y toe
if Settings.Appendix.SSI_springs
    F = fopen(char('AppendixGenerationFiles/ProjectLocation/SSI_p_y_base.tex'),'wt');
    SSI = AppendixData.SSI_py_base;
%     for i = 1:size(SSI,1)
           fprintf(F,'%.2f & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\\hline \n', AppendixData.SSI_elev(end,1),'p', SSI(1,1), SSI(1,2), SSI(1,3), SSI(1,4), SSI(1,5), SSI(1,6), SSI(1,7), SSI(1,8), SSI(1,9),SSI(1,10), SSI(1,11), SSI(1,12), SSI(1,13), SSI(1,14), SSI(1,15), SSI(1,16), SSI(1,17) );
           fprintf(F,'%.2f & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\\hline \n', AppendixData.SSI_elev(end,1),'y', SSI(2,1), SSI(2,2), SSI(2,3), SSI(2,4), SSI(2,5), SSI(2,6), SSI(2,7), SSI(2,8), SSI(2,9), SSI(2,10), SSI(2,11), SSI(2,12), SSI(2,13), SSI(2,14), SSI(2,15), SSI(2,16),SSI(2,17) );
           %     end
    fclose(F);
end

% m-t toe
if Settings.Appendix.SSI_springs
    F = fopen(char('AppendixGenerationFiles/ProjectLocation/SSI_m_t_base.tex'),'wt');
    SSI = AppendixData.SSI_mt_base;
%     for i = 1:size(SSI,1)
           fprintf(F,'%.2f & %s & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f  \\\\\\hline \n', AppendixData.SSI_elev(end,1),'m', SSI(1,1), SSI(1,2), SSI(1,3), SSI(1,4), SSI(1,5), SSI(1,6), SSI(1,7), SSI(1,8), SSI(1,9), SSI(1,10), SSI(1,11));
           fprintf(F,'%.2f & %s & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f & %.8f  \\\\\\hline \n', AppendixData.SSI_elev(end,1),'theta', SSI(2,1), SSI(2,2), SSI(2,3), SSI(2,4), SSI(2,5), SSI(2,6), SSI(2,7), SSI(2,8), SSI(2,9), SSI(2,10), SSI(2,11));
           %     end
    fclose(F);
end

%% Create footer file
%char(strcat(Project,'/',Position,'/Plots/',cell2mat(Position),' footer.tex')), 'w+')
% FID = fopen(char('AppendixGenerationFiles/ProjectLocation/footer.tex'), 'w+');
% fprintf(FID, '\\Appendix\\_for\\_%s', Settings.Appendix.DocumentNoCOWI,Position,Settings.Appendix.RevisionTable{end,1});
% fclose(FID);

