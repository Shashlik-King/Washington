function AppendixGeneration_SiteInput(Settings)
% Appendix - Generation of files with input to documentation
%% Create project specific files for the appendices
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Project_name.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.ProjectName);
    fclose(FID);
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Project_number.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.ProjectNumber);
    fclose(FID);
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Document_number.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.DocumentNoCOWI);
    fclose(FID);
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Document_number_client.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.DocumentNoClient);
    fclose(FID);
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Project_client.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.ProjectClient);
    fclose(FID);
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Document_date.tex'}), 'w+');
    fprintf(FID, '%s', Settings.Appendix.DocumentDate);
    fclose(FID);
    
%     FID = fopen(char({'AppendixGenerationFiles\ProjectSite\Hammer.tex'}), 'w+');
%     fprintf(FID, '%s', 'IHC TEST HAMMER - MAKE REAL INPUT HERE IN CODE');
%     fclose(FID);
    
%% Generate front page table automatically
    value_position = round((size(Settings.Appendix.RevisionTable,1)-2)*0.4,2); % Index for changing layout for front page by moving tables - can be adjusted
    
    FID = fopen(char({'AppendixGenerationFiles\ProjectSite\FrontpageTable.tex'}), 'w+');
    fprintf(FID, '\\put( -4, %1.1f){\\begin{minipage}[b]{18cm}\\scriptsize \n',value_position);
    fprintf(FID, '\\begin{tabular}{llll} \n');
%     if length(Settings.Appendix.DocumentNoClient) > 1 && length(Settings.Appendix.DocumentNoEmployer) > 1
%         fprintf(FID, '\\textcolor{COWI}{\\tiny{PROJECT NO.}} & \\textcolor{COWI}{\\tiny{DOCUMENT NO.}} & \\textcolor{COWI}{\\tiny{%s DOCUMENT NO.}} & \\textcolor{COWI}{\\tiny{%s DOCUMENT NO.}} \\\\ \n',Settings.Appendix.ProjectClient,Settings.Appendix.ProjectEmployer);
%         fprintf(FID, '&  &  &  \\\\ \n');
%         fprintf(FID, '%s & %s & %s & %s \\\\ \n',Settings.Appendix.ProjectNumber,Settings.Appendix.DocumentNoCOWI,Settings.Appendix.DocumentNoClient,Settings.Appendix.DocumentNoEmployer);
%     elseif length(Settings.Appendix.DocumentNoClient) > 1   % Only COWI and client
%         fprintf(FID, '\\textcolor{COWI}{\\tiny{PROJECT NO.}} & \\textcolor{COWI}{\\tiny{DOCUMENT NO.}} & \\multicolumn{2}{l}{\\textcolor{COWI}{\\tiny{%s DOCUMENT NO.}}} \\\\ \n',Settings.Appendix.ProjectClient);
%         fprintf(FID, '&  &  & \\\\ \n');
%         fprintf(FID, '%s & %s & \\multicolumn{2}{l}{%s} \\\\ \n',Settings.Appendix.ProjectNumber,Settings.Appendix.DocumentNoCOWI,Settings.Appendix.DocumentNoClient);
%     elseif length(Settings.Appendix.DocumentNoEmployer) > 1 % Only COWI and eployer
%         fprintf(FID, '\\textcolor{COWI}{\\tiny{PROJECT NO.}} & \\textcolor{COWI}{\\tiny{DOCUMENT NO.}} & \\multicolumn{2}{l}{\\textcolor{COWI}{\\tiny{%s DOCUMENT NO.}}} \\\\ \n',Settings.Appendix.ProjectEmployer);
%         fprintf(FID, '&  &  &  \\\\ \n');
%         fprintf(FID, '%s & %s & \\multicolumn{4}{l}{%s} \\\\ \n',Settings.Appendix.ProjectNumber,Settings.Appendix.DocumentNoCOWI,Settings.Appendix.DocumentNoEmployer);
%     else % If only COWI number is specified
        fprintf(FID, '\\textcolor{COWI}{\\tiny{PROJECT NO.}} & \\textcolor{COWI}{\\tiny{DOCUMENT NO.}} & \\multicolumn{2}{l}{\\textcolor{COWI}{ }} \\\\ \n');
        fprintf(FID, '&  &  &  \\\\ \n');
        fprintf(FID, '%s & %s & & \\\\ \n',Settings.Appendix.ProjectNumber,Settings.Appendix.DocumentNoCOWI);
%     end
    fprintf(FID, '&  &  & \\\\ \n');
    fprintf(FID, '&  &  & \\\\ \n');
    fprintf(FID, '\\end{tabular} \n \n');
    
    
    if any(cellfun('length',Settings.Appendix.RevisionTable(:,3))>82)
        fprintf(FID, '\\begin{tabular}{ll p{8cm} lll} \n');
    else
        fprintf(FID, '\\begin{tabular}{llllll} \n');
    end
    fprintf(FID, '\\textcolor{COWI}{\\tiny{VERSION}} & \\textcolor{COWI}{\\tiny{DATE OF ISSUE}} & \\textcolor{COWI}{\\tiny{DESCRIPTION}} & \\textcolor{COWI}{\\tiny{PREPARED}} & \\textcolor{COWI}{\\tiny{CHECKED}} & \\textcolor{COWI}{\\tiny{APPROVED}}\\\\ \n');
    fprintf(FID, '&  &  &  &  & \\\\ \n');
    fprintf(FID, '%s & %s & %s & %s & %s & %s \\\\ \n',...
        Settings.Appendix.RevisionTable{1,1}, Settings.Appendix.RevisionTable{1,2},...
        Settings.Appendix.RevisionTable{1,3}, Settings.Appendix.RevisionTable{1,4},...
        Settings.Appendix.RevisionTable{1,5}, Settings.Appendix.RevisionTable{1,6});
    fprintf(FID, ' & & & & & \\\\ \n');
    fprintf(FID, '\\end{tabular} \n');
    fprintf(FID, '\\end{minipage}} \n');
    fclose(FID);
    
    
    FID = fopen(char({'AppendixGenerationFiles\TemplateFiles\Main_Body.tex'}), 'w+');
    
    if Settings.Appendix.GeneralInfo ==1;
        fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/General_info} \n');
    end 
    
    if Settings.Appendix.Soil_input ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Soil} \n');
    end    
        
    if Settings.Appendix.Pile_Geometry ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Pile_geometry_Chapter} \n');
    end    
        
    if Settings.Appendix.Lateral_displacement ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Lat_disp_Section} \n');
    end    
        
    if Settings.Appendix.Crit_pile_length ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Crit_pile_length_Section} \n');
    end    
             
    if Settings.Appendix.Natural_frequency ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Natural_frequency_Section} \n');
    end    
    
    if Settings.Appendix.Perm_rotations==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Perm_rotations_Section} \n');
    end     
    
    if Settings.Appendix.Lat_pile_capacity ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Lat_pile_capacity_Section} \n');
    end 
        
    if Settings.Appendix.Axial_pile_capacity ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/Axial_pile_capacity_Section} \n');
    end  
    
    if Settings.Appendix.SSI_springs ==1;
            fprintf(FID, '\\include{AppendixGenerationFiles/TemplateFiles/SSI_Springs_Section} \n');
    end  
    
        fclose(FID);
    
    
    
    
    
    
    
    