function Excel_results_write(Settings,i,AppendixData)

%% Lateral table write
if Settings.Appendix.Perm_rotations == 1 || Settings.Appendix.Lat_pile_capacity == 1
   xlswrite('Results_summary_COSPIN.xlsx', {AppendixData.locationID} ,'Lateral', strcat('A',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.pileLength ,'Lateral', strcat('D',num2str(i+2)))
end

if Settings.Appendix.Perm_rotations == 1
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.perm_rot*180/pi,'Lateral', strcat('E',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', 0.300 ,'Lateral', strcat('F',num2str(i+2)))
end

if Settings.Appendix.Lat_pile_capacity == 1
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.BSH ,'Lateral', strcat('G',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.DNV_GL ,'Lateral', strcat('H',num2str(i+2)))
end



%% Axial capacity table write
if Settings.Appendix.Axial_pile_capacity == 1
   xlswrite('Results_summary_COSPIN.xlsx',  {AppendixData.locationID} ,'Axial Capacity', strcat('A',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,1) ,'Axial Capacity', strcat('D',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,2) ,'Axial Capacity', strcat('E',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,3) ,'Axial Capacity', strcat('F',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,4) ,'Axial Capacity', strcat('G',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,5) ,'Axial Capacity', strcat('H',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_cap(1,6) ,'Axial Capacity', strcat('I',num2str(i+2)))
end

%% Axial settlement table write
if Settings.Appendix.Axial_pile_settlement == 1
   xlswrite('Results_summary_COSPIN.xlsx', {AppendixData.locationID} ,'Axial Settlements', strcat('A',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_set(1,1) ,'Axial Settlements', strcat('D',num2str(i+2)))
   xlswrite('Results_summary_COSPIN.xlsx', AppendixData.results.axial_set(1,2) ,'Axial Settlements', strcat('E',num2str(i+2)))
end