% creation of Excel file with input to soil damping
function soil_damping_input_write(data,Coord,output,loads);
%--------------------------------------------------------------------------
% PURPOSE
% To automate the procedure of writing soil damping input to excel 

% Created by      : FKMV 
% Created         : 03/02/2020
% LAST MODIFIED   : 
%--------------------------------------------------------------------------
sheet_name = [loads.type,' ',data.location];
excelname = '/damping/Soil_damping/2._Input/Soil_Damping_Input.xlsx';
excelname = 'LOC30_FLS_Soil_damping_calculation_mode_1_3.xlsx'

xlswrite(excelname ,Coord(:,2), sheet_name, 'A3')
xlswrite(excelname ,output.hor_defl, sheet_name, 'B3')
xlswrite(excelname ,output.ks(:,1), sheet_name, 'C3')
xlswrite(excelname ,output.ks(:,2), sheet_name, 'D3')
xlswrite(excelname ,cellstr('Coordinate of the node'), sheet_name, 'A1')
xlswrite(excelname ,cellstr('Horizontal deflection'), sheet_name, 'B1')
xlswrite(excelname ,cellstr('p-y secant stiffness'), sheet_name, 'C1')
xlswrite(excelname ,cellstr('p-y secant stiffness'), sheet_name, 'D1')
xlswrite(excelname ,cellstr('[m]'), sheet_name, 'A2')
xlswrite(excelname ,cellstr('[m]'), sheet_name, 'B2')
xlswrite(excelname ,cellstr('[kN/m/m]'), sheet_name, 'C2')
xlswrite(excelname ,cellstr('[kN/m/m]'), sheet_name, 'D2')

xlswrite(excelname ,cellstr('Horizontal Force applied'), sheet_name, 'F1')
xlswrite(excelname ,cellstr('Moment applied'), sheet_name, 'G1')
xlswrite(excelname ,cellstr('kN'), sheet_name, 'F2')
xlswrite(excelname ,cellstr('kNm'), sheet_name, 'G2')
xlswrite(excelname ,loads.H, sheet_name, 'F3') 
xlswrite(excelname ,loads.M, sheet_name, 'G3') 
 
disp(['Damping for location ',data.location,' finished'])
end