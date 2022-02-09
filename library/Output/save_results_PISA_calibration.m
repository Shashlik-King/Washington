function save_results_PISA_calibration(data,element,output,plots,Es)
excelname = 'COSPIN_results.xlsx';

if plots.moment_distribution == 1
    xlswrite(excelname ,Es{1,1}(:,1,3), data.location, 'A2')
    xlswrite(excelname ,Es{1,1}(:,1,2), data.location, 'B2')
    xlswrite(excelname ,element.level(:,1), data.location, 'D2')
    xlswrite(excelname ,cellstr('Moment'), data.location, 'A1')
    xlswrite(excelname ,cellstr('Shear'), data.location, 'B1')
    xlswrite(excelname ,cellstr('Depth'), data.location, 'D1')
    xlswrite(excelname ,output.hor_defl(1:end-1,end), data.location, 'C2')
    xlswrite(excelname ,cellstr('Displacement'), data.location, 'C1')
end

if plots.load_deflection == 1
    xlswrite(excelname ,output.def_calibration', data.location, 'F2')
    xlswrite(excelname ,output.force_calibration', data.location, 'G2')
    xlswrite(excelname ,cellstr('Displacement'), data.location, 'F1')
    xlswrite(excelname ,cellstr('Force'), data.location, 'G1')
end