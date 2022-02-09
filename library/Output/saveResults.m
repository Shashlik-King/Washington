function saveResults(data,element,loads,output,pile,settings,SF)

filename =[data.save_path,'\results.xlsx'];

if settings.lateral_loading
    xlswrite(filename, M, sheet, 'range')
end
end