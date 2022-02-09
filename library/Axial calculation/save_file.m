function [output] = save_file(settings, element, output, p, y, WTG) 
%% CALCULATION OF STATIC AXIAL CAPACITY

if settings.save_file == 1
    for i=1:length(element.level)-1
        output.excel(i*2-1+2,1)=element.level(i,1);
        output.excel(i*2-0+2,1)=element.level(i,1);
        output.excel(i*2-1+2,3:19)=p.top(i,:);
        output.excel(i*2-0+2,3:19)=y.top(i,:);
    end
    output.excel(i*2+1+2,1)=element.level(i+1,1);
    output.excel(i*2+2+2,1)=element.level(i+1,1);
    output.excel(i*2+1+2,3:19)=p.bottom(i,:);
    output.excel(i*2+2+2,3:19)=y.bottom(i,:);
    fname = '\\cowi.net\projects\A125000\A127858\GEO\Calculations\COSPIN_calcs\8. 21oct2019 prelim soil springs\results\';
    filename = '\\cowi.net\projects\A125000\A127858\GEO\Calculations\COSPIN_calcs\8. 21oct2019 prelim soil springs\results\P-y.xlsx';
    xlswrite(filename ,output.excel, WTG{loc});
    saveas(gca, fullfile(fname, WTG{loc}), 'jpeg');
    saveas(gca, fullfile(fname, WTG{loc}), 'fig');
else
end
end