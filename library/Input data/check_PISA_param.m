function [element] = check_PISA_param(element)

for i = 1:element.nelem
    if strcmp(element.model_py(i),'Zero soil') == 0
%% Negative values check
    for j = 1:7
        if element.PISA_param(i,j) < 0 % check fro negative values
%             error([' p-y parameter in node ' num2str(i) ' is negative - see documentation'])
            element.PISA_param(i,j) = 0.001;
        end
    end

    for j = 8:14
        if element.PISA_param(i,j) < 0 % check fro negative values
%             error([' m-t parameter in node ' num2str(i) ' is negative - see documentation'])
            element.PISA_param(i,j) = 0.001;
        end
        
    end
    
    for j = 15:18
        if element.PISA_param(i,j) < 0 % check fro negative values
%             error([' Hb parameter in node ' num2str(i) ' is negative - see documentation'])
            element.PISA_param(i,j) = 0.001;
        end
        
    end
    
    for j = 19:22
        if element.PISA_param(i,j) < 0 % check fro negative values
%             error([' Mb parameter in node ' num2str(i) ' is negative - see documentation'])
            element.PISA_param(i,j) = 0.001;
        end
        
    end
%% Curvature value check    
    if element.PISA_param(i,6) == 1 || element.PISA_param(i,7) == 1 || element.PISA_param(i,13) == 1 || element.PISA_param(i,14) == 1 || element.PISA_param(i,18) == 1 || element.PISA_param(i,22) == 1
        warning([' Curvature parameter in node ' num2str(i) ' is = 1 causing convergence errors - see documentation'])
    end
%% Ultimate lateral displacement value check      
%     if element.PISA_param(i,1) == 0 || element.PISA_param(i,8) == 0 || element.PISA_param(i,15) == 0 || element.PISA_param(i,19) == 0
%         error([' Ultimate lat. disp. parameter in node ' num2str(i) ' is = 0 causing complex results - see documentation'])
%     end
    end
end
