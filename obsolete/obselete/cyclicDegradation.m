function [element] = cyclicDegradation (element)
for i = 1:length(element.cu)
    element.cu(i,:) = element.cu(i)*0.72;
end
end