function [element,soil] = interimRotationalMultipliers (element,soil)

soil.degradation.value_mt_m = ones(size(soil.degradation.value_py_p)); 
soil.degradation.value_mt_t = ones(size(soil.degradation.value_py_p));

%degrdation for Till
soil.degradation.value_mt_m(or(strcmp(soil.description,'4a')', ...
    strcmp(soil.description,'4b')')) = 0.47;
%degradation for Chalk
soil.degradation.value_mt_m(or(strcmp(soil.description,'5a')', ...
    strcmp(soil.description,'5b')')) = 0.21;

element.degradation_mt_m = ones(size(element.degradation_py_p));
element.degradation_mt_m(strcmp(element.model_py,'PISA Till')) = 0.47;
element.degradation_mt_m(strcmp(element.model_py,'PISA Chalk')) = 0.21;

end