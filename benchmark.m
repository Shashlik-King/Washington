%% BENCHMARKING 
% Script comparing the benchmark model of the original version to the checked out version from the repo
clc; close all; clear all;
cprintf('black','Benchmarking initiated \n');
cprintf('black','---------------------------------------------------------------------- \n');

%% Run Inverse model of new version
% system('taskkill /F /IM EXCEL.EXE');
copyfile('input/Initialise.m', 'temp/Initialise_temp.m','f') % move COSPIN input to a temporary folder
% Unit test 1 - reaction curve creator
copyfile('benchmark\Initialise_unit_test_1.m', 'input\Initialise.m','f') % replace normal input by benchmarking input
run_COSPIN % run Inverse model
% Unit test 2 - pile response output
copyfile('benchmark\Initialise_unit_test_2.m', 'input\Initialise.m','f') % replace normal input by benchmarking input
run_COSPIN % run Inverse model
copyfile('temp\Initialise_temp.m', 'input\Initialise.m','f') % copy original input file back to the main reading folder

%% Settings
tolerance = 0.000000001;

%% Benchmarking
path_benchmark_ref  = 'benchmark\';                         % path for the reference datasets
path_benchmark_new  = 'output\';             % path for the newly ran datasets
files = {'deflection','disp','load'}; % list fo files used for comparison
for i = 1:size(files,2)
    %% Read benchmark model and new model
    file_name_ref   = [path_benchmark_ref,files{i},'.txt']; % total path for reference dataset
    file_name_new   = [path_benchmark_new,files{i},'.txt']; % total path for new dataset
    reference       = importdata(file_name_ref);            % import of reference dataset
    new             = importdata(file_name_new);            % import of new dataset
    
    %% Calculate error
    error.response{i,1}       = (new - reference) ./ reference; % calculate error for each discretised element/value
    % tolerance implementation
    for j = 1:size(error.response{i,1},1)
        for k = 1:size(error.response{i,1},2)
            if abs(error.response{i,1}(j,k)) <= tolerance
                error.response{i,1}(j,k) = 0;
            end
        end
    end
    error.response{i,1}(isnan(error.response{i,1}))=0;
    total_error.response{i,1} = sum(sum(error.response{i,1}));  % calculate total error for file analysed
    cprintf('black',['Total error results for ',files{i}, ' are:']);
    if total_error.response{i,1} == 0
        cprintf('green', 'Acceptable. \n');
    else
        cprintf('red', 'Not acceptable. \n');
    end
end

global_error.response = sum([total_error.response{1:end,1}]);
cprintf('black','Benchmarking test results for pile response are:');
if global_error.response == 0
    cprintf('green', 'Acceptable.  \n');
else
    cprintf('red', 'Not acceptable.  \n');
end
cprintf('black','Benchmarking for pile response is complete \n');
cprintf('black','---------------------------------------------------------------------- \n');


writematrix([total_error.response{1:end,1}],'total_error.txt','Delimiter','tab')