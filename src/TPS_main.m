function [] = TPS_main(input_test_file, varargin)

%This function dates ancient DNA by using the tool Temporal Population Structure (TPS) (see Esposito et al. "A genomic dating tool for ancient genomes resolves the origins of hundreds of Eurasian genomes")

% Mandatory input:
%    input_test_file -  xlsx file with the aDNA dataset to date with TPS (it needs to follow a specific format)
%               A test file is provided which can be used to run TPS: "../data/TPS_test_ids.xlsx"
% Optional inputs:
%    input_ref_file - xlsx file with the aDNA dataset (it needs to follow a specific format)
%               Default value: "../data/TPS_reference_ids.xlsx"
%    n_rand_repetitions - Number of random training/testing repetitions
%               Default value: 500
%    training_size - Number of individuals used to build the reference panel for each random iteration
%               Default value: 400
%
% Outputs (stored in the folder ../results/):
%    Output1 - xlsx file with TPS results
%
% How to call this script:
%   - TPS_main('../data/TPS_test_ids.xlsx') to test the tool with the given file
%   - TPS_main('Your_filename.xlsx') - Check the given file TPS_test_ids.xlsx in ../data/ to correctly format your file
%   - NOTE: It is not recommended to apply TPS to DNA older than 14,000 BP
%
% Subfunctions called: 
%   Select_reference_ids.m
%   Build_reference_panel.m (which calls mpFtest.m)
%   Run_TPS.m
%
% Author: Umberto Esposito
% email: umberto.espos@gmail.com
% Created: October 2019
% Last edited: 30 October 2019


%% Assign input arguments
%Default values
inputs_defaults = {'../data/TPS_reference_ids.xlsx', 500, 400};
%'../data/TPS_test_ids.xlsx'

%Assign values to use
inputs = inputs_defaults;
idx = ~cellfun('isempty', varargin);
inputs(idx) = varargin(idx);

%Define input variables
input_ref_file = inputs{1};
n_rand_repetitions = inputs{2};
training_size = inputs{3};


%% Model parameters - Reference panel
%For the reference panel consider only individuals younger than 14,000 BP
ref_thresh = 14000;


%% Read inputs
%Read reference data
opts_ref = detectImportOptions(input_ref_file);
Dataset_ref = readtable(input_ref_file,opts_ref);

%Read test data
opts_test = detectImportOptions(input_test_file);
Dataset_test.data = readtable(input_test_file,opts_test);


%% Variables and data
%Select data to use according to the parameters'value chosen above
IDs_toDiscard_logicalInd = Dataset_ref.DateBP > ref_thresh;
Dataset_ref(IDs_toDiscard_logicalInd,:) = [];

%Admixture components in the input spreadsheet (Columns F-M)
adComponents_names = {'AncientComponent1','AncientComponent2','AncientComponent3','AncientComponent4','AncientComponent5','ModernComponent1','ModernComponent2','ModernComponent3'};

%Group individuals into populations by rounding to the nearest 500
Population = strcat('Date_',string(round(Dataset_ref.DateBP*2,-3)/2),'BP');
Population(ismissing(Population)) = 'No_Date';
Dataset_ref = addvars(Dataset_ref,Population,'After','ID');

%Output folder
output_folder = '../results/';

%Admixtute coefficients of the test individuals
adCoef_test = table2array(Dataset_test.data(:,adComponents_names));

%% TPS bootstrapping routine
fprintf('\nRoutine started:\n\n');

for i_rep = 1:n_rand_repetitions

    %Print progress on screen
    t_iter = tic;
    fprintf('\tIteration #%d\n',i_rep);

    %Randomly choose training_size individuals among the reference ones
    [Dataset_training, ~] = Select_reference_ids(Dataset_ref, training_size);

    %Create the TPS reference panel
    adCoef_training = table2array(Dataset_training(:,adComponents_names));
    [~,exp_txt_GEN,exp_txt_TEM,exp_GEN,exp_TEM,exp_TEM_STD] = Build_reference_panel(Dataset_training, adCoef_training, 'Off');

    %Control statement
    if any(strcmp(exp_txt_GEN,exp_txt_TEM) == 0)
        error('Reference panel populations names are not the same');
    end

    %Run TPS predictive routine
    [predicted_time,predicted_time_STD,min_gen_dist] = Run_TPS(exp_txt_GEN, exp_GEN, exp_TEM, exp_TEM_STD, adCoef_test, 'Off');
    Dataset_test.TPSrawResults.SinglePredictions(:,i_rep) = predicted_time;
    Dataset_test.TPSrawResults.StdSinglePredictions(:,i_rep) = predicted_time_STD;
    Dataset_test.TPSrawResults.MinGenDistSinglePredictions(:,i_rep) = min_gen_dist;

    toc(t_iter)
end

%% TPS individuals performance (averages across single runs)
%Assess TPS performance on the test datasets' individuals
%---TPS average quantities
%Average prediction
Dataset_test.TPSresults.TPSpredictions = nansum(Dataset_test.TPSrawResults.SinglePredictions ./ (Dataset_test.TPSrawResults.StdSinglePredictions.^2),2) ./ nansum(1 ./ (Dataset_test.TPSrawResults.StdSinglePredictions.^2),2);
%Standard deviation
Dataset_test.TPSresults.TPSstdPredictions = sqrt(nansum( (1 ./ (Dataset_test.TPSrawResults.StdSinglePredictions.^2)) .* ((Dataset_test.TPSrawResults.SinglePredictions - Dataset_test.TPSresults.TPSpredictions).^2), 2) ./ (nansum(1 ./ (Dataset_test.TPSrawResults.StdSinglePredictions.^2),2)));
%Confidence intervals
for i=1:size(Dataset_test.TPSrawResults.SinglePredictions,1)
    SinglePredictions = sort(Dataset_test.TPSrawResults.SinglePredictions(i,~isnan(Dataset_test.TPSrawResults.SinglePredictions(i,:))));
    if isempty(SinglePredictions)
        Dataset_test.TPSresults.TPSCIupper(i,:) = nan;
        Dataset_test.TPSresults.TPSCIlower(i,:) = nan;
    else
        margin = (length(SinglePredictions) * 5 / 100);
        upper_element = length(SinglePredictions)-ceil(margin)+1;
        upper_element = (upper_element>length(SinglePredictions)) .* length(SinglePredictions) + (upper_element<=length(SinglePredictions)) .* upper_element;
        lower_element = floor(margin);
        lower_element = (lower_element<1) .* 1 + (lower_element>=1) .* lower_element;
        Dataset_test.TPSresults.TPSCIupper(i,:) = SinglePredictions(upper_element);
        Dataset_test.TPSresults.TPSCIlower(i,:) = SinglePredictions(lower_element);
    end
end
%Minimum genetic distance with the reference panel
Dataset_test.TPSresults.TPSminGenDist = nanmean(Dataset_test.TPSrawResults.MinGenDistSinglePredictions,2);


%% Save & print on file
%Create folder if needed
if exist(output_folder, 'dir') ~= 7
    mkdir(output_folder);
end

%Save mat file
save([output_folder,'TPS_results.mat']);

%Delete existing file
outputFilename = ['TPS_results_',num2str(n_rand_repetitions),'repetitions.xlsx'];
if exist([output_folder,outputFilename], 'file') == 2
    delete([output_folder,outputFilename]);
end

%Write xlsx results file
writetable([Dataset_test.data,struct2table(Dataset_test.TPSresults)],[output_folder,outputFilename],'Sheet','TPSresults')%,'Range','A1:N1000');
