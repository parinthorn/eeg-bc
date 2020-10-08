%%   This script demonstrates generating ground truth models and their
%   corresponding time series data.
%
%   Written by ANAWAT NARTKULPAT, PARINTHORN MANOMAISAOWAPAK, 2020
%
% Licenses
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%

clc;

addpath('data_generation');
addpath('data_generation/haufe_tools');
addpath('data_generation/haufe_tools/matlab_bgl');
addpath('data_generation/haufe_tools/export_fig');

% load head model if it is not already loaded. The head model is then
% contained in varable 'sa'
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('./input_data/nyhead_model/sa')
    load('./input_data/nyhead_model/miscdata')
end

%% User input
n_model = 3;    % number of generated model
nbatch = 2;	% number of time series per one ground truth model
foldername = 'eeg_simulated_data'; % name of folder containing model and time series to be saved
plotting = 0;	% 1 for plot and save figures originally by Haufe, 0 otherwise

is_pink = 1;    % 1 for perturb system with pinknoise, 0 for using white noise instead
% there are 8 ROIs, 4 ROIs are deep areas. We generate active sources from
% 4 ROIs 

% n_source (number of all sources) : now feasible up to hundred
% n_source_active (number of active source): must be less than n_source,
% now feasible up to 100 sources with 5-10% sparsity of underlying GC
% n_deep_roi (number of deep source ROI): from 0 to 4 ROIs

SOURCE_PARAMETER = [50 10 2]; % [n_source n_source_active n_deep_roi]

%% Generate model and time series
for ii=1:n_model
    
    postfix = num2str(ii);	% in 'foldername', model and time series data are saved 
             	% and named model_postfix.mat and eegdata_postfix.mat
                
    fprintf('Generating Model:  %d \n',ii)
    % generated and saved model and saved time series data
    generate_eegdata(nbatch,foldername,plotting,postfix,sa,is_pink,SOURCE_PARAMETER); 
end

