
%% This script generates simulated EEG time series data used for a comparion 
% with source recontruction algorithms in Brainstorm under the setting:
% 50% deep sources and 20% active sources
% Total 100 models are generated (with underlying GC structures) and corresponding time series
% The active sources are randomly chosen to lie in 4 ROIs (out of 8 ROIs).
% Written by Parinthorn Manomaisaowapak
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
%% data generation for brainstorm & our method
clc
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('./input_data/nyhead_model')
end
foldername = 'source_recon_brainstorm';

%% User Input
nbatch = 1; % number of realizations of EEG time series
n_model = 100; % number of generated models

% n_source (number of all sources) : now feasible up to hundred
% n_source_active (number of active source): must be less than n_source,
% now feasible up to 100 sources with 5-10% sparsity of underlying GC
% n_deep_roi (number of deep source ROI): from 0 to 4 ROIs

SOURCE_PARAMETER = [50 10 2]; % [n_source n_source_active n_deep_roi]

%% Generate models
for mm=1:n_model
            postfix = ['_act_',int2str(SOURCE_PARAMETER(2)),'_deep_',int2str(SOURCE_PARAMETER(3)),'_',int2str(mm)];
            disp(postfix)
            % no plot, pink-noise
            generate_eegdata(nbatch,foldername,0,postfix,sa,1,SOURCE_PARAMETER);
end

