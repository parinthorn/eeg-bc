%% Simulation: Varying three factors
% 
% This script runs an experiment of varying three factors: 
% 1) percentage of deep sources: 0%, 50%, 75%
% 2) active source density: 20%, 40%
% 3) number of electrodes: 108, 61, 32, 19
% 
% The number of total combinatorial cases = 3 x 2 x 4 = 24 cases
% In each case, we test with 100-run simulated data with underlying known
% ground-truth GC
% 
% Input: EEGDATA in each of 24 cases saved in 'inpath'
% Output: Estimated State-space parameters: SYSEST; will be saved in
% 'outpath'
% 
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

clear
clc

inpath = './input_data/eeg_simulated_data/';
outpath = './saved_experiment_results/';

%% User Input

nbatch = 1; % number of time series trial per one ground-truth model
n_model = 100; % number of ground-truth model in the generated folder
n_deep = [0 2 4]; % number of deep source ROI (for reading filenames)
n_active_source = [10 20]; % active source density (for reading filenames)

%% Divide works

% L0: original L in 10-5 system (no compensation of signal scaling, transformation, etc.)
% L10: 10-10 system 
% L20_32: 10-20 system with 32 channels 
% L20: 10-20 system with 19 channels

Lname = {'L0','L10','L20_32','L20'}; 
for n_act = n_active_source
    for nd = n_deep
        for mm = 1:n_model
            postfix = ['act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
            disp(postfix)
            load([inpath,'eegdata_',postfix])
            load([inpath,'model_',postfix])
            data = eegdata.EEG_data;
            for LL=Lname
                L = model.(LL{1});
                r = size(L,1);
                if strcmp(LL{1},'L0')
                    ind_chan = 1:1:r;
                elseif strcmp(LL{1},'L20_32')
                    ind_chan = model.L20_ind_32;
                else
                    ind_chan = model.([LL{1},'_ind']);
                end
                Timepoints = 15000; % our generated EEG has 18,000 time points but we trimed the initial segment
                tic;
                for tt=1:nbatch % trial iteration
                    data_in = data(ind_chan,end-Timepoints+1:end,tt);
                    [sys_est(tt).sys,sys_est(tt).C_out] = ...
                        subid_eeg_Lpq(data_in,L,10,-1);  
                    disp([outpath,'result_',postfix,'_',LL{1},int2str(tt)])
                end
                toc;
                save([outpath,'simeeg_varyfactor_sysest_',postfix,'_',LL{1}],'sys_est');
                clear sys_est
            end
        end
    end
end
