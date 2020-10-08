
%% This script runs the experiment of performing GC estimation compared to Brainstorm
% The EEG data is generated under the setting: 50% deep source, 20% active source, 61 No. of electrodes
% Input: EEG data
% Output: Estimated state-space parameters
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

inpath = './input_data/source_recon_brainstorm/eeg_simulated_data_deep50percent_active20percent';
outpath = './saved_experiment_results/simeeg_sourcerecon/';


nbatch = 1;
n_model = 100;
n_deep = 2; % run only the case that there are 2 deep source ROIs
n_active_source = 10; % run the case that there are 10 active sources (from 50 sources)

%% Divide works

for mm = 1:n_model
    postfix = ['act_',int2str(n_active_source),'_deep_',int2str(n_deep),'_',int2str(mm)];
    disp(postfix)
    load([inpath,'eegdata_',postfix])
    load([inpath,'model_',postfix])
    data = eegdata.EEG_data;
    M = reduce_L(eegdata.L);
    L = M.L10;
    ind_chan = M.L10_ind;
    mtilde = size(L,2);
    Timepoints = 15000;
    tic;
    for tt=1:nbatch % trial iteration
        data_in = data(ind_chan,end-Timepoints+1:end,tt);
        [sys_est(tt).sys,sys_est(tt).C_out] = ...
            subid_eeg_Lpq(data_in ...
            ,L,10,-1);
    end
    toc;
    save([outpath,'simeeg_sourcerecon_sysest_',postfix],'sys_est')
    clear sys_est
end

