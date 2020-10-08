%% Simulation data experiment: 50% deep source
%                              20% active source
%                              61 No. of electrodes
% This script generate the permutation index to guarantee Haufe's ordering in the ss method
% Written by: PARINTHORN MANOMAISAOWAPAK
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
inpath = './input_data/source_recon_brainstorm/eeg_simulated_data_deep50percent_active20percent/';
outpath = './saved_experiment_results/simeeg_sourcerecon/';
%%
Lname = {'L10'};
n_model = 100;
n_deep = 2;
n_active_source = 10;

for mm = 1:n_model
    postfix = ['_act_',int2str(n_active_source),'_deep_',int2str(n_deep),'_',int2str(mm)];
    load([inpath,'model_',postfix])
    [~,permind{mm,1}]=sort([unique(model.truth.in_roi),unique(model.truth.augment_data.in_roi)]);
end
save([outpath,'permind'],'permind')
