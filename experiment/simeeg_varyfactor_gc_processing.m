%% Simulation: Varying three factors
% 
% This script runs an experiment of varying three factors: 
% 1) percentage of deep sources: 0%, 50%, 75%
% 2) active source density: 20%, 40%
% 3) number of electrodes: 108, 61, 32, 19
% 
% The number of total combinatorial cases = 3 x 2 x 4 = 24 cases
% In each case, we test on 100-run simulated EEG time series with  underlying known
% ground-truth GC
% 
% Input: SYSEST contains estimated state-space parameters in each of 24 cases
% Output: Estimated GC matrices: F in structure variable M
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

addpath('./input_data/eeg_simulated_data/'); % to read lead-field matrix from ground-truth  model

inpath = './saved_experiment_results/';
outpath = './saved_experiment_results/';


%%
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [0 2 3];
density = [10 20];
density_it = 0;
for n_act = density
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for kk=length(Lname)
            n_elec = n_elec+1;

            for mm = 1:n_model
                postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
                load(['model_',postfix])
                L = model.(Lname{kk});
                disp(postfix)
                load([inpath,'simeeg_varyfactor_sysest_',postfix,'_',Lname{kk}],'sys_est')
                r = size(L,1);
                for tt=1:nbatch                    
                    C_ind = sys_est(tt).C_out.ind_chosen_Lpq.bic;% use bic score
                    A = sys_est(tt).sys.A;
                    C = sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind);
                    E = sys_est(tt).C_out.E_Lpq_CLS(:,:,C_ind);
                    W = sys_est(tt).C_out.W_Lpq_CLS(:,:,C_ind); 
                    [N,V,obj] = noisecovest(E,L,'homo');

                    [Ftmp,~] = calgcss(A,C,W,N,[]);
                    F{mm,tt} = sparse(Ftmp);
                    F_nz_ind{mm,tt} = find(Ftmp);
                end
                F_true{mm,1} = model.F0;
                F_true_ind{mm,1} = find(model.F0);
            end
            M.F_true{deep_it,density_it} = F_true;
            M.F_true_ind{deep_it,density_it} = F_true_ind;
            M.F_nz_ind{deep_it,density_it,n_elec} = F_nz_ind;
            M.F{deep_it,density_it,n_elec} = F;
        end
    end
end

%%
save([outpath,'simeeg_varyfactor_Fresult','M']);

