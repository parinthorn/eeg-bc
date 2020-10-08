% This script demonstrates how to repeat the experiment
% Written by Jitkomut Songsiri, Parinthorn Manomaisaowapak, Anawat Nartkulpat
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
% Assume that we run these codes on the parent folder of eeggc-master

%% Experiment: SIMEEG_VARYFACTOR
% Examine the performance of our method

simeeg_varyfactor_subid

simeeg_varyfactor_gc_processing

simeeg_varyfactor_gcest_calperformance


%% Experiment: SIMEEG_SOURCERECON
% compare the source selection and GC estimation between our method and
% brainstorm

addpath('experiment/');

% Generate source amplitude from head model
% input = head model, ground-truth models
% output = source amplitude of all sources in consideration

simeeg_sourcerecon_source_amplitude

% our method

% prepare ground-truth/estimate source cluster index refering to 4 ROIs 
% and calculate permutation indices to map the index to the 8 ROis with the same order as Hauf's 
% 
% input = ground-truth model
% output = ground-truth/estimate source cluster index and permutation index (both node-wise and ROI-wise)

simeeg_sourcerecon_adjust_roi

% input = eeg data + ground-truth model
% output = estimated state-space parameter
simeeg_sourcerecon_ss_method 

% input = estimated state-space parameter, ground-truth models
% output = estimated C of 100 models with BIC indices

simeeg_sourcerecon_C_processing

% input = estimated state-space paramter
% output = estimated GC (ordered in Haufe ROI format)
simeeg_sourcerecon_ss_method_permuteGC 

% Brainstorm

% input = eeg data
% output = kernel matrices of WMNE, sLORETA, LCMV (to recover reconstructed source)
simeeg_sourcerecon_brainstorm_main 

% input = head model vertices
% output = matrix transform from 75K to 8 ROIs

simeeg_sourcerecon_brainstorm_75Ksource_to_8roi  

% input = kernel matrices of WMNE, sLORETA, LCMV,  AND matrix transform from 75K to 8 ROIs
% output = averaged reconstructed source in ROI
simeeg_sourcerecon_brainstorm_roi_timeseries_extraction

% input = averaged reconstructed source in ROI
% output = Granger causality matrix of reconstucted source data from Brainstorm
gcest_sourcerecon_mvgc_gc

% input = reconstruct source, kernel matrix
% output = heat map of averaged absolute of recontructed source (over time)
% in 75K resolution
simeeg_sourcerecon_brainstorm_sourcesignal_reconstruction 

%% Experiment:  SSVEP EEG (real data set)

% We choose to run either L21 or Lpq formulation

% convex penalty in source selection process
ssvep_data_main_L21;

% non-convex penaly in source selection process
ssvep_data_main_Lpq;

