% Simulation data experiment: 50% deep source
%                             20% active source
%                             61 No. of electrodes
% This script generate a mapping from 75K sources to 8 rois sorting in Haufe's ordering
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

inpath = './input_data/nyhead_model/tess_cortex_01.mat';
outpath = './saved_experiment_results/simeeg_sourcerecon/brainstorm_directory/';
load(inpath,'Vertices')
%%
row_ind = 1:1:74661;
x0 = 0.0095;
y0 = 0.0001;
z0 = 0.060;
roi_bin = double([Vertices(row_ind,1)>x0 Vertices(row_ind,2)>y0 Vertices(row_ind,3)>z0]);

% 000 -> RPI 1
% 001 -> RPS 2
% 010 -> LPI 3
% 011 -> LPS 4
% 100 -> RAI 5
% 101 -> RAS 6
% 110 -> LAI 7
% 111 -> LAS 8

roi_ind = bin2dec(num2str(roi_bin))+1;
roi_bin(sum(roi_bin,2)==3,:)= (roi_bin(sum(roi_bin,2)==3,:))-[0.5 0.5 0.5];
figure(1)
scatter3(Vertices(row_ind,1),Vertices(row_ind,2),Vertices(row_ind,3),[],roi_bin)
figure(2)
histogram(roi_ind)

%% projection

LPS = double(roi_ind==4); % haufe ordering
LPI = double(roi_ind==3);
LAS = double(roi_ind==8);
LAI = double(roi_ind==7);
RPS = double(roi_ind==2);
RPI = double(roi_ind==1);
RAS = double(roi_ind==6);
RAI = double(roi_ind==5);
Kernel_roi = [LPS LPI LAS LAI RPS RPI RAS RAI]';
save([outpath,'Kernel_roi'],'Kernel_roi')
