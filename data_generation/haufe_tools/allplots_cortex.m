function allplots_cortex(sa, data, colorlimits, cm, unit, smooth, printfolder, varargin)
% Stefan Haufe, 2014, 2015
% stefan.haufe@tu-berlin.de

% If you use this code for a publication, please ask Stefan Haufe for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

set(0,'DefaultFigureColor',[1 1 1])
printfolder = [printfolder '/'];
if ~exist(printfolder, 'dir')
    mkdir(printfolder)
end

res = '150';

if smooth 
    vc = sa.cortex75K.vc_smooth;
    sm = '_smooth';
else
    vc = sa.cortex75K.vc;
    sm = '';
end

surface_pars = struct('alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []);

if length(varargin) > 0
   varargin1 = varargin{1};
else
    varargin1 = {};
end

if length(varargin) > 1
    input_pars = varargin{2};
    finames = fieldnames(input_pars);
    for ifi = 1:length(finames)
      surface_pars = setfield(surface_pars, finames{ifi}, getfield(input_pars, finames{ifi}));
    end
end
    
surface_pars.myviewdir = [-1 0 0];
figure; showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left'], ['-r' num2str(res)], '-a2'); 


figure; showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right_inner'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [1 0 0];

figure; showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left_inner'], ['-r' num2str(res)], '-a2'); 

figure; showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [-1e-10 0 1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [0 0 1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top_upright'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [-1e-10 0 -1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [0 1e-10 -1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom_upright'], ['-r' num2str(res)], '-a2'); 


figure; 
hf = imagesc(randn(5)); colormap(cm)
set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
set(hf, 'visible', 'off')
cb = colorbar; 
set(cb, 'fontsize', 30)
ylabel(cb, unit)
export_fig([printfolder 'cortex_cbar'], ['-r' num2str(res)], '-a2')  


% set(0,'DefaultFigureColor','remove')
% figure; 
% hf = imagesc(randn(5)); colormap(cm)
% set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
% set(hf, 'visible', 'off')
% cb = colorbar; 
% set(cb, 'fontsize', 30)
% ylabel(cb, unit)
% export_fig([printfolder 'slices_cbar'], ['-r' num2str(res)], '-a2')  

set(0,'DefaultFigureColor',[1 1 1])
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'axial', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_axial'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'sagittal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_sagittal'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'coronal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_coronal'], '-r300', '-a2'); 

% close all











