function allplots_cortex(sa, headdata, colorlimits, cm, unit, printfolder, varargin)
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
mkdir(printfolder)

res = '150';


% figure; showfield(data, sa.locs_2D, struct('scale', colorlimits, 'cbar', 0, 'cm', cm));
% export_fig([printfolder '2Dtopo'], ['-r' num2str(res)], '-a2'); 
% 


surface_pars = struct('myviewdir', [-1 0 0], 'alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []);


figure; showsurface3(sa.head.vc, sa.head.tri, surface_pars, headdata, varargin{:});
export_fig([printfolder 'scalp_left'], ['-r' num2str(res)], '-a2'); 

surface_pars = struct('myviewdir', [1 0 0], 'alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []);

figure; showsurface3(sa.head.vc, sa.head.tri, surface_pars, headdata, varargin{:});
export_fig([printfolder 'scalp_right'], ['-r' num2str(res)], '-a2'); 


surface_pars = struct('myviewdir', [0 -1e-10 1], 'alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []); %[180 90]

figure; showsurface3(sa.head.vc, sa.head.tri, surface_pars, headdata, varargin{:});
export_fig([printfolder 'scalp_top'], ['-r' num2str(res)], '-a2'); 



% sa.head.vc(sa.head.vc(:, 3) > -30, 3) = -30;
% 
% 
% surface_pars = struct('myviewdir', [-1 -2 3], 'directions', [1 1 1 1 0 0], 'alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []); %[180 90]
% 
% sa.head.vc(sa.head.vc(:, 3) > -30, 3) = -30;
% figure; showsurface3(sa.head.vc, sa.head.tri, surface_pars, headdata, varargin{:});
% export_fig([printfolder 'scalp_angular_left'], ['-r' num2str(res)], '-a2'); 
% 
% 
% surface_pars = struct('myviewdir', [1 -2 3], 'directions', [1 1 1 1 0 0], 'alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', []); %[180 90]
% 
% figure; showsurface3(sa.head.vc, sa.head.tri, surface_pars, headdata, varargin{:});
% export_fig([printfolder 'scalp_angular_right'], ['-r' num2str(res)], '-a2'); 


figure; 
hf = imagesc(randn(5)); colormap(cm)
set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
set(hf, 'visible', 'off')
cb = colorbar; 
set(cb, 'fontsize', 30)
ylabel(cb, unit)
export_fig([printfolder 'scalp_cbar'], ['-r' num2str(res)], '-a2')  


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

% close all











