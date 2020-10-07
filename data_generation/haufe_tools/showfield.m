function h = showfield(z,loc,pars); 
% makes plots of potentials on head
% usage showfield(z,loc,pars);
%
% input:
%  z  values of field/potential to be shown.
%  loc   matrix containing 2D coordinates of channels in second and third column 
% pars is optional
% pars.scale sets the scale of the color map. Either a 1x2 vector
%            corresponding to minimum and a maximum or just a number (say x) 
%            then the scale is from [-x x]. The default is 
%            scale=[ -max(abs(z)) max(abs(z)) ]
% pars.resolution sets the resolution, default is 100
%  pars.cbar=1 draws colorbar, otherwiseno colorbar is not shown. 
%              defaults is 1;
%
% Guido Nolte, 2012-2015
% g.nolte@uke.de

% If you use this code for a publication, please ask Guido Nolte for the
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



ind = find(abs(sum(loc, 2)) > 1e-8); 
z = z(ind, :);
loc = loc(ind, :);

h = 0;

if sign(max(z)) == sign(min(z))
    scal = [min(z) max(z)];
else
    scal=[-max(max(abs(z))),max(max(abs(z)))];
end
cbar=1;
resolution=100;


if nargin>2
    if isfield(pars,'scale')
        scal=pars.scale;
        if length(scal)==1
            scal=[-scal scal];
        end
    end
    if isfield(pars,'cbar')
       cbar=pars.cbar;
    end
    if isfield(pars,'colorbar')
       cbar=pars.colorbar;
    end
    if isfield(pars,'resolution')
       resolution=pars.resolution;
    end
    if ~isfield(pars,'cm')
        cm = jet(64);
    else
        cm=pars.cm;
    end
    if isfield(pars,'marked_electrodes')
        marked_electrodes = pars.marked_electrodes;
        unmarked_electrodes = setdiff(1:size(loc, 1), marked_electrodes);
    else
        unmarked_electrodes = 1:size(loc, 1);
        marked_electrodes = [];
    end
else
    cm = jet(64);
    unmarked_electrodes = 1:size(loc, 1);
    marked_electrodes = [];
end

[n,m]=size(loc);
if m==2;
  x=loc(:,1);
  y=loc(:,2);
else;
  x=loc(:,2);
  y=loc(:,3);
end


xlin = linspace(1.4*min(x),1.4*max(x),resolution);
ylin = linspace(1.4*min(y),1.4*max(y),resolution);
[X,Y] = meshgrid(xlin,ylin);

try
    Z = griddata(x,y,z,X,Y,'invdist');
catch
    Z = griddata(x,y,z,X,Y,'v4');    
end
%Z = griddata(x,y,z,X,Y,'nearest');


  % Take data within head
  rmax=1.02*max(sqrt(x.^2+y.^2));
  mask = (sqrt(X.^2+Y.^2) <= rmax);
  ii = find(mask == 0);
  Z(ii) = NaN;
  
  
surface(X,Y,zeros(size(Z)),Z,'edgecolor','none');shading interp;
colormap(cm);
caxis([ - max(max(abs(z))) max(max(abs(z)))]);
hold on;
plot(x(unmarked_electrodes),y(unmarked_electrodes),'.k','markersize', 5);
plot(x(marked_electrodes),y(marked_electrodes),'ok','markersize', 5, 'linewidth', 2);


%meanx=mean(loc(:,2)*.85+.45);
%meany=mean(loc(:,3)*.85+.45);
scalx=1;
drawhead(0,.0,rmax,scalx);
set(gca,'visible','off');

% axis([-1.2*rmax 1.2*rmax -1.0*rmax 1.4*rmax]);
% axis image;
%axis([-1.4*rmax 1.4*rmax -1.0*rmax 1.4*rmax]);
if cbar==1
  h=colorbar;set(h,'fontweight','bold');
  set(gca, 'linewidth', 2)
  set(h, 'position', [0.85 0.15 0.03 0.7]);
end
set(gca, 'clim', scal);
% set(gca, 'position', [0. 0. 0.75 1]);
axis image


%plot(.985*rmax*sin((0:1000)/1000*2*pi), .985*rmax*cos((0:1000)/1000*2*pi),'linewidth',2,'color','k'); 
return; 



function drawhead(x,y,size,scalx);

cirx=(x+scalx*size*cos((1:1000)*2*pi/1000) )';ciry=(y+size*sin((1:1000)*2*pi/1000))';

plot(cirx,ciry,'k','linewidth', 2);
hold on;

ndiff=20;
plot( [x  cirx(250-ndiff) ],[y+1.1*size ciry(250-ndiff)],'k','linewidth',2);
plot( [x  cirx(250+ndiff) ],[y+1.1*size ciry(250+ndiff)],'k','linewidth',2);


return;
