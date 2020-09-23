%%
load('G:\Shared drives\MASTER_DRIVE\Journal\Source_Localization\NYHead\tess_cortex_01.mat','Vertices')
%% downsample
dsFactor = 2000/(74661); %  sampling 75K to 2K source
[NewTessMat.Faces, NewTessMat.Vertices] = reducepatch(TessMat.Faces, TessMat.Vertices, dsFactor);
[~, I, J] = intersect(TessMat.Vertices, NewTessMat.Vertices, 'rows');
[I, iSort] = sort(I);
% row_ind = 1:1:74661;
row_ind = I;

%%


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
save('Experiment\Kernel_roi2K','Kernel_roi')
