
% this script visualize the result of ss method vs brainstorm method vs groundtruth on the 75K NY head brain mesh

clc
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('G:\Shared drives\MASTER_DRIVE\Journal\DATA_GENERATION\BBCB_code\data\sa')
end
load('G:\Shared drives\MASTER_DRIVE\Journal\Source_Localization\NYHead\tess_cortex_01.mat','Vertices')
%%
model_loc = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';
mets = {'WMNE','LCMV','sLORETA'};
map_loc = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\map_reconstruct_75K_fixtime_no_noise_cov\';
load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\C_result_L_augment.mat')
load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\source_amplitude_augment','s_amp_augmented')

clf
close all

load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\blue_to_red_colorspace')
load('F_result_roi_permind')
% color_space = blue_to_red_colorspace;
% color_space = hot;
color_space = jet;
%%
for iter = 1:100
% iter = 7;
angle = {[0,90],[90,0],[-90,0],[180,-90]};
thresh = 0.4;
thresh_est = 0.4;
view_angle = {'top','right','left','bottom'};
load([map_loc,'map_',int2str(iter)])
load([model_loc,'model__act_10_deep_2_',int2str(iter),'.mat'])
active_index_truth = model.PARAMETER.ind_active;
truth = model.truth;
% clear model
cnt = 0;
plot_sequence = [1,2,3,4];
t = tiledlayout(4,5,'TileSpacing','None','Padding','None');
for vv=[1 2 3 4]
    cnt = cnt+1;
    %plot ground truth
    nexttile
    vX = sa.cortex75K.vc(:,1);
    vY = sa.cortex75K.vc(:,2);
    vZ = sa.cortex75K.vc(:,3);
    muX = mean(vX);
    muY = mean(vY);
    muZ = mean(vZ);
    vX = vX-muX;
    vY = vY-muY;
    vZ = vZ-muZ;
    val = zeros(50,1);
    val(active_index_truth) = 1;
    Source2K = (truth.source_amp*val);
    Source75K = Source2K(sa.cortex2K.in_to_cortex75K_eucl);
    gt_val = norms(model.source_model0.C,2,2);
    gt_val = gt_val/max(gt_val);
    true_active = find(gt_val);

    gt_base = zeros(size(vX,1),1);
    scatter3(vX,vY,vZ,[],gt_base,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor',[0 0 1]) % base ground truth
    hold on


    vX_highlight = vX(Source75K>=thresh);
    vY_highlight = vY(Source75K>=thresh);
    vZ_highlight = vZ(Source75K>=thresh);
    vX_gtcenter = model.truth.centers(:,1)-muX;
    vY_gtcenter = model.truth.centers(:,2)-muY;
    vZ_gtcenter = model.truth.centers(:,3)-muZ;
    gt_colormap = Source75K(Source75K>=thresh);
    scatter3(vX_highlight,vY_highlight,vZ_highlight,[],gt_colormap,'filled','MarkerFaceColor',[1 0 0]) % ground-truth source amplitude
% scatter3(vX_gtcenter(true_active),vY_gtcenter(true_active),vZ_gtcenter(true_active),80*ones(size(vX_gtcenter(true_active),1),1),gt_val(true_active),'filled')

    hold off
    % axis('square')
    axis([min(vX) max(vX) min(vY) max(vY) min(vZ) max(vZ)]*1.2)
    % axis([])
    daspect([1 1 1])

    box off
    grid off
    axis off
    set(gca,'Xticklabel',[],'Xtick',[])
    set(gca,'Yticklabel',[],'Ytick',[])
    set(gca,'Zticklabel',[],'Ztick',[])
    colormap(color_space)
    view(angle{vv})
    if vv==1
        title('Ground truth','FontSize',18)
    end


    %=======================
    nexttile


    active_index_estimate = M.C_nz_ind_bic{1}{iter};
    BIC_val = norms(M.C_bic{1}{iter},2,2);
    BIC_val = full(BIC_val/max(BIC_val));
%     inactive_index_estimate = setdiff(1:1:100,active_index_estimate);
    s_amp = s_amp_augmented{iter};
    val = zeros(100,1);
    val(active_index_estimate) = 1; %BIC_val(active_index_estimate);
    val_inact = 1-val;

    Source2K_inactive = (s_amp*val_inact);
    Source75K_inactive = Source2K_inactive(sa.cortex2K.in_to_cortex75K_eucl);

    Source2K = (s_amp*val);
    Source75K = Source2K(sa.cortex2K.in_to_cortex75K_eucl);

    tmpvar = (s_amp*ones(100,1));
    selected_source = tmpvar(sa.cortex2K.in_to_cortex75K_eucl);
    scatter3(vX,vY,vZ,[],zeros(size(vX,1),1),'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor',[0 0 1]) % base
    hold on
%     vX_selected = vX(selected_source>=thresh);
%     vY_selected = vY(selected_source>=thresh);
%     vZ_selected = vZ(selected_source>=thresh);

    vX_highlight = vX(Source75K>=thresh_est);
    vY_highlight = vY(Source75K>=thresh_est);
    vZ_highlight = vZ(Source75K>=thresh_est);
    inact_thresh = 0.7;
    vX_inact = vX(Source75K_inactive>=inact_thresh);
    vY_inact = vY(Source75K_inactive>=inact_thresh);
    vZ_inact = vZ(Source75K_inactive>=inact_thresh);

    estimated_index = ((selected_source>=thresh));

    % FN_index = setdiff(find(true_index),find(estimated_index));
    % FP_index = setdiff(find(estimated_index),find(true_index));
%     all_source_colormap = zeros(size(vX_selected,1),1);
%     scatter3(vX_selected,vY_selected,vZ_selected,[],all_source_colormap,'filled','MarkerFaceAlpha',0.16,'MarkerFaceColor',[1 1 1]) % all sources in consideration

    % roi_tresh = ((sa.cortex75K.roi_mask==5)|(sa.cortex75K.roi_mask==6)|(sa.cortex75K.roi_mask==7)|(sa.cortex75K.roi_mask==8)); % revise
    % scatter3(vX(roi_tresh),vY(roi_tresh),vZ(roi_tresh),[],roi_tresh(roi_tresh==1),'filled','MarkerFaceAlpha',0.16,'MarkerFaceColor',[1 1 1])
    est_colormap = Source75K(Source75K>=thresh_est)*0;
    est_colormap_inactive = Source75K_inactive(Source75K_inactive>=inact_thresh)*0;
    centers = [model.truth.centers;model.truth.augment_data.centers];
    vX_inact = centers(find(val_inact),1)-muX;
    vY_inact = centers(find(val_inact),2)-muY;
    vZ_inact = centers(find(val_inact),3)-muZ;
    vX_act = centers(find(val),1)-muX;
    vY_act = centers(find(val),2)-muY;
    vZ_act = centers(find(val),3)-muZ;
    color_inact = zeros(90,1);
    MarkerSize = 80*ones(90,1);
%     scatter3(vX_inact,vY_inact,vZ_inact,[],est_colormap_inactive,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',[1 1 1])
%     scatter3(vX_inact,vY_inact,vZ_inact,MarkerSize,color_inact,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',[1 1 1])
%     scatter3(vX_inact,vY_inact,vZ_inact,[],est_colormap_inactive,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',[1 1 1]) % estimated active source location
    figobj = scatter3(vX_highlight,vY_highlight,vZ_highlight,[],est_colormap,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',[1 0 0]); % estimated active source location
%     colormap(figobj,autumn(5))
%     scatter3(vX_act,vY_act,vZ_act,80*ones(size(vX_act,1),1),BIC_val(find(val)),'filled','MarkerFaceAlpha',1) % estimated active source location
    % axis('square')
    axis([min(vX) max(vX) min(vY) max(vY) min(vZ) max(vZ)]*1.2)
    daspect([1 1 1])
    box off
    grid off
    axis off
    set(gca,'Xticklabel',[],'Xtick',[])
    set(gca,'Yticklabel',[],'Ytick',[])
    set(gca,'Zticklabel',[],'Ztick',[])
    view(angle{vv})
    if vv==1
        title('Proposed method','FontSize',18)
    end
    colormap(color_space)
    %=============================

    imtn = zeros(size(Vertices,1),3);
    imtn(:,1) = x.MN;
    imtn(:,2) = x.LCMV;
    imtn(:,3) = x.sLORETA;

    % Normalize/rescale imtn (from min-max to 0-1), use for coloring
    [nvert, nmet] = size(imtn);
    normalized_imtn = zeros(size(imtn));

    for imet = 1:nmet
        im_im = imtn(:,imet);
        normalized_imtn(:,imet) = (im_im - min(im_im))/(max(im_im) - min(im_im));
    end

    % Plotting scatter plot
    vX = Vertices(:,1)*1000; % to mm
    vY = Vertices(:,2)*1000;
    vZ = Vertices(:,3)*1000;
    vX = vX-mean(vX);
    vY = vY-mean(vY);
    vZ = vZ-mean(vZ);
    figure(1)
    tmpXY = ([0 -1;1 0;]*[vX';vY'])'; % rotate to same coordinate
    vX = tmpXY(:,1);
    vY = tmpXY(:,2);
    for imet = 1:nmet
        nexttile
        %     sp(imet) = subplot(4,5,(cnt-1)*5+plot_sequence(imet));
        scatter3(vX,vY,vZ,[],normalized_imtn(:,imet),'filled')
        set(gca,'Xticklabel',[],'Xtick',[])
        set(gca,'Yticklabel',[],'Ytick',[])
        set(gca,'Zticklabel',[],'Ztick',[])

        view(angle{vv})
        ti = strcat(mets{imet});       % Ex: ti = 'LCMV Method'
        if vv==1
            title(ti,'FontSize',18)
        end
        %     axis('square')
        axis([min(vX) max(vX) min(vY) max(vY) min(vZ) max(vZ)]*1.2)
        daspect([1 1 1])
        box off
        grid off
        axis off
        colormap(color_space)
    end


    % print(figure(vv),['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\','source_reconstruct_solid_',view_angle{vv}],'-depsc')
    % saveas(gcf,['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\','source_reconstruct_solid_',view_angle{vv}],'fig')
end
pause(0.1)
print(figure(1),['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_compare_brainstorm_',int2str(iter)],'-depsc')
end

% print(figure(1),'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_compare_brainstorm','-depsc')
% saveas(gca,'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_compare_brainstorm','fig')
%%

permute_roi_index =permind{iter};

groundtruth_source_cluster_index = cell(model.PARAMETER.n_source_cluster,1);
for ii=1:model.PARAMETER.n_source_cluster-1
    groundtruth_source_cluster_index{ii} = [model.PARAMETER.ind_cluster_source(ii):1:model.PARAMETER.ind_cluster_source(ii+1)-1];
end
groundtruth_source_cluster_index{end} = model.PARAMETER.ind_cluster_source(end):1:50;
tmpf = @(x) x+50;
source_cluster_index = [groundtruth_source_cluster_index; ...
    cellfun(tmpf,groundtruth_source_cluster_index, ...
    'UniformOutput',false)];
permute_source_index = (source_cluster_index(permute_roi_index));
centers_reindex = centers([permute_source_index{:}],:);
H = combnk(1:100,2);
dist_mat=zeros(100,100);
for kk = 1:size(H,1)
    ii = H(kk,1);
    jj = H(kk,2);
    dist_mat(ii,jj) = norm(centers_reindex(ii,:)-centers_reindex(jj,:),2);
end
figure(2);
imagesc(dist_mat');
load('miscdata','rois')
roi_name = [model.truth.rois model.truth.augment_data.rois];
roi_name = roi_name([permute_source_index{:}]);
% roi_name = {roi_name{1}};
set(gca,'ytick',1:1:100,'yticklabel',roi_name)
set(gca,'xtick',1:1:100,'xticklabel',roi_name)
xtickangle(-90)
axis('square')
set(gca,'xaxisLocation','top')
% set(gca,'ytick',1:1:8,'yticklabel',rois)
