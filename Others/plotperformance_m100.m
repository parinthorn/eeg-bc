%   This function creates 7 plots relating to performances
%
%   INPUT   inpath_F    =   input estimated F file path
%           inpath_perf =   input performance file path
%   OUTPUT  outpath     =   output file path ('outpath/figx')
%   Written by ANAWAT NARTKULPAT

% function [] = plotperformance(inpath_F,inpath_perf,outpath)
clear
inpath_F = 'data/F_result/F_result_m100';
inpath_perf = 'data/F_result/gcest_performance';
outpath = 'data/F_result';
load(inpath_F);
load(inpath_perf);



%==========================================================================

% plot figure 1: ROC for 20% active with varying number of electrodes
figure(1); hold on;
plot(performance{1,1,1}.FPR_avg,performance{1,1,1}.TPR_avg,'b','LineWidth',1.5);
plot(performance{1,1,2}.FPR_avg,performance{1,1,2}.TPR_avg,'r','LineWidth',1.5);
plot(performance{1,1,3}.FPR_avg,performance{1,1,3}.TPR_avg,'m','LineWidth',1.5);
plot(performance{1,1,4}.FPR_avg,performance{1,1,4}.TPR_avg,'k','LineWidth',1.5);
title('ROC: 20% active sources, 0% deep sources'); xlabel('FPR'); ylabel('TPR');
axis square; axis([0 1 0 1]); 
legend('108 electrodes','61 electrodes','32 electrodes','19 electrodes','Location','southeast');
saveas(gcf,[outpath,'/gcest_20active_roc']); hold off;
print([outpath,'/gcest_20active_roc'], '-depsc')

%==========================================================================

% plot figure 2: ROC for 40% active with varying number of electrodes
figure(2); hold on;
plot(performance{1,2,1}.FPR_avg,performance{1,2,1}.TPR_avg,'b','LineWidth',1.5);
plot(performance{1,2,2}.FPR_avg,performance{1,2,2}.TPR_avg,'r','LineWidth',1.5);
plot(performance{1,2,3}.FPR_avg,performance{1,2,3}.TPR_avg,'m','LineWidth',1.5);
plot(performance{1,2,4}.FPR_avg,performance{1,2,4}.TPR_avg,'k','LineWidth',1.5);
title('ROC: 40% active sources, 0% deep sources'); xlabel('FPR'); ylabel('TPR');
axis square; axis([0 1 0 1]); 
legend('108 electrodes','61 electrodes','32 electrodes','19 electrodes','Location','southeast');
saveas(gcf,[outpath,'/gcest_40active_roc']); hold off;
print([outpath,'/gcest_40active_roc'], '-depsc')

%==========================================================================

% plot figure 3: ROC for 20% active, 61 electrodes with varying shallowness of sources
figure(3); hold on;
plot(performance{1,1,2}.FPR_avg,performance{1,1,2}.TPR_avg,'b','LineWidth',1.5);
plot(performance{2,1,2}.FPR_avg,performance{2,1,2}.TPR_avg,'r','LineWidth',1.5);
plot(performance{3,1,2}.FPR_avg,performance{3,1,2}.TPR_avg,'m','LineWidth',1.5);
title('ROC: 20% active sources, 61 electrodes'); xlabel('FPR'); ylabel('TPR');
axis square; axis([0 1 0 1]); 
legend('0% deep sources','50% deep sources','75% deep sources','Location','southeast');
saveas(gcf,[outpath,'/gcest_20active61electrodes_roc']); hold off;
print([outpath,'/gcest_20active61electrodes_roc'], '-depsc')

%==========================================================================

% plot figure 4: group bar plot of TPR
figure(4);
subplot(121);
thresh_chosen_ind = 2;
x = [1:4];
y = [performance{1,1,1}.TPR_avg(thresh_chosen_ind),performance{2,1,1}.TPR_avg(thresh_chosen_ind),performance{3,1,1}.TPR_avg(thresh_chosen_ind);...
     performance{1,1,2}.TPR_avg(thresh_chosen_ind),performance{2,1,2}.TPR_avg(thresh_chosen_ind),performance{3,1,2}.TPR_avg(thresh_chosen_ind);
     performance{1,1,3}.TPR_avg(thresh_chosen_ind),performance{2,1,3}.TPR_avg(thresh_chosen_ind),performance{3,1,3}.TPR_avg(thresh_chosen_ind);
     performance{1,1,4}.TPR_avg(thresh_chosen_ind),performance{2,1,4}.TPR_avg(thresh_chosen_ind),performance{3,1,4}.TPR_avg(thresh_chosen_ind)];

e = [performance{1,1,1}.TPR_std(thresh_chosen_ind),performance{2,1,1}.TPR_std(thresh_chosen_ind),performance{3,1,1}.TPR_std(thresh_chosen_ind);...
     performance{1,1,2}.TPR_std(thresh_chosen_ind),performance{2,1,2}.TPR_std(thresh_chosen_ind),performance{3,1,2}.TPR_std(thresh_chosen_ind);
     performance{1,1,3}.TPR_std(thresh_chosen_ind),performance{2,1,3}.TPR_std(thresh_chosen_ind),performance{3,1,3}.TPR_std(thresh_chosen_ind);
     performance{1,1,4}.TPR_std(thresh_chosen_ind),performance{2,1,4}.TPR_std(thresh_chosen_ind),performance{3,1,4}.TPR_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0 1.125]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
% legend(b,'3 deep sources','2 deep sources','No deep source');
ylabel('TPR'); xlabel('Number of electrodes'); title('20% active sources');

subplot(122);
y = [performance{1,2,1}.TPR_avg(thresh_chosen_ind),performance{2,2,1}.TPR_avg(thresh_chosen_ind),performance{3,2,1}.TPR_avg(thresh_chosen_ind);...
     performance{1,2,2}.TPR_avg(thresh_chosen_ind),performance{2,2,2}.TPR_avg(thresh_chosen_ind),performance{3,2,2}.TPR_avg(thresh_chosen_ind);
     performance{1,2,3}.TPR_avg(thresh_chosen_ind),performance{2,2,3}.TPR_avg(thresh_chosen_ind),performance{3,2,3}.TPR_avg(thresh_chosen_ind);
     performance{1,2,4}.TPR_avg(thresh_chosen_ind),performance{2,2,4}.TPR_avg(thresh_chosen_ind),performance{3,2,4}.TPR_avg(thresh_chosen_ind)];

e = [performance{1,2,1}.TPR_std(thresh_chosen_ind),performance{2,2,1}.TPR_std(thresh_chosen_ind),performance{3,2,1}.TPR_std(thresh_chosen_ind);...
     performance{1,2,2}.TPR_std(thresh_chosen_ind),performance{2,2,2}.TPR_std(thresh_chosen_ind),performance{3,2,2}.TPR_std(thresh_chosen_ind);
     performance{1,2,3}.TPR_std(thresh_chosen_ind),performance{2,2,3}.TPR_std(thresh_chosen_ind),performance{3,2,3}.TPR_std(thresh_chosen_ind);
     performance{1,2,4}.TPR_std(thresh_chosen_ind),performance{2,2,4}.TPR_std(thresh_chosen_ind),performance{3,2,4}.TPR_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0 1.125]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
legend(b,'0% deep sources','50% deep sources','75% deep sources','Location','northeast');
ylabel('TPR'); xlabel('Number of electrodes'); title('40% active sources');
set(gcf,'Position',[488 342 560 320])
saveas(gcf,[outpath,'/gcest_groupedbar_tpr']); hold off;
print([outpath,'/gcest_groupedbar_tpr'], '-depsc')

%==========================================================================

% plot figure 5: group bar plot of FPR
figure(5);
subplot(121);
x = [1:4];
y = [performance{1,1,1}.FPR_avg(thresh_chosen_ind),performance{2,1,1}.FPR_avg(thresh_chosen_ind),performance{3,1,1}.FPR_avg(thresh_chosen_ind);...
     performance{1,1,2}.FPR_avg(thresh_chosen_ind),performance{2,1,2}.FPR_avg(thresh_chosen_ind),performance{3,1,2}.FPR_avg(thresh_chosen_ind);
     performance{1,1,3}.FPR_avg(thresh_chosen_ind),performance{2,1,3}.FPR_avg(thresh_chosen_ind),performance{3,1,3}.FPR_avg(thresh_chosen_ind);
     performance{1,1,4}.FPR_avg(thresh_chosen_ind),performance{2,1,4}.FPR_avg(thresh_chosen_ind),performance{3,1,4}.FPR_avg(thresh_chosen_ind)];

e = [performance{1,1,1}.FPR_std(thresh_chosen_ind),performance{2,1,1}.FPR_std(thresh_chosen_ind),performance{3,1,1}.FPR_std(thresh_chosen_ind);...
     performance{1,1,2}.FPR_std(thresh_chosen_ind),performance{2,1,2}.FPR_std(thresh_chosen_ind),performance{3,1,2}.FPR_std(thresh_chosen_ind);
     performance{1,1,3}.FPR_std(thresh_chosen_ind),performance{2,1,3}.FPR_std(thresh_chosen_ind),performance{3,1,3}.FPR_std(thresh_chosen_ind);
     performance{1,1,4}.FPR_std(thresh_chosen_ind),performance{2,1,4}.FPR_std(thresh_chosen_ind),performance{3,1,4}.FPR_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0 0.12]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
% legend(b,'3 deep sources','2 deep sources','No deep source');
ylabel('FPR'); xlabel('Number of electrodes'); title('20% active sources');
legend(b,'0% deep sources','50% deep sources','75% deep sources','Location','northeast');

subplot(122);
y = [performance{1,2,1}.FPR_avg(thresh_chosen_ind),performance{2,2,1}.FPR_avg(thresh_chosen_ind),performance{3,2,1}.FPR_avg(thresh_chosen_ind);...
     performance{1,2,2}.FPR_avg(thresh_chosen_ind),performance{2,2,2}.FPR_avg(thresh_chosen_ind),performance{3,2,2}.FPR_avg(thresh_chosen_ind);
     performance{1,2,3}.FPR_avg(thresh_chosen_ind),performance{2,2,3}.FPR_avg(thresh_chosen_ind),performance{3,2,3}.FPR_avg(thresh_chosen_ind);
     performance{1,2,4}.FPR_avg(thresh_chosen_ind),performance{2,2,4}.FPR_avg(thresh_chosen_ind),performance{3,2,4}.FPR_avg(thresh_chosen_ind)];

e = [performance{1,2,1}.FPR_std(thresh_chosen_ind),performance{2,2,1}.FPR_std(thresh_chosen_ind),performance{3,2,1}.FPR_std(thresh_chosen_ind);...
     performance{1,2,2}.FPR_std(thresh_chosen_ind),performance{2,2,2}.FPR_std(thresh_chosen_ind),performance{3,2,2}.FPR_std(thresh_chosen_ind);
     performance{1,2,3}.FPR_std(thresh_chosen_ind),performance{2,2,3}.FPR_std(thresh_chosen_ind),performance{3,2,3}.FPR_std(thresh_chosen_ind);
     performance{1,2,4}.FPR_std(thresh_chosen_ind),performance{2,2,4}.FPR_std(thresh_chosen_ind),performance{3,2,4}.FPR_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0 0.12]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
ylabel('FPR'); xlabel('Number of electrodes'); title('40% active sources');
set(gcf,'Position',[488 342 560 320])
saveas(gcf,[outpath,'/gcest_groupedbar_fpr']); hold off;
print([outpath,'/gcest_groupedbar_fpr'], '-depsc')

%==========================================================================

% plot figure 6: group bar plot of ACC
figure(6);
subplot(121);
x = [1:4];
y = [performance{1,1,1}.ACC_avg(thresh_chosen_ind),performance{2,1,1}.ACC_avg(thresh_chosen_ind),performance{3,1,1}.ACC_avg(thresh_chosen_ind);...
     performance{1,1,2}.ACC_avg(thresh_chosen_ind),performance{2,1,2}.ACC_avg(thresh_chosen_ind),performance{3,1,2}.ACC_avg(thresh_chosen_ind);
     performance{1,1,3}.ACC_avg(thresh_chosen_ind),performance{2,1,3}.ACC_avg(thresh_chosen_ind),performance{3,1,3}.ACC_avg(thresh_chosen_ind);
     performance{1,1,4}.ACC_avg(thresh_chosen_ind),performance{2,1,4}.ACC_avg(thresh_chosen_ind),performance{3,1,4}.ACC_avg(thresh_chosen_ind)];

e = [performance{1,1,1}.ACC_std(thresh_chosen_ind),performance{2,1,1}.ACC_std(thresh_chosen_ind),performance{3,1,1}.ACC_std(thresh_chosen_ind);...
     performance{1,1,2}.ACC_std(thresh_chosen_ind),performance{2,1,2}.ACC_std(thresh_chosen_ind),performance{3,1,2}.ACC_std(thresh_chosen_ind);
     performance{1,1,3}.ACC_std(thresh_chosen_ind),performance{2,1,3}.ACC_std(thresh_chosen_ind),performance{3,1,3}.ACC_std(thresh_chosen_ind);
     performance{1,1,4}.ACC_std(thresh_chosen_ind),performance{2,1,4}.ACC_std(thresh_chosen_ind),performance{3,1,4}.ACC_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0.85 1]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
% legend(b,'3 deep sources','2 deep sources','No deep source');
ylabel('ACC'); xlabel('Number of electrodes'); title('20% active sources');

subplot(122);
y = [performance{1,2,1}.ACC_avg(thresh_chosen_ind),performance{2,2,1}.ACC_avg(thresh_chosen_ind),performance{3,2,1}.ACC_avg(thresh_chosen_ind);...
     performance{1,2,2}.ACC_avg(thresh_chosen_ind),performance{2,2,2}.ACC_avg(thresh_chosen_ind),performance{3,2,2}.ACC_avg(thresh_chosen_ind);
     performance{1,2,3}.ACC_avg(thresh_chosen_ind),performance{2,2,3}.ACC_avg(thresh_chosen_ind),performance{3,2,3}.ACC_avg(thresh_chosen_ind);
     performance{1,2,4}.ACC_avg(thresh_chosen_ind),performance{2,2,4}.ACC_avg(thresh_chosen_ind),performance{3,2,4}.ACC_avg(thresh_chosen_ind)];

e = [performance{1,2,1}.ACC_std(thresh_chosen_ind),performance{2,2,1}.ACC_std(thresh_chosen_ind),performance{3,2,1}.ACC_std(thresh_chosen_ind);...
     performance{1,2,2}.ACC_std(thresh_chosen_ind),performance{2,2,2}.ACC_std(thresh_chosen_ind),performance{3,2,2}.ACC_std(thresh_chosen_ind);
     performance{1,2,3}.ACC_std(thresh_chosen_ind),performance{2,2,3}.ACC_std(thresh_chosen_ind),performance{3,2,3}.ACC_std(thresh_chosen_ind);
     performance{1,2,4}.ACC_std(thresh_chosen_ind),performance{2,2,4}.ACC_std(thresh_chosen_ind),performance{3,2,4}.ACC_std(thresh_chosen_ind)];

b = bar(x,y); hold on; xlim([0.25 4.75]); ylim([0.85 1]); axis square;
set(gca, 'XTickLabel', {'108','61','32','19'});
groupwidth = 3/4.5;
for i = 1:3
    xx = (1:4) - groupwidth/2 + (2*i-1)*groupwidth/6;
    er = errorbar(xx,y(:,i),e(:,i),'.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
legend(b,'0% deep sources','50% deep sources','75% deep sources','Location','northeast');
ylabel('ACC'); xlabel('Number of electrodes'); title('40% active sources');
set(gcf,'Position',[488 342 560 320])
% ----- fixed legend location
saveas(gcf,[outpath,'/gcest_groupedbar_acc']); hold off;
print([outpath,'/gcest_groupedbar_acc'], '-depsc')

%==========================================================================

% plot figure 7: example of ground truth F with estimated F from varying number of electrodes
figure(7); 
[n_model,n_realization] = size(M.F{1,1,1});
ground_truth_ind = 12;
realization_ind = 1;
ax1 = subplot(131); imagesc(M.F_true{2,1}{ground_truth_ind}); colormap(flipud(hot)); axis square; title('Ground truth GC');
% caxis([0 1]);
% colorbar;
TPR = performance{2,1,2}.TPR(ground_truth_ind,thresh_chosen_ind);
FPR = performance{2,1,2}.FPR(ground_truth_ind,thresh_chosen_ind);
ACC = performance{2,1,2}.ACC(ground_truth_ind,thresh_chosen_ind);
ax2 = subplot(132); imagesc(M.F{2,1,2}{ground_truth_ind,realization_ind}); colormap(flipud(hot)); axis square;
title('Estimated GC, 61 electrodes');% caxis([0 1]);
xlabel(sprintf('TPR=%.3f, FPR=%.3f, ACC=%.3f',TPR,FPR,ACC),'Interpreter','None');
% colorbar;
TPR = performance{2,1,4}.TPR(ground_truth_ind,thresh_chosen_ind);
FPR = performance{2,1,4}.FPR(ground_truth_ind,thresh_chosen_ind);
ACC = performance{2,1,4}.ACC(ground_truth_ind,thresh_chosen_ind);
ax3 = subplot(133); imagesc(M.F{2,1,4}{ground_truth_ind,realization_ind}); colormap(flipud(hot)); axis square;
title('Estimated GC, 19 electrodes');% caxis([0 1]);
xlabel(sprintf('TPR=%.3f, FPR=%.3f, ACC=%.3f',TPR,FPR,ACC),'Interpreter','None');
% colorbar;
set(gcf,'Position',[198.6000 451.4000 962.4000 277.6000]);

linkaxes([ax1,ax2,ax3],'y');
 ha=get(gcf,'children');
set(findobj(ha,'type','line'),'linew',1.4)
 set(ha(1),'position',[.65 .11 .2 .815])
 set(ha(2),'position',[.4 .11 .2 .815])
 set(ha(3),'position',[.15 .11 .2 .815])
 set(ha(2),'Yticklabel',[])
 set(ha(1),'Yticklabel',[])

saveas(gcf,[outpath,'/gcest_fhat61and19electrode']); hold off;
print([outpath,'/gcest_fhat61and19electrode'], '-depsc')




% end