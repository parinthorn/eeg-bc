%% raw accuracy: full compare mhat = 100
clear
clc
load('result_C_augment_L_100realization_with_bic_index.mat')



RR = compare_C(M.C_nz_ind{1},M.C_true_ind{1},100,75);
FPR_raw{1,1,1} = RR.FPR_all;
TPR_raw{1,1,1} = RR.TPR_all;
ACC_raw{1,1,1} = RR.ACC_all;
% FPR_array = mean(mean(RR.FPR_all,1),3);
% TPR_array = mean(mean(RR.TPR_all,1),3);
% ACC_array = mean(mean(RR.ACC_all,1),3);
clear RR
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\result_C_augment_L_100realization_raw_accuracy_fullcompare','RR')
%% raw accuracy: only groundtruth(m=50) compare mhat = 100
clear
clc
load('result_C_augment_L_100realization_with_bic_index.mat')
for ii=1:100
    for jj=1:75
        tmp = M.C_nz_ind{1,1,1}{ii}{jj};
        tmp(tmp>50) = [];
        M.C_nz_ind{1,1,1}{ii}{jj} = tmp;
    end

end
RR = compare_C(M.C_nz_ind{1},M.C_true_ind{1},50,75);
FPR_raw{1,1,1} = RR.FPR_all;
TPR_raw{1,1,1} = RR.TPR_all;
ACC_raw{1,1,1} = RR.ACC_all;
clear RR
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\result_C_augment_L_100realization_raw_accuracy_partial_compare','RR')
%% raw accuracy: only groundtruth(m=50) compare mhat = 50
clear
clc
% load('result_C_augment_L_100realization_with_bic_index.mat')
load('result_C_augment_L_m50_100realization_with_bic_index')
RR = compare_C(M.C_nz_ind{1},M.C_true_ind{1},50,75);
FPR_raw{1,1,1} = RR.FPR_all;
TPR_raw{1,1,1} = RR.TPR_all;
ACC_raw{1,1,1} = RR.ACC_all;
clear RR
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\result_C_augment_L_m50_100realization_raw_accuracy','RR')
%%

figure(1)

% clear
load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\C_result_augment_L_m50.mat')
Lname = {'L10'};
nbatch = 1;
n_model = 100;
n_deep = [2];
density = [10];

density_it = 0;
for n_act = density
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for LL=Lname
            n_elec = n_elec+1;
            RR = compare_C(M.C_nz_ind{deep_it,density_it,n_elec},M.C_true_ind{deep_it,density_it},50,75);
            FPR_array_2{deep_it,density_it,n_elec} = mean(mean(RR.FPR_all,1),3);
            TPR_array_2{deep_it,density_it,n_elec} = mean(mean(RR.TPR_all,1),3);
            ACC_array_2{deep_it,density_it,n_elec} = mean(mean(RR.ACC_all,1),3);
        end
    end
end

%%
plot(FPR_array_2{1},TPR_array_2{1},'color',[0 0.5 0],'LineWidth',1.4)
hold on
plot(FPR_array,TPR_array,'b','LineWidth',1.4)
hold off

axis([0 1 0 1])
axis('square')
box on
grid on
xlabel('FPR')
ylabel('TPR')
legend('$$\mathrm{\tilde{m}}$$ = 50','$$\mathrm{\tilde{m}}$$ = 100','Interpreter','latex','Location','southeast')

set(gca,'FontSize',18)

%% result C compare
clear
clc
load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\C_result_L_augment.mat')
M_100 = M;
load('G:\Shared drives\MASTER_DRIVE\Journal\result_C_augment_L.mat')
M_50 = M;
clear M
%%
mm=6;
figure(1)
subplot(1,3,1)
imagesc(((abs(M_100.C_true{1}{mm}))))
title('ground-truth $$C$$','Interpreter','latex')
subplot(1,3,2)
imagesc((abs(M_50.C_bic{1}{mm})))

title('$$\hat{C}_{\mathrm{bic}}$$, $$\tilde{m}=50$$','Interpreter','latex')
subplot(1,3,3)
imagesc((abs(M_100.C_bic{1}{mm})))


title('$$\hat{C}_{\mathrm{bic}}$$, $$\tilde{m}=100$$','Interpreter','latex')
ax1 = subplot(1,3,1);
colormap(flipud(hot))
ylabel('$$\mathrm{n^{th} ~source}$$','Interpreter','latex')
grid on
box on

ax2 = subplot(1,3,2);
colormap(flipud(hot))
xlabel('$$\mathrm{n^{th} ~state}$$','Interpreter','latex')
grid on
box on
ax3 = subplot(1,3,3);
colormap(flipud(hot))
grid on
box on

linkaxes([ax1,ax2,ax3],'y');
 ha=get(gcf,'children');
set(findobj(ha,'type','line'),'linew',1.4)
 set(ha(1),'position',[.6 .11 .2 .815],'FontSize',24)
 set(ha(2),'position',[.4 .11 .2 .815],'FontSize',24)
 set(ha(3),'position',[.18 .11 .2 .815],'FontSize',24)
 set(ha(2),'Yticklabel',[])
 set(ha(1),'Yticklabel',[])


%% bic boxplot

RR = compare_C(M_100.C_nz_ind_bic{1},M_100.C_true_ind{1},100,-1);
FPR_array = RR.FPR_all;
TPR_array = RR.TPR_all;
ACC_array = RR.ACC_all;

%% boxplot
clf
Lname = {'L10'};
nbatch = 5;
n_model = 20;
n_deep = [0 2 3];
DD = {'50%'};
density = [10];

density_it = 0;
cnt = 0;
ind_plot = [3 1;4 2];
% 221 ->set(gca, 'Position',[0.1300 0.55 0.3347 0.3370])
% 224 ->p = get(gca,'Position'); p(2) = 0.55;

for n_act = density
    density_it = density_it+1;
    deep_it = 0;

    for nd = n_deep
        deep_it = deep_it+1;
        figure(deep_it)
%         clf;
        for kk=1:1
            FPR_elec(:,kk) = FPR_array;
            TPR_elec(:,kk) = TPR_array;
            ACC_elec(:,kk) = ACC_array;
        end
        subplot(2,2,ind_plot(density_it,1))

        boxplot(FPR_elec);
        set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',18)

%         set(gca, 'Position',[0.1300 0.55 0.3347 0.3370])
        xlabel('Number of electrodes')
        if density_it ==1
            ylabel('FPR')
        end
%         title(['m_active ',int2str(density(density_it)*2),'%'],'Interpreter','None')
        ylim([0 0.2])
        axis('normal')
        grid on


        subplot(2,2,ind_plot(density_it,2))
        boxplot(TPR_elec);
        set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',18)

%         xlabel('Number of electrodes')
        if density_it ==1
            ylabel('TPR')
        end

        title([int2str(density(density_it)*2),'% active sources'],'Interpreter','None')
        ylim([0 1])
        axis('normal')
        grid on


%         subplot(1,3,3)
%         boxplot(ACC_elec);
%         set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',16)
%         xlabel('Number of electrodes')
%         ylabel('TPR (100 models)')
% %         title(['m_active ',int2str(density(density_it)*2),'%'],'Interpreter','None')
%         ylim([0 1])
%         axis('normal')
%         grid on

%         sgtitle(['m_active ',int2str(density(density_it)*2),'%',', ',DD{deep_it},' deep sources'],'Interpreter','None')
    end

end
HH = {'0','50','75'};
for ii=1:3
        figure(ii)
 ha=get(gcf,'children');
set(findobj(ha,'type','line'),'linew',1.4)
 set(ha(2),'position',[.5 .08 .4 .4],'FontSize',18)
 set(ha(4),'position',[.1 .08 .4 .4],'FontSize',18)
 set(ha(1),'position',[.5 .5 .4 .4],'FontSize',18)
 set(ha(3),'position',[.1 .5 .4 .4],'FontSize',18)
 set(ha(1),'Xticklabel',[],'Yticklabel',[])
 set(ha(3),'Xticklabel',[])
 set(ha(2),'Yticklabel',[])
%  set(gcf, 'Position', get(0, 'Screensize'));
%  print(figure(ii),['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_bic_deep_',HH{ii}],'-depsc')
end
% close all
% subplot(2,2,1,'Position',[ 0 0 0])
