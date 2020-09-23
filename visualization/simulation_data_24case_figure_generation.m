%% plot ROC C
clear
clc
% load('result_C_100realization_with_bic_index')
load('debug_result_C_100realization_with_bic_index.mat')
%% ROC

Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
% n_deep = [0 2 3];
n_deep = [2];
density = [10 20];

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
%             FPR_array{deep_it,density_it,n_elec} = mean(mean(RR.FPR_all,1),3);
%             TPR_array{deep_it,density_it,n_elec} = mean(mean(RR.TPR_all,1),3);
%             ACC_array{deep_it,density_it,n_elec} = mean(mean(RR.ACC_all,1),3);
            FPR_raw{deep_it,density_it,n_elec} = RR.FPR_all;
            TPR_raw{deep_it,density_it,n_elec} = RR.TPR_all;
            ACC_raw{deep_it,density_it,n_elec} = RR.ACC_all;
        end
    end
end
clear RR
% RR.FPR_array = FPR_array;
% RR.TPR_array = TPR_array;
% RR.ACC_array = ACC_array;

RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
% save('result_C_100realization_raw_accuracy','RR')
save('debug_result_C_100realization_raw_accuracy_onlydeep2','RR')


%% fig 1 fig2
% load('result_C_m100_accuracy')
load('debug_result_C_100realization_raw_accuracy_onlydeep2')
FPR_array = RR.FPR;
TPR_array = RR.TPR;
ACC_array = RR.ACC;
% n_deep = {'deep_0','deep_50','deep_75'};
n_deep = {'deep_50'};
c = [1 0 0;0 0.5 0;0 0 1;0 0 0];
TT = [20,40];
for density=1:size(FPR_array,2)
    for nd =1:size(FPR_array,1)
        figure(1)

        for n_elec=1:4
%             plot(FPR_array{nd,density,n_elec},TPR_array{nd,density,n_elec},'LineWidth',1.2,'Color',c(n_elec,:))
            plot(mean(FPR_array{nd,density,n_elec},1),mean(TPR_array{nd,density,n_elec},1),'LineWidth',1.2,'Color',c(n_elec,:))
            if n_elec==1
                hold on
            end

        end
        legend('108 electrodes','61 electrodes','32 electrodes','19 electrodes','Location','southeast')
        box on
        grid on
        hold off
        axis([0 1 0 1])
        axis('square')
        xlabel('FPR')
        ylabel('TPR')
        set(gca,'FontSize',18)

%         fig_dir = 'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\';
        % mkdir(fig_dir)
%         print(figure(1),[fig_dir,'source_selection_active',int2str(TT(density)),'_',(n_deep{nd})],'-depsc')
%         saveas(gcf,[fig_dir,'source_selection_active',int2str(TT(density)),'_',(n_deep{nd})],'fig')
    end
end

%% fig 3

FPR_array = RR.FPR_array;
TPR_array = RR.TPR_array;
ACC_array = RR.ACC_array;
figure(1)
c = [1 0 0;0 0.5 0;0 0 1;0 0 0];
% n_deep = 3;
n_elec = 2;
density = 1;

for nd=1:3
    plot(FPR_array{nd,density,n_elec},TPR_array{nd,density,n_elec},'LineWidth',1.2,'Color',c(nd,:))
    if nd==1
        hold on
    end

end
legend('0% deep source','50% deep sources','75% deep sources','Location','southeast')
box on
grid on
hold off
axis([0 1 0 1])
axis('square')
xlabel('FPR')
ylabel('TPR')
set(gca,'FontSize',14)

fig_dir = 'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\';
% mkdir(fig_dir)
print(figure(1),[fig_dir,'source_selection_active10_nelec_61'],'-depsc')
saveas(gcf,[fig_dir,'source_selection_active10_nelec_61'],'fig')

%% bic boxplot data
load('result_C_m100')
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
        for LL=Lname
            n_elec = n_elec+1;
            RR = compare_C(M.C_nz_ind_bic{deep_it,density_it,n_elec},M.C_true_ind{deep_it,density_it},50,-1);
            E.FPR_array{deep_it,density_it,n_elec} = RR.FPR_all;
            E.TPR_array{deep_it,density_it,n_elec} = RR.TPR_all;
            E.ACC_array{deep_it,density_it,n_elec} = RR.ACC_all;
            E.FPR_array_mean{deep_it,density_it,n_elec} = mean(RR.FPR_all);
            E.TPR_array_mean{deep_it,density_it,n_elec} = mean(RR.TPR_all);
            E.ACC_array_mean{deep_it,density_it,n_elec} = mean(RR.ACC_all);
            E.FPR_array_std{deep_it,density_it,n_elec} = std(RR.FPR_all);
            E.TPR_array_std{deep_it,density_it,n_elec} = std(RR.TPR_all);
            E.ACC_array_std{deep_it,density_it,n_elec} = std(RR.ACC_all);
        end
    end
end
%% boxplot with acc
clf
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [0 2 3];
DD = {'0%', '50%', '75%'};
density = [10 20];

density_it = 0;
cnt = 0;
ind_plot = [3 1 5;4 2 6];
% 221 ->set(gca, 'Position',[0.1300 0.55 0.3347 0.3370])
% 224 ->p = get(gca,'Position'); p(2) = 0.55;

for n_act = density
    density_it = density_it+1;
    deep_it = 0;

    for nd = n_deep
        deep_it = deep_it+1;
        figure(deep_it)
%         clf;
        for kk=1:4
            FPR_elec(:,kk) = FPR_array{deep_it,density_it,kk}(:);
            TPR_elec(:,kk) = TPR_array{deep_it,density_it,kk}(:);
            ACC_elec(:,kk) = ACC_array{deep_it,density_it,kk}(:);
        end
        subplot(3,2,ind_plot(density_it,1))

        boxplot(FPR_elec);
%         set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',18)

%         set(gca, 'Position',[0.1300 0.55 0.3347 0.3370])
%         xlabel('Number of electrodes')
        if density_it ==1
            ylabel('FPR')
        end
%         title(['m_active ',int2str(density(density_it)*2),'%'],'Interpreter','None')
        ylim([0 0.2])
        axis('normal')
        grid on


        subplot(3,2,ind_plot(density_it,2))
        boxplot(TPR_elec);
%         set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',18)

%         xlabel('Number of electrodes')
        if density_it ==1
            ylabel('TPR')
        end

        title([int2str(density(density_it)*2),'% active sources'],'Interpreter','None')
        ylim([0 1])
        axis('normal')
        grid on


        subplot(3,2,ind_plot(density_it,3))
        boxplot(ACC_elec);
        set(gca,'Xticklabel',{'108','61','32','19'},'FontSize',16)
        xlabel('Number of electrodes')
        if density_it ==1
            ylabel('ACC')
        end
%         title(['m_active ',int2str(density(density_it)*2),'%'],'Interpreter','None')
        ylim([0.5 1])
        axis('normal')
        grid on
        box on

%         sgtitle(['m_active ',int2str(density(density_it)*2),'%',', ',DD{deep_it},' deep sources'],'Interpreter','None')
    end

end
HH = {'0','50','75'};
for ii=1:3
        figure(ii)
 ha=get(gcf,'children');
set(findobj(ha,'type','line'),'linew',1.4)
%  set(ha(2),'position',[.5 .08+.045 .4 .4],'FontSize',18)
%  set(ha(4),'position',[.1 .08+.045 .4 .4],'FontSize',18)
%  set(ha(1),'position',[.5 .5+.055 .4 .4],'FontSize',18)
%  set(ha(3),'position',[.1 .5+.055 .4 .4],'FontSize',18)
%
 set(ha(1),'position',[.5 .1-0.05+0.05 .4 .2],'FontSize',20)
 set(ha(2),'position',[.5 .7-0.15+0.04+0.05 .4 .3],'FontSize',20)
 set(ha(3),'position',[.5 .4-0.15+0.02+0.05 .4 .3],'FontSize',20)
 set(ha(4),'position',[.1 .1-0.05+0.05 .4 .2],'FontSize',20)
 set(ha(5),'position',[.1 .7-0.15+0.04+0.05 .4 .3],'FontSize',20)
 set(ha(6),'position',[.1 .4-0.15+0.02+0.05 .4 .3],'FontSize',20)

  set(ha(3),'Xticklabel',[]);
  set(ha(6),'Xticklabel',[]);
   set(ha(1),'Yticklabel',[])
 set(ha(2),'Yticklabel',[],'Xticklabel',[])
  set(ha(5),'Xticklabel',[])
 set(ha(3),'Yticklabel',[])
%  print(figure(ii),['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_bic_deep_',HH{ii}],'-depsc')
end
% close all
% subplot(2,2,1,'Position',[ 0 0 0])

%% boxplot without acc

clf
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [0 2 3];
DD = {'0%', '50%', '75%'};
density = [10 20];

density_it = 0;
cnt = 0;
ind_plot = [3 1;4 2];
for n_act = density
    density_it = density_it+1;
    deep_it = 0;

    for nd = n_deep
        deep_it = deep_it+1;
        figure(deep_it)
%         clf;
        for kk=1:4
            FPR_elec(:,kk) = E.FPR_array{deep_it,density_it,kk}(:);
            TPR_elec(:,kk) = E.TPR_array{deep_it,density_it,kk}(:);
            ACC_elec(:,kk) = E.ACC_array{deep_it,density_it,kk}(:);
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
        ylim([0 0.3])
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
    end
end
HH = {'0','50','75'};
for ii=1:3
        figure(ii)
 ha=get(gcf,'children');
set(findobj(ha,'type','line'),'linew',1.4)
 set(ha(2),'position',[.5 .07+0.05 .4 .4],'FontSize',18)
 set(ha(4),'position',[.1 .07+0.05 .4 .4],'FontSize',18)
 set(ha(1),'position',[.5 .5+0.05 .4 .4],'FontSize',18)
 set(ha(3),'position',[.1 .5+0.05 .4 .4],'FontSize',18)
 set(ha(1),'Xticklabel',[],'Yticklabel',[])
 set(ha(3),'Xticklabel',[])
 set(ha(2),'Yticklabel',[])
%  set(gcf, 'Position', get(0, 'Screensize'));
%  print(figure(ii),['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_bic_deep_',HH{ii}],'-depsc')
%  saveas(gcf,['G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\figures\source_selection_bic_deep_',HH{ii}],'fig')
end
% close all
% subplot(2,2,1,'Position',[ 0 0 0])
