clear;
load('data/parinthorn/F_result_108');

F_true = M.F_true;
F_pn = M.bic.F_pn;
F_wn = M.bic.F_wn;
ind_F_true = cell(20,1);
for i = 1:20
    ind_F_true{i} = find(F_true{i} > 0);
end

range = 0:0.001:2;
n = length(range);

ind_F_pn = cell(20,5);
ind_F_wn = cell(20,5);
FPRpn108 = zeros(n,1);
TPRpn108 = zeros(n,1);
FPRwn108 = zeros(n,1);
TPRwn108 = zeros(n,1);
for j = 1:n
    for i = 1:20
        th = range(j)*max(max(F_true{i}));
        ind_F_pn{i,j} = find(F_pn{i,1} >= th);
        ind_F_wn{i,j} = find(F_wn{i,1} >= th);
    end
    Mpn108{j} = compare_F(ind_F_pn(:,j),ind_F_true,50);
    Mwn108{j} = compare_F(ind_F_wn(:,j),ind_F_true,50);
    FPRpn108(j) = Mpn108{j}.FPR;
    TPRpn108(j) = Mpn108{j}.TPR;
    FPRwn108(j) = Mwn108{j}.FPR;
    TPRwn108(j) = Mwn108{j}.TPR;
end
figure(1); plot(FPRpn108,TPRpn108,'b','LineWidth',1.8); axis square; hold on;
plot(FPRwn108,TPRwn108,'r','LineWidth',1.8); axis square; xlabel('FPR'); ylabel('TPR');

figure(2); 
subplot(131); imagesc(F_true{1}); colormap(1-gray); axis square; title('Ftrue');
subplot(132); imagesc(F_pn{1,1}); colormap(1-gray); axis square; title('F pinknoise');
subplot(133); imagesc(F_wn{1,1}); colormap(1-gray); axis square; title('F pinknoise');

load('data/parinthorn/F_result_61');

F_true = M.F_true;
F_pn = M.bic.F_pn;
F_wn = M.bic.F_wn;
ind_F_true = cell(20,1);
for i = 1:20
    ind_F_true{i} = find(F_true{i} > 0);
end

range = 0:0.001:2;
n = length(range);

ind_F_pn = cell(20,5);
ind_F_wn = cell(20,5);
FPRpn61 = zeros(n,1);
TPRpn61 = zeros(n,1);
FPRwn61 = zeros(n,1);
TPRwn61 = zeros(n,1);
for j = 1:n
    for i = 1:20
        th = range(j)*max(max(F_true{i}));
        ind_F_pn{i,j} = find(F_pn{i,1} >= th);
        ind_F_wn{i,j} = find(F_wn{i,1} >= th);
    end
    Mpn61{j} = compare_F(ind_F_pn(:,j),ind_F_true,50);
    Mwn61{j} = compare_F(ind_F_wn(:,j),ind_F_true,50);
    FPRpn61(j) = Mpn61{j}.FPR;
    TPRpn61(j) = Mpn61{j}.TPR;
    FPRwn61(j) = Mwn61{j}.FPR;
    TPRwn61(j) = Mwn61{j}.TPR;
end
figure(1);plot(FPRpn61,TPRpn61,'-.b','LineWidth',1.8); axis square; hold on;
plot(FPRwn61,TPRwn61,'-.r','LineWidth',1.8); axis square; xlabel('FPR'); ylabel('TPR');
legend('108 probes pinknoise','108 probes whitenoise','61 probes pinknoise','61 probes whitenoise');
for i = 1:20
F_pn_th = F_pn{i,1};
F_pn_th((setdiff(1:2500,ind_F_pn{i,2}))) = 0;
F_wn_th = F_wn{i,1};
F_wn_th((setdiff(1:2500,ind_F_wn{i,2}))) = 0;

Mpn(i) = compare_F(ind_F_pn(i,2),ind_F_true(i),50);
Mwn(i) = compare_F(ind_F_wn(i,2),ind_F_true(i),50);
TPR(i) = Mpn(i).TPR;
end
figure(3); 
subplot(131); imagesc(F_true{i}); colormap(1-gray); axis square; title('Ftrue');
subplot(132); imagesc(F_pn_th); colormap(1-gray); axis square; title('F pinknoise');
title(sprintf('Fpinknoise , TPR =%.3f, FPR =%.3f, ACC =%.3f',Mpn.TPR,Mpn.FPR,Mpn.ACC) ...
    ,'Interpreter','None')
subplot(133); imagesc(F_wn_th); colormap(1-gray); axis square; title('F whitenoise');
title(sprintf('Fwhitenoise , TPR =%.3f, FPR =%.3f, ACC =%.3f',Mwn.TPR,Mwn.FPR,Mwn.ACC) ...
    ,'Interpreter','None')






