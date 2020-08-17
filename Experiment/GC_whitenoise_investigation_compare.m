clear
clc
wn_folder='G:\Shared drives\MASTER_DRIVE\Journal\DATA\WHITE_NOISE_DATASET\';
pn_folder='G:\Shared drives\MASTER_DRIVE\Journal\DATA\PINK_NOISE_DATASET\';
ii=1;

cnt = 0;
source_location = ["_shallow","_deep"];
for sloc= source_location
    cnt = cnt+1;
    for ii=1:10
        load([pn_folder,'model_',int2str(ii),sloc{1}])
        load([pn_folder,'eegdata_',int2str(ii),sloc{1}])
        data_pn = eegdata.EEG_data;
        load([wn_folder,'eegdata_',int2str(ii),sloc{1}])
        data_wn = eegdata.EEG_data;
        
        L = model.L0;
        r = size(L,1);
        ind_chan = 1:1:r;
        
        nstate = (model.PARAMETER.lag+ ...
            2*model.PARAMETER.nbutter)*model.PARAMETER.m_active;
        ind_Ctrue = find(norms(model.source_model0.C,2,2));
        Timepoints = 15000;
        tic;
        for tt=1:5 % trial iteration
            data_in = data_pn(ind_chan,end-Timepoints+1:end,tt);
            [sys_pn(tt).sys,sys_pn(tt).C_out] = ...
                subid_eeg_Lpq(data_in ...
                ,L,ind_Ctrue,2*ceil(nstate/r),nstate);
            
            data_in = data_wn(ind_chan,end-Timepoints+1:end,tt);
            [sys_wn(tt).sys,sys_wn(tt).C_out] = ...
                subid_eeg_Lpq(data_in ...
                ,L,ind_Ctrue,2*ceil(nstate/r),nstate);
            ind = sys_pn(tt).C_out.ind_chosen_Lpq.cheated;
            A = sys_pn(tt).sys.A;
            % C = sys_pn.C_out.C_Lpq_CLS(:,:,ind);
            C = sys_pn(tt).C_out.C_cheated;
            E = sys_pn(tt).C_out.E_Lpq_CLS(:,:,ind);
            W = sys_pn(tt).C_out.W_Lpq_CLS(:,:,ind);
            [N,V,obj] = noisecovest(E,L,'homo');
            [F_pn,~,~] = calgcss(A,C,W,N,[]);
            ind = sys_wn(tt).C_out.ind_chosen_Lpq.cheated;
            A = sys_wn(tt).sys.A;
            % C = sys_wn.C_out.C_Lpq_CLS(:,:,ind);
            C = sys_wn(tt).C_out.C_cheated;
            E = sys_wn(tt).C_out.E_Lpq_CLS(:,:,ind);
            W = sys_wn(tt).C_out.W_Lpq_CLS(:,:,ind);
            [N,V,obj] = noisecovest(E,L,'homo');
            [F_wn,~,~] = calgcss(A,C,W,N,[]);
            M(cnt).F_pn{ii,tt} = find(F_pn);
            M(cnt).F_wn{ii,tt} = find(F_wn);
        end
        toc;
        save([wn_folder,'result',int2str(ii),sloc{1}],'sys_wn')
        save([pn_folder,'result',int2str(ii),sloc{1}],'sys_pn')
        clear sys_wn sys_pn
    end
    
end
%%



%%
% clc
% tmp_true{1} = ind_Ctrue;
% tmp_pn{1} = sys_pn.C_out.nz_ind_C_Lpq;
% tmp_wn{1} = sys_wn.C_out.nz_ind_C_Lpq;
% acc_pn = compare_C(tmp_pn,tmp_true,size(L,2),75);
% acc_wn = compare_C(tmp_wn,tmp_true,size(L,2),75);
% figure(1)
% subplot(1,2,1)
% plot(acc_pn.FPR,acc_pn.TPR)
% subplot(1,2,2)
% plot(acc_wn.FPR,acc_wn.TPR)
% pause(0.1)
%%

%%

%%
% figure(3)
% subplot(1,3,1)
% imagesc((F_pn))
% axis('square')
% 
% subplot(1,3,2)
% imagesc((F_wn))
% axis('square')
% subplot(1,3,3)
% imagesc(model.F0)
% axis('square')