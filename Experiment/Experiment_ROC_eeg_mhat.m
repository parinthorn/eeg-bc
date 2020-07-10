clear
clc
EEG_NOISE = {'10'};
ACTIVE_CHANNEL = [5,25,40];
SELECT_CHANNEL = [100,100,100];
% EEG_CHANNEL = [10,30,60,90];
% EEG_CHANNEL = [40,50];
EEG_CHANNEL = [40];
MODEL_SAVE_DIRECTORY = 'MODEL';
TIMESERIES_SAVE_DIRECTORY = 'DATA';
RESULT_SAVE_DIRECTORY = 'RESULT';
variance = [5];
variance_name = {'5'};
% EEG_NOISE = [10];


% ACTIVE_CHANNEL = 5;
% EEG_CHANNEL = 35;
% EEG_NOISE = 0.1;
for hh=1:1
for ii=3:length(ACTIVE_CHANNEL)
    MODEL_NAME = strcat(MODEL_SAVE_DIRECTORY,'\ss_model_',int2str(ii));
    nchan = SELECT_CHANNEL(ii);
    load(MODEL_NAME) %load model
    disp(MODEL_NAME)
    for kk=1:length(EEG_CHANNEL)
        for jj=1:length(EEG_NOISE)
            
            EEG_NAME = strcat(TIMESERIES_SAVE_DIRECTORY,'\ss_model_',int2str(ii),'_eeg_',int2str(EEG_CHANNEL(kk)),'_',EEG_NOISE{jj});
            disp(EEG_NAME)
            load(EEG_NAME) % load eeg, y.r, y.noise, y.data
            
            TP_Lpq = 0;
            FP_Lpq = 0;
            TN_Lpq = 0;
            FN_Lpq = 0;
            TP_L21 = 0;
            FP_L21 = 0;
            TN_L21 = 0;
            FN_L21 = 0;
            ROC_GRID_Lpq = zeros(75,4);
            ROC_GRID_L21 = zeros(75,4);
            for TT=1:20
                
                nstate = 5;
                Timepoints = 1000;
                
                Lhat = augmentL(y.L,nchan,variance(hh));
%                 [out_sys.ss.A, out_sys.ss.B, out_sys.ss.H, out_sys.ss.D, out_sys.ss.K, ...
%                     out_sys.est.C, out_sys.est.C_L21_bic, out_sys.est.C_Lpq, ...
%                     out_sys.est.C_L21, out_sys.est.bic_ind_Lpq, out_sys.est.bic_ind_L21] =  ...
%                     subid_eeg_Lpq(y.data(:,end-Timepoints+1:end,TT),Lhat,[],2*ceil(nstate/y.r),nstate);
%                 
                [out_sys.ss.A,out_sys.ss.B,out_sys.ss.H,out_sys.ss.D,out_sys.ss.K,out_sys.est] = ...
                    subid_eeg_Lpq(y.data(:,end-Timepoints+1:end,TT),Lhat,[],2*ceil(nstate/y.r),nstate);
                


                [TP_1,TN_1,FP_1,FN_1] = ROC_eeg(out_sys.est.C_Lpq_bic,SYSTEM(TT).sys.sys_ss.C);
                [TP_2,TN_2,FP_2,FN_2] = ROC_eeg(out_sys.est.C_L21_bic,SYSTEM(TT).sys.sys_ss.C);
                
                [val_Lpq] = batch_ROC(out_sys.est.C_Lpq,SYSTEM(TT).sys.sys_ss.C);
                [val_L21] = batch_ROC(out_sys.est.C_L21,SYSTEM(TT).sys.sys_ss.C);
                ROC_GRID_Lpq = ROC_GRID_Lpq+val_Lpq;
                ROC_GRID_L21 = ROC_GRID_L21+val_L21;
                
                TP_Lpq = TP_Lpq+TP_1;
                FP_Lpq = FP_Lpq+FP_1;
                TN_Lpq = TN_Lpq+TN_1;
                FN_Lpq = FN_Lpq+FN_1;
                TP_L21 = TP_L21+TP_2;
                FP_L21 = FP_L21+FP_2;
                TN_L21 = TN_L21+TN_2;
                FN_L21 = FN_L21+FN_2;
                estimated_model(TT) = out_sys;
            end
            FPR_L21 = FP_L21/(FP_L21+TN_L21);
            TPR_L21 = TP_L21/(TP_L21+FN_L21);
            FPR_Lpq = FP_Lpq/(FP_Lpq+TN_Lpq);
            TPR_Lpq = TP_Lpq/(TP_Lpq+FN_Lpq);
            
            TPR_CVX = (ROC_GRID_L21(:,1)./(ROC_GRID_L21(:,4)+ROC_GRID_L21(:,1)));
            FPR_CVX = (ROC_GRID_L21(:,3)./(ROC_GRID_L21(:,2)+ROC_GRID_L21(:,3)));
            TPR_nCVX = (ROC_GRID_Lpq(:,1)./(ROC_GRID_Lpq(:,4)+ROC_GRID_Lpq(:,1)));
            FPR_nCVX = (ROC_GRID_Lpq(:,3)./(ROC_GRID_Lpq(:,2)+ROC_GRID_Lpq(:,3)));
            
            sys_est.model = estimated_model;
            sys_est.T = Timepoints;
            sys_est.nstate = nstate;
            sys_est.stat.cvx_bic = [FPR_L21 TPR_L21];
            sys_est.stat.ncvx_bic = [FPR_Lpq TPR_Lpq];
            sys_est.stat.cvx_ROC = [FPR_CVX TPR_CVX];
            sys_est.stat.ncvx_ROC = [FPR_nCVX TPR_nCVX];
            sys_est.r = y.r;
            sys_est.sigma_v = y.sigma_v;
            out_file = strcat(RESULT_SAVE_DIRECTORY,'\Lvariance',variance_name{hh},'_mhat',int2str(nchan),'_ss5_estimate_model_',int2str(ii),'_eeg_',int2str(EEG_CHANNEL(kk)),'_',EEG_NOISE{jj},'T',int2str(sys_est.T));
            save(out_file,'sys_est')
            figure(2)
            plot(sys_est.stat.ncvx_ROC(:,1),sys_est.stat.ncvx_ROC(:,2),'b','LineWidth',2)
            hold
            plot(sys_est.stat.cvx_ROC(:,1),sys_est.stat.cvx_ROC(:,2),'--r','LineWidth',2)
            hold
            axis('square')
            grid on
            legend('NON-CONVEX','CONVEX')
            set(gca,'FontSize',12);
            pause(0.1)
            fprintf('FPR L21: %f, TPR L21: %f, \n FPR Lpq: %f, TPR Lpq: %f \n \n' ...
    ,FPR_L21,TPR_L21,FPR_Lpq,TPR_Lpq)


            clear sys_est
        end
    end
end
end