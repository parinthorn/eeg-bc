clc
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('G:\Shared drives\MASTER_DRIVE\Journal\DATA_GENERATION\BBCB_code\data\sa')
end
foldername = 'PINK_NOISE_DATASET';
for ii=1:10
    fprintf('Generating Model:  %d \n',ii)
    generate_eegdata(5,foldername,0,[int2str(ii),'_shallow'],sa,1);
end

%%
for ii=1:10
    fprintf('loading Model:  %d \n',ii)
    load(['G:\Shared drives\MASTER_DRIVE\Journal\DATA\',foldername,'\model_',int2str(ii),'_shallow'])
    generate_eegdata_from_model(5,model,'WHITE_NOISE_DATASET',0,[int2str(ii),'_shallow'],sa,0);
%     generate_eegdata_from_model(nbatch,foldername,model,plotting,postfix)
end