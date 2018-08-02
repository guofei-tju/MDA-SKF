clear;
load('You_dataset.mat');
y = miRNA_disease_Y;

K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;

nfolds =5; nruns=1;
crossval_idx = crossvalind('Kfold',y(:),nfolds);
for fold = 1:nfolds
    y_train = miRNA_disease_Y;
    test_idx  = find(crossval_idx==fold);
    y_train(test_idx) = 0;
    K1(:,:,3)=kernel_gip(y_train,1, 1);
    K2(:,:,3)=kernel_gip(y_train,2, 1);
    K_COM1=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},50,10,0.1);		
    K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},50,10,0.1);		
    [F_1] = LapRLS_mb(K_COM1,K_COM2,y_train, 2^(-5),50,1);
    y(test_idx)= F_1(test_idx);
end
        
[X_1,Y_1,tpr,aupr_F_1] = perfcurve(miRNA_disease_Y(:),y(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_F_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(miRNA_disease_Y(:),y(:),1);
