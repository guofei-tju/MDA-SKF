clear;
load('You_dataset.mat');
MU = 0.005;GAMMA=0.002;
y_train = miRNA_disease_Y;
[II,JJ] = find(miRNA_disease_Y == 1);
K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;
[n_miRNA,n_Diseas]=size(y_train);

Pre_value = [];
for i = 1:n_Diseas  
    y_train = miRNA_disease_Y;
    y_train(:,i) = 0;    
     K1(:,:,3)=kernel_gip(y_train,1, 1);
     K2(:,:,3)=kernel_gip(y_train,2, 1);

    K_COM1=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},241,10,0.1);		
    K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},241,10,0.1);		
    [F_1] = LapRLS_mb(K_COM1,K_COM2,y_train, 2^(-2),241,1);
    Pre_value = [Pre_value,F_1(:,i)];
     if mod(i,50) == 0
         i             
     end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1);

