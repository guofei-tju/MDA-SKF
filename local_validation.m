function [Pre_value]=local_validation(A) 
%tju cs for bioinformatics 
    load('You_dataset.mat');
    y_train = miRNA_disease_Y;
    [II,JJ] = find(miRNA_disease_Y == 1);
    K1 = [];
    K1(:,:,1)=miRNA_Function_S;
    K1(:,:,2)=miRNA_Sequences_Needle_S;
    K2 = [];
    K2(:,:,1)=disease_Function_S;
    K2(:,:,2)=disease_Sem_S;

    y_train = miRNA_disease_Y;
    y_train(:,A) = 0;    
    K1(:,:,3)=kernel_Hamads(y_train,1 ); %Hamming distance
    K2(:,:,3)=kernel_Hamads(y_train,2 );
    K_COM1=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},192,10,0.1);		
    K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},192,10,0.1);		
    [F_1] = LapRLS_mb(K_COM1,K_COM2,y_train, 2^(-1),1);
    Pre_value = F_1(:,A);
end