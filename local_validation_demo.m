% colon cancer  No.91
%gastric cancer	No.140	
%pancreatic cancer	No.307	
%colorectal cancer	No.93	
%esophageal cancer No.126	
%kidney cancer	No.205
%breast cancer No.50
%lymphoma No.240


%%%The local validation 
load('Three_database_8.mat');
load('miRNA_name.mat');
A=240;
name = lymphoma;
Pre_value = local_validation(A);
[B index]=sort(Pre_value,'descend');


a = cell(50,2);

for i= 1:50
    a{i,1} = miRNA_name{index(i)};
    [I,J]=find(upper(name(:,1:2))==upper(a{i,1}));    
    ll=name(I,3);
    ll=unique(ll);
    if length(ll)==1
        a{i,2}=[ll{1}];
    end
    if length(ll)==2
        a{i,2}=[ll{1},'&',ll{2}];
    end
    if length(ll)==3
        a{i,2}=[ll{1},'&',ll{2},'&',ll{3}];
    end    
end


function [Pre_value]=local_validation(A)
load('You_dataset.mat');
y_train = miRNA_disease_Y;
[II,JJ] = find(miRNA_disease_Y == 1);
K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;
Pre_value = [];

y_train = miRNA_disease_Y;
y_train(:,A) = 0;    
K1(:,:,3)=kernel_gip(y_train,1, 1);
K2(:,:,3)=kernel_gip(y_train,2, 1);
K_COM1=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},241,10,0.1);		
K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},241,10,0.1);		
[F_1] = LapRLS_mb(K_COM1,K_COM2,y_train, 2^(-2),241,1);
Pre_value = F_1(:,A);
end



