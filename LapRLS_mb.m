function [LapA] = LapRLS_mb(W1,W2,y, lambda,p_nearest_neighbor,type)
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
%[1] Xia Z, Wu L Y, Zhou X, et al. 
%      Semi-supervised drug-protein interaction prediction from heterogeneous biological spaces[J]. 
%           Bmc Systems Biology, 2010, 4(S2):1-16.
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% y  : binary adjacency matrix, (m-by-n)
%lambda: Regularized item (4)
%p_nearest_neighbor: the p nearest neighbor samples (26)

%Network of Laplacian Regularized Least Square
[num_1,num_2] = size(y);
%1.Sparsification of the similarity matrices
%1.1
fprintf('Caculating nearest neighbor graph 1\n');
%N_1 = neighborhood_Com(W1,p_nearest_neighbor);
%fprintf('Sparsification of the similarity matrix 1\n');
%S_1 = N_1.*W1;
S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;

d_tmep_1=eye(num_1)/(D_1^(1/2));
L_D_11 = d_tmep_1*L_D_1*d_tmep_1;

%1.2
%fprintf('Caculating nearest neighbor graph 2\n');
%N_2 = neighborhood_Com(W2,p_nearest_neighbor);
%fprintf('Sparsification of the similarity matrix 2\n');
%S_2 = N_2.*W2;
S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;

d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;

fprintf('Laplacian Regularized Least Square\n');
%A_1 = (W1/(W1 + lambda*L_D_1*W1))*y;
%A_2 = (W2/(W2 + lambda*L_D_2*W2))*y';

A_1 = W1*pinv(W1 + lambda*L_D_11*W1)*y;
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*y';

LapA =[];
if type==1
	LapA = (A_1 + A_2')/2;
else
	A_2 = A_2';
	LapA=A_1.*(A_1>=A_2)+A_2.*(A_1<A_2);
end
end


function similarities_N = neighborhood_Com(similar_m,kk)

similarities_N=zeros(size(similar_m));

mm = size(similar_m,1);

for ii=1:mm
	
	for jj=ii:mm
		iu = similar_m(ii,:);
		iu_list = sort(iu,'descend');
		iu_nearest_list_end = iu_list(kk);
		
		ju = similar_m(:,jj);
		ju_list = sort(ju,'descend');
		ju_nearest_list_end = ju_list(kk);
		if similar_m(ii,jj)>=iu_nearest_list_end & similar_m(ii,jj)>=ju_nearest_list_end
			similarities_N(ii,jj) = 1;
			similarities_N(jj,ii) = 1;
		elseif similar_m(ii,jj)<iu_nearest_list_end & similar_m(ii,jj)<ju_nearest_list_end
			similarities_N(ii,jj) = 0;
			similarities_N(jj,ii) = 0;
		else
			similarities_N(ii,jj) = 0.5;
			similarities_N(jj,ii) = 0.5;
		end
	
	
	end


end

end
