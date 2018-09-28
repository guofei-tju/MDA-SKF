function d = kernel_Hamads(adjmat,dim)
%tju cs for bioinformatics 
% Calculates the Hamming distance kernel from a graph adjacency 
% matrix. If tha graph is unipartite, ka = kb.
%INPUT: 
% adjmat : binary adjacency matrix
% dim    : dimension (1 - rows, 2 - cols)
%OUTPUT:
% d : kernel matrix for adjmat over dimension 'dim'
    y = adjmat;
    % Graph based kernel
	if dim == 1
        ga = squareform(pdist(y,'hamming'));
    else
        ga = squareform(pdist(y','hamming'));
    end
    d=ones(size(ga))-ga;  
end
