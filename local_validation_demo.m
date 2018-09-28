%tju cs for bioinformatics 
load('database_all.mat');
load('miRNA_name.mat');
load('number.mat');
local_case_all=[];

for ii=1:length(number)
    A=str2num(number(ii,3));
    eval(['name = ' char(number(ii,2)) ';'])
    Pre_value = local_validation(A);
    [B index]=sort(Pre_value,'descend');
    
    a = cell(50,2);

    for i= 1:50
        a{i,1} = number(ii,4);
        a{i,2} = miRNA_name{index(i)};
        [I,J]=find(upper(name(:,1:2))==upper(a{i,2}));    
        ll=unique(name(I,4));
        if length(ll)==0
            a{i,3}='unconfirmed';
        end
        if length(ll)==1
            a{i,3}=[ll{1}];
        end
        if length(ll)==2
            a{i,3}=[ll{1},'&',ll{2}];
        end
        if length(ll)==3
            a{i,3}=[ll{1},'&',ll{2},'&',ll{3}];
        end    
    end
    local_case_all=[local_case_all;a];     
end

num=[];
for ii=1:length(number)  
    tempa = local_case_all((ii-1)*50+1:ii*50,3);
    [I,J] = find(ismember(tempa, cellstr('unconfirmed')));    
    B=length(I);
    num=[num;[ii,50-B]];   
end
