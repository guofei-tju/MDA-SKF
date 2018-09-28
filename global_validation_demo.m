%tju cs for bioinformatics 
load('database_all.mat');
load('miRNA_name.mat');
load('number.mat');

Pre_value = global_validation;
gloabl_case_all=[];
for ii=1:length(number)
    A=str2num(number(ii,3));
    eval(['name = ' char(number(ii,2)) ';'])
    [I]=find(name(:,end)== 'dbDEMC');
    [J]=find(name(:,end)== 'mir2disease');
    name=name([I;J],:);
    [B index]=sort(Pre_value(:,A),'descend');
    a = cell(50,3);
    for i= 1:50
        a{i,1} = number(ii,4);
        a{i,2} = miRNA_name{index(i)};
        [I,J]=find(upper(name(:,1:2))==upper(a{i,2}));
        ll=unique(name(I,4));
        if length(ll)==1
            a{i,3}=[ll{1}];
        end
        if length(ll)==2
            a{i,3}=[ll{1},'&',ll{2}];
        end
        if length(ll)==0
            a{i,3}='unconfirmed';
        end
    end
    gloabl_case_all=[gloabl_case_all;a];
end
num=[];
for ii=1:length(number)
tempa = gloabl_case_all((ii-1)*50+1:ii*50,3);
%     [I,J] = find(tempa ==  cellstr('unconfirmed'));
[I,J] = find(ismember(tempa, cellstr('unconfirmed')));
B=length(I);
num=[num;[ii,50-B]];
end