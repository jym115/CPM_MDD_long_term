%% 3.nCPM LOOCV of whole brain and 11 networks
    
network_flg=2; %%%%  whole brain(network_flg=1) or network(network_flg=2)
thresh=0.1 %%%%  edge selection threshold (0.005/0.01/0.05/0.1)
model_flg =4; %%%%  edge set (5=positive/6=negtive/4=whole edge set)
class_flag=2;%%%%  classification outcome 1=TRD 2=relapse

    pos=ones(1,length(all_conn_valid));
    neg=ones(1,length(all_conn_valid));
    pos_all=zeros(1,length(all_conn_valid));
    weight=zeros(1,length(all_conn_valid));
    weight_all=zeros(1,length(all_conn_valid));
    neg_all=zeros(1,length(all_conn_valid));
 if class_flag==1
    score=group_remit;
 elseif class_flag==2
    score=group_relapse;
     elseif class_flag==3
    score=aremit;

 end

%  ffd(find (isnan(ffd)==1) )=0.1;
%  med=treatment( (id110-2000) ,1);
 
 ynan=zeros(length(n),1);  
 sensitivity=zeros(length(n),1);  
 specificity=zeros(length(n),1);  
 false_positive=zeros(length(n),1);  
 false_negative=zeros(length(n),1);  
 accuracy=zeros(length(n),1);  
 a_all=zeros(length(n),1);  
 aa_all=zeros(length(n),1);  
 aaa_all=zeros(length(n),1);  
 aaaa_all=zeros(length(n),1);  
 aaaaa_all=zeros(length(n),1);  

for n= 4
    if( network_flg ==1)
        if class_flag==1
        all_conn_valid=mdd_delta_conn_110;
        elseif class_flag==2
         all_conn_valid=mdd_delta_conn_relapse;
        end
    elseif( network_flg==2)
    all_conn_valid=cell2mat(all_within_network(n));
    end
    total=size(all_conn_valid,1);
    Yfit = zeros(total,1);  
 
 for i = 1:total
    train = all_conn_valid;
    train(i, :) = [];
    train_score = score;
    train_score(i) = [];
%     train_ffd = med; %ffd med
%     train_ffd(i,:)=[];
    
  [r, p] = corr( train, train_score, 'type','spearman');%    );%kendall
%     [r, p] = partialcorr( train, train_score,train_ffd,'type','spearman');%);%, %    );%

    pos_edge = find( p<thresh & r>0);
    neg_edge = find( p<thresh & r<0);
    
    %svm model
    if( model_flg ==4) 
        if isempty([train(:, pos_edge),train(:, neg_edge)])
        Yfit(i,:)=nan;
        else
    	mdl = fitcsvm([train(:, pos_edge),train(:, neg_edge)],train_score);          
        Yfit(i)=predict(mdl,[all_conn_valid(i, pos_edge),all_conn_valid(i, neg_edge)]);
        end
    elseif( model_flg ==5) %只取正相关
        if isempty([train(:, pos_edge)])
        Yfit(i,:)=nan;
        else
    	mdl = fitcsvm([train(:, pos_edge)],train_score);          
        Yfit(i)=predict(mdl,[all_conn_valid(i, pos_edge)]);
        end
    elseif( model_flg ==6) %%只取负相相关
        if isempty([train(:, neg_edge)])
        Yfit(i,:)=nan;
        else
    	mdl = fitcsvm([train(:, neg_edge)],train_score);          
        Yfit(i)=predict(mdl,[all_conn_valid(i, neg_edge)]);
        end    
    end
    
    pos_mask = zeros(1,length(all_conn_valid)); %(table3-6 svm weight and frequencies)
    pos_mask(pos_edge) = 1;
    neg_mask = zeros(1,length(all_conn_valid)); 
    neg_mask(neg_edge) = 1;

%     weight(pos_edge) = mdl.Beta;
%     weight_all=weight+weight_all;
    pos=pos_mask.*pos;
    pos_all=pos_all+pos_mask;
    neg=neg_mask.*neg;
    neg_all=neg_all+neg_mask;

 end
 
 
 Yfit(find(isnan(Yfit)==1))=4;
 a=find((score-Yfit)==0);%accuracy
 aa=find((score-Yfit)==-2);%-1-1 false pos
 aaa=find((score-Yfit)==2);%1-(-1)  false neg
 aaaa=find((score-Yfit)==-5)%-1-4 false pos
 aaaaa=find((score-Yfit)==-3)%1-4 false neg

 a_all(n,:)=length(a);
 aa_all(n,:)=length(aa)+length(aaaa);
 aaa_all(n,:)=length(aaa)+length(aaaaa);
 aaaa_all(n,:)=length(aaaa);
 aaaaa_all(n,:)=length(aaaaa);
 
 if class_flag==1
 sensitivity(n,:)=(26-aa_all(n,:))./26*100
 specificity(n,:)=(84-aaa_all(n,:))./84*100
  accuracy(n,:)=a_all(n,:)/110*100

 elseif class_flag==2
specificity(n,:)=(40-aa_all(n,:))./40*100;
 sensitivity(n,:)=(28-aaa_all(n,:))./28*100;
 accuracy(n,:)=a_all(n,:)/68*100;

 end
 false_positive(n,:)=100-specificity(n,:);
 false_negative(n,:)=100-sensitivity(n,:);

end


%  ynan=aaaa_all+aaaaa_all;
%   [accuracy,sensitivity,specificity,ynan]

%% 
specificity_spearman_005=specificity
sensitivity_spearman_005=sensitivity
accuracy_spearman_005=accuracy

%% 
specificity_pearson_005=specificity
sensitivity_pearson_005=sensitivity
accuracy_pearson_005=accuracy

%% table LOOCV positive edge svm frequencies
%  weight_mean=weight_all./pos_all;
%  a=weight_mean(find (weight_mean>-100))';

net=cell2mat(all_within_network(4));
 b=net(:,find(pos_all>0));
 c=mean( b(find(group_remit==1),:) ) 
 d= std( b(find(group_remit==1),:) )
 
 e=mean( b(find(group_remit==-1),:) )
 f=std( b(find(group_remit==-1),:) )
 g=pos_all(find(pos_all>0)');
%  e=mean(mdd_followup_conn_110(:,find(pos_all>0)) )';
table4=[g',c',d',e',f'];


%% table LOOCV negative edge svm frequencies

 b=all_conn_valid(:,find(neg_all>0));
 c=mean( b(find(group_remit==1),:) )'
 d= std( b(find(group_remit==1),:) )'
 
 e=mean( b(find(group_remit==-1),:) )'
 f=std( b(find(group_remit==-1),:) )'
 g=neg_all(find(neg_all>0))';
%  e=mean(mdd_followup_conn_110(:,find(neg_all>0)) )';
table4=[g,c,d,e,f];
 %% plot back   %% 3.creat positive edge visualization
                  pos_mask=find(pos>0) ;
                  % pos_all are edge appear in all LOOCV
                  % pos are edges appear in each LOOCV
                  no_node=50;
                  aa=ones(no_node,no_node);
                  aaa=triu(aa,1);
                  upp_id=find (aaa>0);
                  pos_stable = zeros(no_node, no_node);
                  pos_stable(upp_id(pos_mask)) = 1;
                  pos_stable = pos_stable+pos_stable';
                  % %
                  node_id=find(thirteen_network(:,2)==5);

                  insert_id=setdiff(1:268,thirteen_network(node_id,1));
                  rearrange=pos_stable;
                  for i=1:length(insert_id)
                  h=zeros(1,length(rearrange));
                  rearrange=[rearrange(1:(insert_id(i)-1),:) ; h ; rearrange(insert_id(i):end,:)];
                  l=zeros(length(rearrange(:,1)),1);
                  rearrange=[rearrange(:,1:(insert_id(i)-1)), l, rearrange(:,insert_id(i):end)];
                  end
                  length(find(rearrange>0))
                  save (['pos_1relapseconsensus_edge.txt'], 'rearrange','-ascii');
 % visualize in bioimagesuite 
  %% plot back   %% 3.creat negative edge visualization
                  neg_mask=find(neg>0) ;
                  % neg_all are edge appear in all LOOCV
                  % neg are edges appear in each LOOCV
                  no_node=50;
                  aa=ones(no_node,no_node);
                  aaa=triu(aa,1);
                  upp_id=find (aaa>0);
                  neg_stable = zeros(no_node, no_node);
                  neg_stable(upp_id(neg_mask)) = 1;
                  neg_stable = neg_stable+neg_stable';
                  % %
                  node_id=find(thirteen_network(:,2)==5);

                  insert_id=setdiff(1:268,thirteen_network(node_id,1));
                  rearrange=neg_stable;
                  for i=1:length(insert_id)
                  h=zeros(1,length(rearrange));
                  rearrange=[rearrange(1:(insert_id(i)-1),:) ; h ; rearrange(insert_id(i):end,:)];
                  l=zeros(length(rearrange(:,1)),1);
                  rearrange=[rearrange(:,1:(insert_id(i)-1)), l, rearrange(:,insert_id(i):end)];
                  end
                  length(find(rearrange>0))
                  save (['neg_1relapseconsensus_edge.txt'], 'rearrange','-ascii');
 % visualize in bioimagesuite  
 
 %% table for svm edges
mask = textread(['pos_01_edge.txt']);%268*268
load  /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/cpm_thre005.mat 'namelistS3'
figure; imagesc(mask)

id=find (mask>0);
length(id)
no_node=268;
aa=ones(no_node,no_node);
 aaa=triu(aa,1);
upp_id=find (aaa>0);

id=intersect(id,upp_id);
mask=zeros(268,268);
posnode_id=zeros(length(id),2);
posnetwork_1=zeros(length(id),5);
posnetwork_2=zeros(length(id),5);

for i = 1:length(id)
 mask(id(i))=1;   
 
 for ii=1:268   
   for a=1:268
  if isempty(find (mask(a,ii)>0))
    b=0;
 else
    
 posnode_id(i,1)=a;
  posnode_id(i,2)=ii;

  id_node=(find (namelistS3(:,1)==a));
  posnetwork_1(i,:)=namelistS3(a,:);
  id_node=(find (namelistS3(:,1)==ii));
  posnetwork_2(i,:)=namelistS3(ii,:);

 end
 end
     
  end
 end

setdiff(posnetwork_1(:,1),thirteen_network(node_id,1)) %  check if all nodes were in dmn network
%%%%copy posnetwork_1 posnetwork_2 network names to excels and then to tables

                  %% table for svm edges
mask = textread(['neg_01_edge.txt']);%268*268
load  /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/cpm_thre005.mat 'namelistS3'
figure; imagesc(mask)

id=find (mask>0);
length(id)
no_node=268;
aa=ones(no_node,no_node);
 aaa=triu(aa,1);
upp_id=find (aaa>0);

id=intersect(id,upp_id);
mask=zeros(268,268);
negnode_id=zeros(length(id),2);
negnetwork_1=zeros(length(id),5);
negnetwork_2=zeros(length(id),5);

for i = 1:length(id)
 mask(id(i))=1;   
 
 for ii=1:268   
   for a=1:268
  if isempty(find (mask(a,ii)>0))
    b=0;
 else
    
 negnode_id(i,1)=a;
  negnode_id(i,2)=ii;

  id_node=(find (namelistS3(:,1)==a));
  negnetwork_1(i,:)=namelistS3(a,:);
  id_node=(find (namelistS3(:,1)==ii));
  negnetwork_2(i,:)=namelistS3(ii,:);

 end
 end
     
  end
 end

setdiff(negnetwork_1(:,1),thirteen_network(node_id,1)) %  check if all nodes were in dmn network
%%%%copy negnetwork_1 negnetwork_2 network names to excels and then to tables
 





 