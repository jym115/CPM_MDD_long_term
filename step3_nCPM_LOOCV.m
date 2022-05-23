%% 3.nCPM LOOCV of whole brain and 11 networks
    
network_flg=2; %%%%  whole brain(network_flg=1) or network(network_flg=2)
thresh=0.1 %%%%  edge selection threshold (0.005/0.01/0.05/0.1)
class_flag=2;%%%%  classification outcome 1=TRD 2=relapse

    pos=ones(1,length(all_conn_valid));
    neg=ones(1,length(all_conn_valid));
    pos_all=zeros(1,length(all_conn_valid));
    neg_all=zeros(1,length(all_conn_valid));

 if class_flag==1
    score=group_remit;
 elseif class_flag==2
    score=group_relapse;
     elseif class_flag==3
    score=aremit;
 end
 
 ynan=zeros(length(n),1);  
 sensitivity=zeros(length(n),1);  
 specificity=zeros(length(n),1);  
 FP=zeros(length(n),1);  
 FN=zeros(length(n),1);  
 accuracy=zeros(length(n),1);  
 FP_all=zeros(length(n),1);  
 FN_all=zeros(length(n),1);  

for n= 4    % n=1:11
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
   %select edge correlated to treatment outcome as features
  [r, p] = corr( train, train_score, 'type','spearman');%    );%kendall
% [r, p] = partialcorr( train, train_score,train_ffd,'type','spearman');%);

    pos_edge = find( p<thresh & r>0);
    neg_edge = find( p<thresh & r<0); 
    %fit svm model
        if isempty([train(:, pos_edge),train(:, neg_edge)])
        Yfit(i,:)=nan;
        else
        cbest=findbestc(all_edges, all_behav, thresh, k,cmin,cstep,cmax);
    	mdl = fitcsvm([train(:, pos_edge),train(:, neg_edge)],train_score,'Cost',[0 cbest;1 0]);          
        Yfit(i)=predict(mdl,[all_conn_valid(i, pos_edge),all_conn_valid(i, neg_edge)]);
        end
    
    pos_mask = zeros(1,length(all_conn_valid)); %(table3-6 svm weight and frequencies)
    pos_mask(pos_edge) = 1;
    neg_mask = zeros(1,length(all_conn_valid)); 
    neg_mask(neg_edge) = 1;

    pos=pos_mask.*pos;
    pos_all=pos_all+pos_mask;
    neg=neg_mask.*neg;
    neg_all=neg_all+neg_mask;

 end

 Yfit(find(isnan(Yfit)==1))=4;
 a=find((score-Yfit)==0);%accuracy
 FP=find((score-Yfit)==-2);%-1-1 false pos
 FN=find((score-Yfit)==2);%1-(-1)  false neg

 a_all(n,:)=length(a);
 FP_all(n,:)=length(FP);
 FN_all(n,:)=length(FN);
 if class_flag==1
 sensitivity(n,:)=(26-FP_all(n,:))./26*100
 specificity(n,:)=(84-FN_all(n,:))./84*100
 accuracy(n,:)=a_all(n,:)/110*100

 elseif class_flag==2
 specificity(n,:)=(40-FP_all(n,:))./40*100;
 sensitivity(n,:)=(28-FN_all(n,:))./28*100;
 accuracy(n,:)=a_all(n,:)/68*100;

 end


end


%% table LOOCV positive edge svm frequencies

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
 %% 3.creat positive edge visualization
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
 %% 3.creat negative edge visualization
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
 

 
