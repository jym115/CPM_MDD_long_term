
network_flg=2; %%%%  whole brain(network_flg=1) or network(network_flg=2)
thresh=0.01 %%%%  edge selection threshold (0.005/0.01/0.05/0.1)
model_flg =4; %%%%  edge set (5=positive/6=negtive/4=whole edge set)
class_flag=2;%%%%  classification outcome 1=TRD 2=relapse

N_iteration=1000
n_a_all=zeros (11,N_iteration);%accuracy
n_FP_all=zeros (11,N_iteration);%false positive rate
n_TP_all=zeros (11,N_iteration);%true positive rate
n_FN_all=zeros (11,N_iteration);%false negative rate
n_TN_all=zeros (11,N_iteration);%true negative rate

%final cost parameter was generated by averaging the optimized   
% cost parameters in all the nested LOOCV to reduce overfitting
cbest=1.3; 

for n=4
    network_flg=1;
    if( network_flg ==1)
       all_conn_valid=mdd_delta_conn_110;
    elseif( network_flg==2)
       all_conn_valid=cell2mat(all_within_network(n));
    end
model_flg =4;

total=size(all_conn_valid,1);
score= group_remit;

a_all=zeros (1,N_iteration);%accuracy
FP_all=zeros (1,N_iteration);%false positive rate
TP_all=zeros (1,N_iteration);%true positive rate
TP_all=zeros(1,N_iteration);
TN_all=zeros(1,N_iteration);
 

for it = 1:N_iteration
    no_per = [randperm(110)]';
    score=group_remit(no_per,:);
    Yfit = zeros(1, total);    

 for i = 1:total
    
    train = all_conn_valid;
    train(i, :) = [];
    
    train_score = score;
    train_score(i) = [];
    
    [r, p] = corr( train, train_score,'type','spearman');
    pos_edge = find( p<thresh & r>0);
    neg_edge = find( p<thresh & r<0);
    
   
    %svm regression model
      
    	mdl = fitcsvm([train(:, pos_edge),train(:, neg_edge)],train_score,'Cost',[0 cbest;1 0]);%与假阴性相比，对假阳性应用双倍的惩罚        
        Yfit(i)=predict(mdl,[all_conn_valid(i, pos_edge),all_conn_valid(i, neg_edge)]);
        
a=find((group_remit-Yfit)==0);%accuracy       
FP=find((score-Yfit)==-2);%-1-1 false pos rate 
TP=find((group_remit==1)&(Yfit'==1));
FN=find((score-Yfit)==2);%1-(-1)  false neg rate  
TN=find((group_remit==-1)&(Yfit'==-1));
 a_all(it)=lengtha(a);
 FP_all(it)=length(FP);%false pos
 TP_all(it)=length(TP);%false neg
 FN_all(it)=length(FN);
 TN_all(it)=length(TN);
end
n_a_all(n,:)=a_all;
n_FP_all(n,:)=FP_all;
n_TP_all(n,:)=TP_all;
n_FN_all(n,:)=FN_all;
n_TN_all(n,:)=TN_all;
it
end
n
end

%% 

n_sensitivity_all=zeros (11,N_iteration);
n_specificity_all=zeros (11,N_iteration);
n_accuracy_all=zeros (11,N_iteration);
for n=1:11
 n_sensitivity_all(n,:)=(26-n_FP_all(n,:))./26*100;
 n_specificity_all(n,:)=(84-n_FN_all(n,:))./84*100;
  n_accuracy_all(n,:)=n_a_all(n,:)/110*100;
end






