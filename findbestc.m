

function [cbest]=findbestc(all_edges, all_behav, thresh, k,cmin,cstep,cmax)

c_all=[cmin:cstep:cmax];
cfold=length(c_all);
acc_c=zeros(1,cfold);
total=length(all_behav);

for m = 1:cfold
    ctest=c_all(m) ; %penalized index      
    Yfit_c=zeros(1,total); 
    indices = crossvalind('Kfold', total-1, k);


          for ii = 1:k
             
            test_idx = (indices==ii);
            train_idx = (indices~=ii);
            train_c = all_edges(train_idx,:);
            train_score_c = all_behav(train_idx, :);
            test_c =all_edges(test_idx,:);
            
                [r, p] = corr( train_c, train_score_c);%feature selection
                pos_edge_c = find( p<thresh & r>0); %functional connectivity
                neg_edge_c = find( p<thresh & r<0);
                
                if isempty([train_c(:, pos_edge_c),train_c(:, neg_edge_c)])
                    Yfit_c(test_idx)=nan;
                else
                
                mdl_c = fitcsvm([train_c(:, pos_edge_c),train_c(:, neg_edge_c)],train_score_c,'Cost',[0 ctest;1 0]);%与假阴性相比，对假阳性应用ctest倍的惩罚        
                Yfit_c(test_idx)=predict(mdl_c,[test_c(:, pos_edge_c),test_c(:, neg_edge_c)]); 
                end
           end

                acc=find(all_behav==Yfit_c');%accuracy
                acc_c(m)=length(acc) ;  %store accuracy for different c  
end

        for mm= 1:cfold
            if mean(acc_c)==max(acc_c)
             cbest=1;
            elseif length( find( acc_c==max(acc_c)) )>1%if cbest reaches max  
                cbestid=find( acc_c==max(acc_c));
                cbest=c_all(cbestid(1));
            elseif acc_c(mm)==max(acc_c)
             cbest=c_all(mm);
            end
        end

