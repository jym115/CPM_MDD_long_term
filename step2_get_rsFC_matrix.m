%% 2.get rsFC matrix of whole brain and 11 networks
  mdd_len=110;
  no_node = 256;
aa = ones(no_node, no_node);
aa_upp = triu(aa, 1);
upp_id = find( aa_upp>0);
upp_len = length(upp_id);

mdd_baseline_conn_110=zeros(mdd_len,upp_len);
mdd_followup_conn_110=zeros(mdd_len,upp_len);
mdd_delta_conn_110=zeros(mdd_len,upp_len);
count = 0;
cur_d = dir('/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/dmn110/mdd_2time_matrices');
cd /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/dmn110/mdd_2time_matrices
for i = 1:mdd_len
    cur=cur_d(i+3,:);
            cur_mat = textread([cur.name]);
            count = count+1;
            mdd_baseline_conn_110(count,:) = cur_mat(upp_id);
           
            cur=cur_d(i+3+110,:);
            cur_mat = textread([cur.name]);
            mdd_followup_conn_110(count,:) = cur_mat(upp_id);
            mdd_delta_conn_110(count,:) =mdd_baseline_conn_110(count,:)-mdd_followup_conn_110(count,:);
            

end
%% get delta rsFC matrix of whole brain and 11 networks
  network_flg=3;  
  if( network_flg ==1)
       all_conn_valid=mdd_conn_all;
    elseif( network_flg==2)
       all_conn_valid=mdd_followup_conn_110;
     elseif(network_flg==3) 
         all_conn_valid=mdd_delta_conn_110;
  end
no_networks=11;
load /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/cpm.mat 'ten_net' 'thirteen_net' 'thirteen_network'
ten_net=thirteen_net;
no_sub=length(all_conn_valid(:,1))

% within network

for mm = 1:no_networks
    for k=1:no_sub
        matrix=all_conn_valid(k,:);
        matr = zeros(no_node, no_node);
        matr (upp_id)=matrix;
        matr=matr+matr';
   
        new_assignments_matr = matr(ten_net(:,2),:);%regroup nodes into network 
        new_assignments_final_matr = new_assignments_matr(:,ten_net(:,2));%regroup nodes into network 
        
        [indices] = find( ten_net(:, 3)==mm);
        
            network = new_assignments_final_matr(indices,indices);
            upper_network = triu(network, 1);
            upp_id1 = find(upper_network~=0);
            matr_upp = network(upp_id1);
            all_se_network(:,k) = matr_upp;
        end
      all_within_network{mm}=all_se_network';
      clear all_se_network
end



%% 

no_sub=length(all_conn_valid(:,1))

for mm = 1:no_networks
    for i=1:mm 
     for k=1:no_sub
%          matrix=pos_all;
        matrix=all_conn_valid(k,:);
        matr = zeros(no_node, no_node);
        matr (upp_id)=matrix;
        matr=matr+matr';
   
        new_assignments_matr = matr(ten_net(:,2),:);%regroup nodes into network 
        new_assignments_final_matr = new_assignments_matr(:,ten_net(:,2));%regroup nodes into network 
 
            [indices] = find( ten_net(:, 3)==i);
            [indices2] = find( ten_net(:, 3)==mm);
            network = new_assignments_final_matr(indices,indices2);
            
            if (i==4) && (mm ==4)
             network = new_assignments_final_matr(indices,indices2);
            end
            
            aa = ones(length(indices), length(indices2));
            upp_id1 = find( aa>0);
             matr_upp = network(upp_id1);
             
              all_se_network(k,:) = matr_upp;           
        end
         all_network{i,mm} = all_se_network;
         clear all_se_network
    end   
    
end

%% 

cd /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final
SaveFile = '/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final/cpm_relapse.mat';

save(SaveFile)


