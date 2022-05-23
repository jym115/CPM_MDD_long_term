%% 2021-08-14
load /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final/cpm_TRD.mat
DataDir = '/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/110mddvs143hc';
SaveFile = '/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final/cpm_TRD.mat';
%1.get behavior data
load /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/110mddvs143hc/DMN_remission.mat 'remit_235' 'dmn20210717' 
remit_235=dmn20210717(177:(176+235),[3:24]);

    remit=remit_235(find (remit_235(:,1)>-2),:);
    group_remit=remit(:,1);%remitter = 1 non-remitter = -1
    demomdd=remit(:,[3:5,7]); %age gender education FFD,ect...
save(SaveFile, 'demomdd', 'group_remit', 'dmn20210717')


