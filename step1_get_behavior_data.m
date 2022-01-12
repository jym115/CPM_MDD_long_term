%% 2021-08-14
load /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final/cpm_TRD.mat
DataDir = '/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/110mddvs143hc';
SaveFile = '/Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/matlab_paper_final/cpm_TRD.mat';
%1.get behavior data
load /Users/yumeng/Desktop/1/1zhumadian/network/1relapse_matlab/dmn_nbs/110mddvs143hc/DMN_remission.mat 'remit_235' 'dmn20210717' 
remit_235=dmn20210717(177:(176+235),[3,5:7,23,24]);
    remit=remit_235(find (remit_235(:,1)>-2),:);
    id110=1:235;
    id110=id110(find (remit_235(:,1)>-2));
    id110=2000+id110;
    group=remit(:,1);
    group_remit=group;
    design=[group_remit,group_remit,remit(:,2:end)];
    demomdd=design(:,[3:5,7]);
save(SaveFile, 'demomdd', 'group_remit', 'dmn20210717')


