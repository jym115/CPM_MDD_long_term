function create_roc(xn, yn)


%%we have run all the LOOCV and copy into a excel
%%load the excel data 
thr_label=4
if thr_label==1
    yn=results(1:12,4) %sensitivty . +7 +14
    xn=100-results(1:12,3) %false positive 

elseif thr_label==2
    yn=results(17:17+11,4)
    xn=100-results(17:17+11,3)
        

elseif thr_label==3
    yn=results(33:33+11,4)
    xn=100-results(33:33+11,3)
       

elseif thr_label==4
    yn=results(49:49+11,4)
    xn=100-results(49:49+11,3)
   
end
label=['SMN';'CON';'AUN';'DMN';'VN ';'FPN';'SN ';'SCN';'VAN';'DAN';'CBN';'WB '];

%% 

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot (0,0);     %先把图标框生成
hold on;    
for i=1:length(xn)     %我是逐点画图的，因为要生成不同的颜色
    if yn(i)>0 
        plot(xn(i),yn(i),'sk','MarkerSize',60,'Marker','.','LineStyle','none','color','black');%,,'color','red'
    else
        a=1; plot(xn(i),yn(i),'sg','MarkerSize',60,'Marker','.','LineStyle','none','color','black');%,'color
    end
%     text(xn(i)+1.5,yn(i)+0.5,label(i,:),'fontsize',20);%给每个点加上标注,为了不被点本身遮挡住，需要对画标注的位置做个偏移
end
hold off;     
% title('ROC','fontsize',20);
% xlabel('False Positive Rate','fontsize',20);
% ylabel('Specificity','fontsize',20);
ylim(axes1,[0 100]);
xlim(axes1,[0 100]);

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16,'XGrid','on','YGrid','on');

annotation(figure1,'line',[0.131837307152875 0.904628330995792],...
    [0.111167300380228 0.918250950570342],'LineWidth',2,'LineStyle','--');
