function  Output=Plot_ModRes_Plot_V1(ModRes,ParaSet)
% code for studies entitled:
% "Extended SEIR model for death and cure population of COVID-19 in China 
% under public intervention and clinical progress"
% Su Feng, Yuan Peijiang, Li Jianmin
% modified in 2020/03/29

Output=[];
N=ParaSet.T;
% %% confirmed
ci=1;
Output{ci,1}='Newly confirm';
Output{ci,2}=[(1:N)'  ModRes.DeltaI];

ci=2;
Output{ci,1}='Cummulative confirm';
Output{ci,2}=[(1:N)'  ModRes.CumsumI];

% %% death
TR2=ModRes.Death;
TR2(1:ParaSet.Death(1)-1)=0;
ci=3;
Output{ci,1}='Cummulative death';
Output{ci,2}=[(1:N)' TR2];

% %% cure
TR2=ModRes.Cure;
ci=4;
Output{ci,1}='Cummulative cure';
Output{ci,2}=[(1:N)' TR2];

% visualization
figure;
Nplot=size(Output,1);
XTickInd=1:30:N;
XTickLabel=ModRes.TimeSeq(XTickInd);
for xi=1:length(XTickLabel)
    XTickLabel{xi}=strrep(XTickLabel{xi},'2020/','');
end
for ci=1:Nplot
subplot(1,Nplot,ci); hold on; grid on
plot(Output{ci,2}(:,2),'k')
title(Output{ci,1},'FontName','Deng')
set(gca,'XTick',XTickInd)
set(gca,'xtickLabel',XTickLabel)

set(gcf,'Position', [37 387 1154 286])
end

end

