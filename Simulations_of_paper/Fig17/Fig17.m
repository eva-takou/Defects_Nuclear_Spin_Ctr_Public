%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

%%
clear;
close all;

tolHF=3;
nn=1e4; %5e5 was used in the paper

tic
OutCPMG1=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG2=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG3=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG4=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG5=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG6=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG7=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
OutCPMG8=Fidelity_reduction_subroutine(tolHF,nn,'CPMG');
toc
%%

OutUDD3_1=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_2=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_3=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_4=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_5=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_6=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_7=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');
OutUDD3_8=Fidelity_reduction_subroutine(tolHF,nn,'UDD3');

%%

OutUDD4_1=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_2=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_3=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_4=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_5=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_6=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_7=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');
OutUDD4_8=Fidelity_reduction_subroutine(tolHF,nn,'UDD4');


%% To take the average first of the random generations run the following:

clc;
%for ii=1:6 %random generations
    
    for kk=1:6 %1 through 6 spins in the bath
    
    
      for jj=1:length( OutCPMG1.yData{kk})    %intervals of one-tangles        
            
            OutCPMGav.ydata{kk}(jj) = mean([OutCPMG1.yData{kk}(jj),OutCPMG2.yData{kk}(jj),...
                OutCPMG3.yData{kk}(jj),OutCPMG4.yData{kk}(jj),...
                OutCPMG5.yData{kk}(jj),OutCPMG6.yData{kk}(jj),...
                OutCPMG7.yData{kk}(jj),OutCPMG8.yData{kk}(jj)] );
            
            
            OutUDD3av.ydata{kk}(jj) = mean([OutUDD3_1.yData{kk}(jj),OutUDD3_2.yData{kk}(jj),...
                OutUDD3_3.yData{kk}(jj),OutUDD3_4.yData{kk}(jj),...
                OutUDD3_5.yData{kk}(jj),OutUDD3_6.yData{kk}(jj),...
                OutUDD3_7.yData{kk}(jj),OutUDD3_8.yData{kk}(jj)] );
            
            
            OutUDD4av.ydata{kk}(jj) = mean([OutUDD4_1.yData{kk}(jj),OutUDD4_2.yData{kk}(jj),...
                OutUDD4_3.yData{kk}(jj),OutUDD4_4.yData{kk}(jj),...
                OutUDD4_5.yData{kk}(jj),OutUDD4_6.yData{kk}(jj),...
                OutUDD4_7.yData{kk}(jj),OutUDD4_8.yData{kk}(jj)] );
            
      end
    end
        
    for ii=1:length(OutCPMGav.ydata{1})
    
       AV_GPMG(ii) = mean([OutCPMGav.ydata{1}(ii),OutCPMGav.ydata{2}(ii)...
                         ,OutCPMGav.ydata{3}(ii),OutCPMGav.ydata{4}(ii),...
                          OutCPMGav.ydata{5}(ii),OutCPMGav.ydata{6}(ii)]);
                      
      
       AV_UDD3(ii) = mean([OutUDD3av.ydata{1}(ii),OutUDD3av.ydata{2}(ii)...
                         ,OutUDD3av.ydata{3}(ii),OutUDD3av.ydata{4}(ii),...
                          OutUDD3av.ydata{5}(ii),OutUDD3av.ydata{6}(ii)]);   
     
       AV_UDD4(ii) = mean([OutUDD4av.ydata{1}(ii),OutUDD4av.ydata{2}(ii)...
                         ,OutUDD4av.ydata{3}(ii),OutUDD4av.ydata{4}(ii),...
                          OutUDD4av.ydata{5}(ii),OutUDD4av.ydata{6}(ii)]);                   
    end
    
%Do also an average over all k in 1,6


    
%end

%% PLOT THE RESULTS (average)



close all
markers = {'o','d','*','^','>','+','s','p'};


listcolors ={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};


subplot(4,1,1)
for ii=1:6
    
    %semilogy(x_c,1-yData{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii})
    
    semilogy(OutCPMG2.xData,OutCPMGav.ydata{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii},...
        'markersize',8)
    hold on
end



ylabel('1-$F$','interpreter','latex')
xlabel('$\epsilon_{p|q}/\epsilon^*$','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',22,'fontname','Mircosoft Sans Serif')

hold on

for jj=1:1
errorbar(OutCPMG2.xData,OutCPMGav.ydata{jj},OutCPMG2.ErrBarStart,OutCPMG2.ErrBarEnd,'horizontal','linewidth',1,'color','k','linestyle','none')
hold on
end




set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
ylim([1e-2,1])
line([0.08,1],[1e-1,1e-1],'color','k','linewidth',0.7,'linestyle','--')
line([0.08,1],[5e-1,5e-1],'color','k','linewidth',0.7,'linestyle','--')


legend('1spin','2spins','3spins','4spins','5spins','6spins','','NumColumns',6,'color','none','edgecolor','none')
xlim([0.08,0.7])
set(gca,'xtick',0.1:0.1:0.7,'xticklabel',{'','','','','','',''})
set(gca,'ytick',[1e-2,1e-1,5e-1,1],'yticklabel',{'10^{-2}','10^{-1}','5x10^{-1}','10^{0}'})


subplot(4,1,2)
for ii=1:6
    
    %semilogy(x_c,1-yData{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii})
    
    semilogy(OutUDD3_1.xData,OutUDD3av.ydata{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii},...
        'markersize',8)
    hold on
end



ylabel('1-$F$','interpreter','latex')
xlabel('$\epsilon_{p|q}/\epsilon^*$','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',22,'fontname','Mircosoft Sans Serif')

hold on

for jj=1:1
errorbar(OutUDD3_1.xData,OutUDD3av.ydata{jj},OutUDD3_1.ErrBarStart,OutUDD3_1.ErrBarEnd,'horizontal','linewidth',1,'color','k','linestyle','none')
hold on
end



set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
ylim([1e-2,1])

line([0.08,1],[1e-1,1e-1],'color','k','linewidth',0.7,'linestyle','--')
line([0.08,1],[5e-1,5e-1],'color','k','linewidth',0.7,'linestyle','--')
%legend('1spin','2spins','3spins','4spins','5spins','6spins','','NumColumns',2,'color','none','edgecolor','none')
xlim([0.08,0.7])
set(gca,'xtick',0.1:0.1:0.7,'xticklabel',{'','','','','','',''})
set(gca,'ytick',[1e-2,1e-1,5e-1,1],'yticklabel',{'10^{-2}','10^{-1}','5x10^{-1}','10^{0}'})

subplot(4,1,3)
for ii=1:6
    
    %semilogy(x_c,1-yData{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii})
    
    semilogy(OutUDD4_1.xData,OutUDD4av.ydata{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii},...
        'markersize',8)
    hold on
end



ylabel('1-$F$','interpreter','latex')
xlabel('$\epsilon_{p|q}/\epsilon^*$','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',22,'fontname','Mircosoft Sans Serif')
xlim([0.08,0.7])
set(gca,'xtick',0.1:0.1:0.7,'xticklabel',{'','','','','','',''})

hold on

for jj=1:1
errorbar(OutUDD4_1.xData,OutUDD4av.ydata{jj},OutUDD4_1.ErrBarStart,OutUDD4_1.ErrBarEnd,'horizontal','linewidth',1,'color','k','linestyle','none')
hold on
end

line([0.08,1],[1e-1,1e-1],'color','k','linewidth',0.7,'linestyle','--')
line([0.08,1],[5e-1,5e-1],'color','k','linewidth',0.7,'linestyle','--')

ylim([1e-2,1])
%legend('1spin','2spins','3spins','4spins','5spins','6spins','','NumColumns',2,'color','none','edgecolor','none')

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)

set(gca,'ytick',[1e-2,1e-1,5e-1,1],'yticklabel',{'10^{-2}','10^{-1}','5x10^{-1}','10^{0}'})

subplot(4,1,4)

    
    %semilogy(x_c,1-yData{ii},'linewidth',2,'marker',markers{ii},'markerfacecolor',listcolors{ii},'markeredgecolor',listcolors{ii})
    plot(OutCPMG2.xData,AV_GPMG,'linewidth',2)
    hold on
    plot(OutUDD3_1.xData,AV_UDD3,'linewidth',2)
    hold on
    plot(OutUDD4_1.xData,AV_UDD4,'linewidth',2)
    
    
    legend('CPMG','UDD\fontsize{16}3','UDD\fontsize{16}4','color','none','edgecolor','none','NumColumns',3)
ylabel('1-$F$','interpreter','latex')
xlabel('$\epsilon_{p|q}/\epsilon^*$','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',22,'fontname','Mircosoft Sans Serif')
xlim([0.08,0.7])
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
%set(gca,'ytick',[1e-1,5e-2,1e-1,5e-1,1],'yticklabel',{'10^{-2}','5x10^{-2}','10^{-1}','5x10^{-1}','10^{0}'})
set(gca,'xtick',0.1:0.1:0.7)


