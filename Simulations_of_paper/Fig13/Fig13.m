%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%%

clc;
clearvars;
close all
%% Fix a constant B for unwanted spins

At = 60;  %Target spin's parallel HF
Bt = 30;  %Target spin's perpendicular HF
wL = 314; %Larmor frequency
s0 =  0;  %Electron spin projection
s1 = -1;  %Electron spin projection
%% Get the maxima of Ep for the target spin for Time_Max = 2 ms

B=[30,30.05,31,50]; %Unwanted spin's perpendicular HF
DeltaA=cell(1,length(B));
Ep_Opt=cell(1,length(B));

parfor jj=1:length(B)

    [DeltaA{jj},Ep_Opt{jj}]=subroutine_Fig13(wL,s0,s1,At,Bt,B(jj));

end


%% Plot the results

close all
plot(DeltaA{1},Ep_Opt{1},'linewidth',2.5,'marker','o')

hold on
plot(DeltaA{1},Ep_Opt{2},'linewidth',1.5,'linestyle','-.')

hold on
plot(DeltaA{1},Ep_Opt{3},'linewidth',1,'marker','*')

hold on
plot(DeltaA{1},Ep_Opt{4},'linewidth',1.5,'marker','.')


xlabel('$|\Delta A|$','interpreter','latex')
ylabel('$\epsilon_{p}^{unw.}/\epsilon_p^*$','interpreter','latex')

legend('$B=2\pi \cdot 30$ kHz','$B=2\pi \cdot 30.05$ kHz',...
    '$B=2\pi \cdot 31$ kHz','$B=2\pi \cdot 50$ kHz',...
    'interpreter','latex','color','none','edgecolor','none','location','best')

set(gca,'fontsize',25,'fontname','Microsoft Sans Serif')
set(gcf,'color','w')

axes('Position',[.7 .7 .2 .2])
box on
semilogy(DeltaA{1},Ep_Opt{1},'marker','o','linewidth',1.5)
set(gca,'ytick',[1e-5,1e-3,1e0],'yticklabels',{'10^{-5}','10^{-3}','10^0'},...
     'xtick',0:2.5:5,'xticklabels',{'0','2.5','5'})
%set(gca,'xtick',0:2.5:5,'xticklabels',{'0','2.5','5'})

xlim([0,5])
xlabel('$|\Delta A|$','interpreter','latex')
ylabel('$\epsilon_{p}^{unw.}/\epsilon_p^*$','interpreter','latex')
set(gca,'fontsize',18,'fontname','Microsoft Sans Serif')





