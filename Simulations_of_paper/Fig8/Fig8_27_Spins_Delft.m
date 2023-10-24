function Out=Fig8_27_Spins_Delft(Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

clc;

switch Option
    
    case 'Load'
        
        load('27SpinsData.mat')
        

    case 'Run_Again'
        
        Out=Get_Optimal_Cases;
end



%Collect the data into arrays:
cnt=0;
for k=1:5
    
    for Spin_Indx=1:27
        cnt=cnt+1;
        temp=Out{Spin_Indx,k};
        if isstruct(temp)
            
            EpT_mean(cnt) = mean(temp.EpT);
            EpU_mean(cnt) = mean(temp.EpUnw);
            Target_Spins(cnt) = length(temp.TargetNuc);
            Gate_Time(cnt)    = temp.Gate_Time;
            N(cnt)            = temp.N;
            Infid(cnt)        = temp.Infid;
            
        else
            EpT_mean(cnt)=nan;
            EpU_mean(cnt)=nan;
            Target_Spins(cnt)=nan;
            Gate_Time(cnt)=nan;
            N(cnt)=nan;
            Infid(cnt)=nan;
        end
        
    end
end

yData={EpT_mean,EpU_mean,Target_Spins,N,Gate_Time,Infid};

close all;
rows  = 3;
cols  = 2;
Markers  = {'o','s','v','>','<'};
FntSize  = 22;

ylims={[0.85,1],[0,0.04],[2,7],[0,400],[0,1500],[0,0.2]};
ylabels={'$\bar{\epsilon}_p^T$','$\bar{\epsilon}_p^U$',...
    '# of nuclei','# of pulses','Gate time (\mus)','$1-F$'};
interp = {'latex','latex','none','none','tex','latex'};


figure(1)

for ii=1:length(yData)
   
    subplot(rows,cols,ii)
    
    for k=1:5
        cnst = (k-1)*27;
        stem(1+cnst:27+cnst,yData{ii}(1+cnst:27+cnst),"filled",...
            'marker',Markers{k},'linewidth',2,'markersize',8)
        hold on
        
    
    end
    
        fig_defaults(FntSize)
        xlabel('Case #')
        ylabel(ylabels{ii},'interpreter',interp{ii})
        ylim(ylims{ii})
        xlim([1,135])
        if ii==1
            
            legend({'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$'},...
                'location','best','interpreter','latex','color','none',...
                'edgecolor','none','NumColumns',5)
            
        end
       
    
end






end


function Out=Get_Optimal_Cases

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,1}=Optimize_Case_27_spin_Register(Spin_Indx,1)
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,2}=Optimize_Case_27_spin_Register(Spin_Indx,2)
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,3}=Optimize_Case_27_spin_Register(Spin_Indx,3)
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,4}=Optimize_Case_27_spin_Register(Spin_Indx,4)
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,5}=Optimize_Case_27_spin_Register(Spin_Indx,5)
   
end
toc



end






