function Fig3_Makhlin
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
clc;
close all


c=2*pi*1e-3;
wL=314; A=60; B=30; s0=1/2; s1=-s0;  k=1;
W0=sqrt((wL+s0*A)^2+(s0*B)^2)*c;
W1=sqrt((wL+s1*A)^2+(s1*B)^2)*c;
W=W0+W1;
tk=@(k) 4*pi*(2*k-1)/W;
%Optimize the time a bit
tt=linspace(tk(k)-0.2,tk(k)+0.2,1e4);

%t_CPMG= Get_Opt_Time("CPMG",wL,A,B,s0,s1,tt);
t_CPMG=tk(k);
t_UDD3= Get_Opt_Time("UDD3",wL,A,B,s0,s1,tt);
t_UDD4= Get_Opt_Time("UDD4",wL,A,B,s0,s1,tt);

%Get expected maxima
maxNumber=6;

Maxima_CPMG = Get_Iters_of_Maxima("CPMG",wL,A,B,s0,s1,t_CPMG,maxNumber);
Maxima_UDD3 = Get_Iters_of_Maxima("UDD3",wL,A,B,s0,s1,t_UDD3,maxNumber);
Maxima_UDD4 = Get_Iters_of_Maxima("UDD4",wL,A,B,s0,s1,t_UDD4,maxNumber);
All_Maxima  = {Maxima_CPMG,Maxima_UDD3,Maxima_UDD4};
%%===================== Plot the results =================================

Max_Iters    = max([Maxima_CPMG,Maxima_UDD3,Maxima_UDD4]);

Num_Subplots = 4;

Sequences = ["CPMG","UDD3","UDD4"];
Titles    = {'CPMG','UDD_3','UDD_4'};
cnt=0;
times     = {t_CPMG,t_UDD3,t_UDD4};
xlimits   = {[0,100],[0,300],[0,300]};
FNTsize   = 22;

for ii=1:Num_Subplots
    if ii~=2
        cnt=cnt+1;
        subplot(2,2,ii)
        
        [Ep,G1,G2]=Get_EpVsN(Sequences(cnt),wL,A,B,s0,s1,times{cnt},Max_Iters);
        
        plot(1:Max_Iters,[Ep;G1;G2],'linewidth',2)
        hold on
        for jj=1:maxNumber
            
            line([All_Maxima{cnt}(jj),All_Maxima{cnt}(jj)],[0,1],'color','k','linestyle','--','linewidth',2)
            
        end
        
        xlim(xlimits{cnt})
        ylim([0,3])
        xlabel('$N$','interpreter','latex')
        title(Titles{cnt})
        fig_defaults(FNTsize)
        
        if ii==1
           legend({'$\epsilon_p/\epsilon_p^*$','$G_1$','$G_2$'},'interpreter','latex'...
               ,'location','best','color','none','edgecolor','none','NumColumns',3) 
            
        end
        set(gca,'YMinorTick', 'off','XMinorTick', 'off')
    
    
    end


end




end


function t_Opt=Get_Opt_Time(Option,wL,A,B,s0,s1,tt)

Nnuc=1;
k=1; %doesn't matter since we provide time
N=1;

if strcmp(Option,"CPMG")

for ii=1:length(tt)
    
   Spin =  SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
   Spin =  Spin.CPMG(tt(ii));
   Spin =  Spin.Rot_Angles_And_Axes;
   n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
    
end

elseif strcmp(Option,"UDD3")
 
for ii=1:length(tt)
    
   Spin =  SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
   Spin =  Spin.UDD(3,tt(ii));
   Spin =  Spin.Rot_Angles_And_Axes;
   n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
    
end
    
    
    
elseif strcmp(Option,"UDD4")

for ii=1:length(tt)
    
   Spin =  SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
   Spin =  Spin.UDD(4,tt(ii));
   Spin =  Spin.Rot_Angles_And_Axes;
   n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
    
end
        
    

end


[~,indx]=min(n0n1,[],'all','linear');

t_Opt = tt(indx);



end


function Maxima_N=Get_Iters_of_Maxima(Option,wL,A,B,s0,s1,t,maxNumber)

k=1;
N=1;
Nnuc=1;

Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);

if strcmp(Option,"CPMG")
Spin = Spin.CPMG(t);
elseif strcmp(Option,"UDD3")
Spin = Spin.UDD(3,t);    
elseif strcmp(Option,"UDD4")
Spin = Spin.UDD(4,t);    
end

Spin    = Spin.Expected_Maxima(maxNumber);

Maxima_N = Spin.Nmax;




end

function [Ep,G1,G2]=Get_EpVsN(Option,wL,A,B,s0,s1,t,Max_Iters)

Nnuc=1;
k=1;

if strcmp(Option,"CPMG")
    
    for N=1:Max_Iters
        
        
        Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
        Spin = Spin.CPMG(t);
        Spin = Spin.Makhlin_Inv;
        Ep(N) = Spin.Ep/(2/9);
        G1(N) = real(Spin.Makhlin{1}{1});
        G2(N) = real(Spin.Makhlin{1}{2});
        
    end
    
elseif strcmp(Option,"UDD3")    
    
    for N=1:Max_Iters
        
        
        Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
        Spin = Spin.UDD(3,t);
        Spin = Spin.Makhlin_Inv;
        Ep(N) = Spin.Ep/(2/9);
        G1(N) = real(Spin.Makhlin{1}{1});
        G2(N) = real(Spin.Makhlin{1}{2});
        
    end
    
elseif strcmp(Option,"UDD4")    
    
    for N=1:Max_Iters
        
        
        Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
        Spin = Spin.UDD(4,t);
        Spin = Spin.Makhlin_Inv;
        Ep(N) = Spin.Ep/(2/9);
        G1(N) = real(Spin.Makhlin{1}{1});
        G2(N) = real(Spin.Makhlin{1}{2});
        
    end
    
    
end



end
