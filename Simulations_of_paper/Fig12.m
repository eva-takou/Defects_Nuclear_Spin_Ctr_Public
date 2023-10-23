function Fig12
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------



close all; clc;

s0=0; s1=-1; %Fig12 is for s0=0 and s1=-1, and Fig5 is for s0=1/2 and s1=-1/2

wL=314; At=60; Bt=30; c=2*pi*1e-3; k=1; N=1; Nnuc=1;

W0 =  sqrt((wL+s0*At)^2+(s0*Bt)^2)*c; 
W1 =  sqrt((wL+s1*At)^2+(s1*Bt)^2)*c;
W  =  W0+W1;
tk =  @(k) 4*pi*(2*k-1)/W;

Target_Spin = SubClass_U4Operations(wL,At,Bt,s0,s1,Nnuc,k,N);
Target_Spin = Target_Spin.CPMG(tk(k));
Target_Spin = Target_Spin.Expected_Maxima(10);
N           = Target_Spin.Nmax(2); %select for example the 2nd maximum

%===== Generate unwanted spins and find their one-tangle ==================

Au = 10:5:200; Bu = Au; 
numA = numel(Au); numB=numel(Bu);
 
parfor ii=1:numA
    for jj=1:numB    
        Unwanted_Spin = SubClass_U4Operations(wL,Au(ii),Bu(jj),s0,s1,Nnuc,k,N);
        Unwanted_Spin = Unwanted_Spin.CPMG(tk(k));
        Unwanted_Spin = Unwanted_Spin.Makhlin_Inv;
        Ep_Unw(ii,jj) = Unwanted_Spin.Ep/(2/9); %The nuclear tangle can be
        %evaluated either with entangling power of 2-qubit gate or
        %with one-tangle expression.   
    end
end

%=== Find which nuclei satisfy the condition of trivial Evol ==============

tol_Tangle  = 0.1;
tol_Trivial = 0.1;

INDs = find(Ep_Unw<tol_Tangle);
[rows,cols] = ind2sub(size(Ep_Unw),INDs);

A_red = Au(rows);
B_red = Bu(cols);

%Check if trivial evol is satisfied
w0 =@(A,B)  sqrt((wL+s0*A)^2+(s0*B)^2)*c; w1 =@(A,B)  sqrt((wL+s1*A)^2+(s1*B)^2)*c;

Trivial_1 =@(A,B) mod(tk(k)/4*w1(A,B)*c,2*pi);
Trivial_2 =@(A,B) mod(tk(k)/4*w0(A,B)*c,2*pi);
cnt=0;
for ii=1:length(Au)
    for jj=1:length(Bu)
        
        cond1=Trivial_1(Au(ii),Bu(jj));
        cond2=Trivial_2(Au(ii),Bu(jj));
        
        disp(['Cond1:',num2str(cond1)])
        disp(['Cond2:',num2str(cond2)])
        
      if  ( cond1 <tol_Trivial || cond2<tol_Trivial ) && Ep_Unw(ii,jj)< tol_Tangle 
          cnt=cnt+1;
          A_store(cnt)=Au(ii);
          B_store(cnt)=Bu(jj);
      else
          cnt=cnt+1;
          A_store(cnt)=nan;
          B_store(cnt)=nan;
      end
        
        
    end
end


%=== Optimize in k=[1,5] of target spin to minimize unwanted tangle =======

kmax=5;

Out=minimize_unwanted_Tangle(wL,s0,s1,Au,Bu,At,Bt,kmax,tk);


indxA=find(Au==60);
indxB=find(Bu==30);

Min_EPU = Out.EPU;
Min_EPU(indxA,indxB)=nan; %exclude the target spin


%=== Plot everything =====================================================
FNTsize=22;


subplot(2,2,1)

surf(Au,Bu,Ep_Unw.')
properties_of_surf_plot(FNTsize)
xlabel('A/2\pi (kHz)'); ylabel('B/2\pi (kHz)'); 
xlim([10,200]); ylim([10,200])
title('$\epsilon_{12|3}/\epsilon_p^*$','interpreter','latex')

hold on

scatter3(A_red,B_red,ones(1,length(A_red)),20,'filled','g')
hold on
scatter3(A_store,B_store,ones(1,length(A_store)),25,'r','linewidth',1.5)


subplot(2,2,2)

surf(Au,Bu,2/9*1e3*Min_EPU.')
text(200,250,'x10^{-3}','fontsize',FNTsize)
xlabel('A/2\pi (kHz)'); ylabel('B/2\pi (kHz)'); 
xlim([10,200]); ylim([10,200])
title('min($\epsilon_{12|3}$)','interpreter','latex')
properties_of_surf_plot(FNTsize)

subplot(2,2,3)
surf(Au,Bu,Out.NOpt.')
xlabel('A/2\pi (kHz)'); ylabel('B/2\pi (kHz)'); 
xlim([10,200]); ylim([10,200])
title('$N^*$','interpreter','latex')
properties_of_surf_plot(FNTsize)

subplot(2,2,4)
surf(Au,Bu,Out.kOpt.')
xlabel('A/2\pi (kHz)'); ylabel('B/2\pi (kHz)'); 
xlim([10,200]); ylim([10,200])
title('$k^*$','interpreter','latex')
properties_of_surf_plot(FNTsize)













end

function [Iters,Times,Res]=find_all_maxima_of_target_ep(wL,At,Bt,s0,s1,kmax,tk)

N=1; Nnuc=1; Maxima_Cutoff = 18;

for k=1:kmax
    
     Target_Spin = SubClass_U4Operations(wL,At,Bt,s0,s1,Nnuc,k,N);    
     Target_Spin = Target_Spin.CPMG(tk(k));
     Target_Spin = Target_Spin.Expected_Maxima0(Maxima_Cutoff);
     temp        = Target_Spin.Nmax;
     Nmax{k}     = temp(temp<301);
     disp(Nmax{k})
end

%Arrange the times and iters in terms of cases

Iters=[];Times=[]; Res=[];
for k=1:kmax
     L=length(Nmax{k});
     Iters=[Nmax{k},Iters];
     Times=[repmat(tk(k),[1,L]),Times];
     Res  =[repmat(k,[1,L]),Res]; 
end

end

function Out=minimize_unwanted_Tangle(wL,s0,s1,Au,Bu,At,Bt,kmax,tk)

[Iters,Times,Res] = find_all_maxima_of_target_ep(wL,At,Bt,s0,s1,kmax,tk);
Max_Cases         = length(Times);

Nnuc = 1; k = 1; numA = numel(Au);  numB = numel(Bu);  

parfor Cases=1:Max_Cases
    
    for ii=1:numA
        
       for jj=1:numB
           
           Unwanted_Spin       = SubClass_U4Operations(wL,Au(ii),Bu(jj),s0,s1,Nnuc,k,Iters(Cases));
           Unwanted_Spin       = Unwanted_Spin.CPMG(Times(Cases));
           Unwanted_Spin       = Unwanted_Spin.Makhlin_Inv;
           Ep_Unw(Cases,ii,jj) = Unwanted_Spin.Ep/(2/9);
           
       end
        
    end
    
end

%For each {Au,Bu} find which case is the optimal, yielding lowest unwanted
%one-tangle

for ii=1:numA
    
    for jj=1:numB
        
    temp = squeeze(Ep_Unw(:,ii,jj));
    
    [min_Val,indx]    = min(temp,[],'all','linear');
    t_Opt(ii,jj)      = Times(indx);
    k_Opt(ii,jj)      = Res(indx);
    N_Opt(ii,jj)      = Iters(indx);
    min_Ep_Unw(ii,jj) = min_Val;
    
        
    end
    
end

%Make sure as a sanity check that the target one-tangle remains maximal:


for ii=1:numA
    for jj=1:numB
        
      Target_Spin = SubClass_U4Operations(wL,At,Bt,s0,s1,Nnuc,k_Opt(ii,jj),N_Opt(ii,jj));
      Target_Spin = Target_Spin.CPMG(t_Opt(ii,jj));
      Target_Spin = Target_Spin.Makhlin_Inv;
      EpT         = Target_Spin.Ep/(2/9);
      if EpT<0.98
          error(['Found target one tangle <',num2str(0.98)])
          
      end
      
    end
end

Out.tOpt = t_Opt;
Out.kOpt = k_Opt;
Out.NOpt = N_Opt;
Out.EPU  = min_Ep_Unw;


end




function properties_of_surf_plot(FNTsize)

h.EdgeColor='none';    view(0,90)
shading interp
colormap(brewermap([],'PuBu')); cbh=colorbar('NorthOutside');  %set(cbh,'YTick',[0:0.5:1]); 

fig_defaults(FNTsize)
set(gca,'YMinorTick', 'off','XMinorTick', 'off')


end




