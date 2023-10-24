%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%%

clc; 
close all; 
clearvars;
%%
t = linspace(0,20,1e6);

A=70; B=10;  wL=314;  c=2*pi*1e-3;
s0=1/2;  s1=-1/2;  Nnuc=1;  k=1;  N=1;

kmax=20;  n=4;  Nmax=5000;

parfor ii=1:length(t)
    
   Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N); 
    
   Spin = Spin.UDD(n,t(ii));
    
   Spin = Spin.Rot_Angles_And_Axes;
   
   n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
    
end

W0 = sqrt((wL+s0*A)^2+(s0*B)^2)*c;
W1 = sqrt((wL+s1*A)^2+(s1*B)^2)*c;
W  = W0+W1;

tk_1  = @(k) 4*pi*(2*k-1)/W;   
tk_2  = @(k) 8*pi*(2*k-1)/W;  
tk_3  = @(k) 12*pi*(2*k-1)/W;  
tk_4  = @(k) 16*pi*(2*k-1)/W;  
tk_5  = @(k) 18*pi*(2*k-1)/W;  
tk_6  = @(k) 20*pi*(2*k-1)/W;  
tk_7  = @(k) 22*pi*(2*k-1)/W;  
tk_8  = @(k) 10*pi*(2*k-1)/W; 
tk_9  = @(k) 14*pi*(2*k-1)/W;  
tk_10 = @(k) 6*pi*(2*k-1)/W;   


plot(t,n0n1,'linewidth',2)

hold on

for k=1:kmax
    
   line([tk_1(k),tk_1(k)],[-1,1],'color','k','linestyle','--','linewidth',2) 
   hold on
   line([tk_2(k),tk_2(k)],[-1,1],'color','r','linestyle','--','linewidth',2) 
   hold on
   line([tk_3(k),tk_3(k)],[-1,1],'color','g','linestyle','--','linewidth',2) 
   hold on
   line([tk_4(k),tk_4(k)],[-1,1],'color','m','linestyle','--','linewidth',2) 
   hold on
   line([tk_5(k),tk_5(k)],[-1,1],'color','y','linestyle','--','linewidth',2) 
   hold on
   line([tk_6(k),tk_6(k)],[-1,1],'color','b','linestyle','--','linewidth',2) 
   hold on
   line([tk_7(k),tk_7(k)],[-1,1],'color','c','linestyle','--','linewidth',2) 
   hold on
   line([tk_8(k),tk_8(k)],[-1,1],'color',[0.4660 0.6740 0.1880],'linestyle','--','linewidth',2) 
   hold on
   
   line([tk_9(k),tk_9(k)],[-1,1],'color',[0.8500 0.3250 0.0980],'linestyle','--','linewidth',2) 
end

set(gcf,'color','w')
set(gca,'fontsize',20,'fontname','microsoft sans serif')
xlabel('$t$ ($\mu$ s)','interpreter','latex')
ylabel('$\textbf{n}_0 \cdot \textbf{n}_1$','interpreter','latex')

ylim([-1,1])
xlim([0,max(t)])


%What about "false" resonances? 

T_false_UDD4 = [tk_1(1),7.96076,11.144623,12.7372,14.3295,17.5136]; %T_false 
LT           = length(T_false_UDD4);

parfor ii=1:Nmax
    for kk=1:LT
    
    
   Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,ii); 
    
   Spin = Spin.UDD(n,T_false_UDD4(kk));
    
   Spin = Spin.Makhlin_Inv;
   Spin = Spin.Rot_Angles_And_Axes;
   
   phi0(ii,kk) = Spin.angles{1};
   phi1(ii,kk) = Spin.angles{2};
   Ep(ii,kk) = Spin.Ep;
   
    
    
    end

end

%%
figure (2)

for kk=1:length(T_false_UDD4)
plot([1:Nmax],Ep(:,kk)/(2/9),'linewidth',2)
hold on
legs{kk} = strcat(num2str(T_false_UDD4(kk)),'$\mu$ s');
end
legend(legs,'color','none','interpreter','latex','edgecolor','none')
set(gca,'fontsize',20,'fontname','microsoft sans serif')
set(gcf,'color','w')
xlabel('$N$','interpreter','latex')
ylabel('$\epsilon_p^{nuc.}/\epsilon_p^* $','interpreter','latex')

figure (4)

subplot(3,1,1)

for kk=1:2
plot([1:Nmax],phi0(:,kk)/pi,'linewidth',2)
hold on
legs{kk} = strcat(num2str(T_false_UDD4(kk)),'$\mu$s');
end

legend(legs,'color','none','interpreter','latex','edgecolor','none')
set(gca,'fontsize',20,'fontname','microsoft sans serif')
set(gcf,'color','w')
xlabel('$N$','interpreter','latex')
ylabel('$\phi_0/\pi$','interpreter','latex')
%xlim([0,1000])

subplot(3,1,2)
for kk=1:2
plot([1:Nmax],phi1(:,kk)/pi,'linewidth',2)
hold on
legs{kk} = strcat(num2str(T_false_UDD4(kk)),'$\mu$ s');
end

legend(legs,'color','none','interpreter','latex','edgecolor','none')
set(gca,'fontsize',20,'fontname','microsoft sans serif')
set(gcf,'color','w')
xlabel('$N$','interpreter','latex')
ylabel('$\phi_1/\pi$','interpreter','latex')



subplot(3,1,3)
for kk=1:2
plot([1:Nmax],abs(phi0(:,kk)/pi-phi1(:,kk)/pi),'linewidth',2)
hold on
legs{kk} = strcat(num2str(T_false_UDD4(kk)),'$\mu$ s');
end

legend(legs,'color','none','interpreter','latex','edgecolor','none')
set(gca,'fontsize',20,'fontname','microsoft sans serif')
set(gcf,'color','w')
xlabel('$N$','interpreter','latex')
ylabel('$|\Delta_\phi|/\pi$','interpreter','latex')


%xlim([0,1000])







