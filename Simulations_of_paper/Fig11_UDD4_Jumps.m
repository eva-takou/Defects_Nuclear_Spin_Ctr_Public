function Fig11_UDD4_Jumps
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
close all;

Nnuc=1;
wL=314;
A=60;
B=30;
c=2*pi*1e-3;
s0=1/2; s1=-1/2;

W0=sqrt((wL+s0*A)^2+(s0*B)^2)*c;
W1=sqrt((wL+s1*A)^2+(s1*B)^2)*c;

W=W0+W1;

k=1; N=1;
Nmax=400;


tk = @(k) 4*pi*(2*k-1)/W;

%Optimize the time a bit:
tt = linspace(tk(k)-0.15,tk(k)+0.15,1e4);

for ii=1:length(tt)
    
   Spin=SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
   Spin=Spin.UDD(4,tt(ii));
   Spin=Spin.Rot_Angles_And_Axes;
   n0n1(ii)=dot(Spin.axes{1},Spin.axes{2});
   
end

[~,indx]=min(n0n1);

tOpt = tt(indx);
clear n0n1
disp(['t_Opt:',num2str(tOpt)])




for N=1:Nmax

Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N);
Spin = Spin.UDD(4,tOpt);
Spin = Spin.Rot_Angles_And_Axes;

DeltaPhi(N) = abs(Spin.angles{1}-Spin.angles{2});

Spin = Spin.Makhlin_Inv;
Ep(N)  = Spin.Ep/(2/9);

n0n1(N) = dot(Spin.axes{1},Spin.axes{2});

end

plot(1:Nmax,[n0n1;DeltaPhi;Ep],'linewidth',2)
line([0,Nmax],[0,0],'color','k','linewidth',2,'linestyle','--')
FntSize=22;
legend({'$\textbf{n}_0\cdot\textbf{n}_1$','$\Delta \phi$','$\epsilon_p/\epsilon_{p}^*$'},'interpreter','latex',...
    'color','none','edgecolor','none','NumColumns',3,'location','best')
fig_defaults(FntSize)

set(gca,'YMinorTick','off','XMinorTick','off')
xlabel('$N$','interpreter','latex')
ylim([-1,1])
xlim([0,Nmax])




end
