%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------
%%
clc;
clear;
close all;

A  = 30;
B  = 60;
wL = 314;
t  = 3;
s0=0;
s1=-1;
Nmax=300;


G1=@(phi0,phi1,n0n1) (cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2))^2;
G2=@(phi0,phi1,n0n1) 1+n0n1*sin(phi0)*sin(phi1)+...
                     2*( cos(phi0/2)^2*cos(phi1/2)^2+n0n1^2*sin(phi0/2)^2*sin(phi1/2)^2);

for jj=1:Nmax
    
    Spin=SubClass_U4Operations(wL,A,B,s0,s1,1,1,jj);
    Spin=Spin.CPMG(t);
    Spin=Spin.Makhlin_Inv;
    Spin=Spin.Rot_Angles_And_Axes;
    
    phi0 = Spin.angles{1};
    phi1 = Spin.angles{1};
    n0   = Spin.axes{1};
    n1   = Spin.axes{2};
    
    G1_numer(jj) = Spin.Makhlin{1}{1};
    G2_numer(jj) = Spin.Makhlin{1}{2};
    
    G1_analy(jj) = G1(phi0,phi1,dot(n0,n1));
    G2_analy(jj) = G2(phi0,phi1,dot(n0,n1));
    
end


figure(1)
subplot(2,1,1)
plot(1:Nmax,G1_numer,'linewidth',2)
hold on
plot(1:Nmax,G1_analy,'linewidth',2,'linestyle','--')
set(gcf,'color','w')
set(gca,'fontsize',20,'fontname','Microsoft Sans Serif')
legend('Numerical','Analytical')
ylabel('$G_1$','interpreter','latex')
xlabel('$N$','interpreter','latex')

subplot(2,1,2)
plot(1:Nmax,G2_numer,'linewidth',2)
hold on
plot(1:Nmax,G2_analy,'linewidth',2,'linestyle','--')
set(gcf,'color','w')
set(gca,'fontsize',20,'fontname','Microsoft Sans Serif')
legend('Numerical','Analytical')
ylabel('$G_2$','interpreter','latex')
xlabel('$N$','interpreter','latex')

