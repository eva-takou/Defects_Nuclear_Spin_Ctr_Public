function Fig2_Trivial_Evol
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

close all


listcolors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','#0048BA',...
    '#B284BE','#F19CBB','#9966CC','#3DDC84','#CD9575','#915C83','#FBCEB1','#8F9779','#27346F','#9F8170',...
    '#3D2B1F','#CAE00D','#3D0C02','#BFAFB2','#1F75FE','#0D98BA','#8A2BE2',...
    '#006A4E','#D891EF','#5F9EA0','#EFBBCC','#C95A49','#6D9BC3','#232B2B','#D2691E','#6F4E37','#FBEC5D','#483D8B'};

s0=1/2; s1=-1/2; wL=314;  MaxHF=300;

times_1 = [9.8 ,13 ,19, 22, 30, 30];
kappa_1 = [1   ,1  ,2,  2,  2,  3];


colors_1 = {listcolors{1},[],listcolors{3},listcolors{4},[],listcolors{6}};
colors_2 = {[],listcolors{2},[],[],listcolors{5},[]};

FntSize = 22;

subplot(1,2,1)

for jj=1:length(times_1)
    
   Get_Trivial_Evol(kappa_1(jj),times_1(jj),s0,s1,wL,MaxHF,colors_1{jj},colors_2{jj}) 
    hold on
end

xlabel('A/(2\pi) (kHz)')
ylabel('B/(2\pi) (kHz)')
title('S=1/2')
fig_defaults(FntSize)
set(gca,'XminorTick','off','YMinorTick','off','YGrid','off')



s0=3/2; s1=-1/2;

times_2 = [6 ,8 ,9.5 ,14 ,14, 30];
kappa_2 = [1 ,1 ,1   ,1  ,2,  2];

colors_1 = {listcolors{1},listcolors{2},listcolors{3},[],listcolors{5},[]};
colors_2 = {[],[],[],listcolors{4},[],listcolors{6}};


subplot(1,2,2)

for jj=1:length(times_2)
    
   Get_Trivial_Evol(kappa_2(jj),times_2(jj),s0,s1,wL,MaxHF,colors_1{jj},colors_2{jj}) 
    hold on
end
xlabel('A/(2\pi) (kHz)')
ylabel('B/(2\pi) (kHz)')
title('S=3/2')
fig_defaults(FntSize)
set(gca,'XminorTick','off','YMinorTick','off','YGrid','off')




end

function Get_Trivial_Evol(kappa,t,s0,s1,wL,MaxHF,color1,color2)

A = 2:1:200;
c = 2*pi*1e-3;

%find B's (1st condition)

cnt=0;

for ii=1:length(A)
        
       if  (8*kappa*pi/t)^2 > (wL*c+s0*A(ii)*c)^2 && 1/abs(s0)*sqrt(   (8*kappa*pi/t)^2 - (wL*c+s0*A(ii)*c)^2 )/c <= MaxHF
          
         cnt=cnt+1;
         Bsol(cnt)=  1/abs(s0)*sqrt(   (8*kappa*pi/(t))^2 - (wL*c+s0*A(ii)*c)^2 );  %B is in MHz
         BsolkHz(cnt)=  Bsol(cnt)/c;
         AsolkHz(cnt) = A(ii);
         
         
       end
       
    
end


if exist('AsolkHz')

x0 = -wL*c/s0;  y0 = 0;
R = kappa*8*pi/(t*s0);
th=linspace(0,2*pi,1e3);
x = R*cos(th)+x0;
y = R*sin(th)+y0;

hold on

plot(x/c,y/c,'linewidth',2,'color',color1,'linestyle','-')
hold on

scatter(AsolkHz,BsolkHz,20,'filled','markerfacecolor',color1,...
        'markeredgecolor',color1)

hold on


axis equal
xlim([0 MaxHF])
ylim([0 MaxHF])

end


%2nd condition:

cnt=0;
for ii=1:length(A)
        
       if  (8*kappa*pi/t)^2 > (wL*c+s1*A(ii)*c)^2 && 1/abs(s1)*sqrt(   (8*kappa*pi/t)^2 - (wL*c+s1*A(ii)*c)^2 )/c <= MaxHF
          
         cnt=cnt+1;
         Bsol2(cnt)=  1/abs(s1)*sqrt(   (8*kappa*pi/(t))^2 - (wL*c+s1*A(ii)*c)^2 );  %B is in MHz
         BsolkHz2(cnt)=  Bsol2(cnt)/c;
         AsolkHz2(cnt) = A(ii);
         
         
       end
       
    
end



if exist('AsolkHz2')

x0 = -wL*c/s1;  y0 = 0;
R = kappa*8*pi/(t*s1);
th=linspace(0,2*pi,1e3);
x = R*cos(th)+x0;
y = R*sin(th)+y0;

hold on

plot(x/c,y/c,'linewidth',2,'color',color2,'linestyle','-')
hold on

scatter(AsolkHz2,BsolkHz2,20,'filled','markerfacecolor',color2,...
        'markeredgecolor',color2)

hold on


axis equal
xlim([0 MaxHF])
ylim([0 MaxHF])

end


% if exist('AsolkHz')
%     
%     out.Acond1 = AsolkHz;
%     out.Bcond1 = BsolkHz;
% end
% 
% if exist('AsolkHz2')
%     
%     out.Acond2 = AsolkHz2;
%     out.Bcond2 = BsolkHz2;
% end



end