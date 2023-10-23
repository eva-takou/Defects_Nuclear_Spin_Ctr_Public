function Fig4_One_Tangles
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
close all;

phi0 = linspace(-pi,pi,1e3);
phi1 = phi0;

G1  = @(phi0,phi1,n0n1) (cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2))^2;

for ii=1:length(phi0)
    for jj=1:length(phi1)
        G1_plus1(ii,jj)  = G1(phi0(ii),phi1(jj),1);
        G1_minus1(ii,jj) = G1(phi0(ii),phi1(jj),-1);
        
        if G1_plus1(ii,jj)<1e-5
           G1_plus1(ii,jj)=0; 
        end
        
        if G1_minus1(ii,jj)<1e-5
           G1_minus1(ii,jj)=0; 
        end
        
    end
end

%==== Get the bound (not-tight) of U(n) ===================================

Empty=SubClass_Ent_and_Fid;

for n=3:10
    
    ep_bound(n-2)   = Empty.one_tangles_General_U_Max_Bound(n).Max_Bound;
    ep_el_max(n-2)  = 1/3 -1/3^n;
    ep_nuc_max(n-2) = 2/9;
    
end

%=== Find max one-tangle for U(n) =========================================

Trials=3; %100 used in the paper
tic
Ep_4=Get_Ep_Random(4,Trials);
Ep_5=Get_Ep_Random(5,Trials);
toc

Trials=1; %in the paper used 5 trials
tic
Ep_6=Get_Ep_Random(6,Trials);
toc


Ep_Random_Max = [Ep_4,Ep_5,Ep_6];

Ep_Random_Max=[epAME,Ep_Random_Max];

FNTsize=20;
titles={'$\textbf{n}_0\cdot\textbf{n}_1=1$','$\textbf{n}_0\cdot\textbf{n}_1=-1$'};

subplot(2,2,1)

h=surf(phi0,phi1,G1_plus1.');  

properties_of_surf_plot(FNTsize,titles{1})

subplot(2,2,2)

h=surf(phi0,phi1,G1_minus1.');  
properties_of_surf_plot(FNTsize,titles{2})

subplot(2,2,3)

plot(3:10,ep_el_max,'marker','s','markerfacecolor','k','linewidth',2,'markersize',9)
hold on
plot(3:10,ep_nuc_max,'marker','d','markerfacecolor','k','linewidth',2,'markersize',9)
hold on
plot(3:10,ep_bound,'marker','o','markerfacecolor','k','linewidth',2,'markersize',9)
hold on
plot(3:6,Ep_Random_Max,'marker','^','markerfacecolor','k','linewidth',2,'markersize',9)

xlabel('Number of qubits (n)')
ylabel('$\epsilon_{p|q}^*$','interpreter','latex')
fig_defaults(FNTsize)
line([3,10],[1/3,1/3],'color','k','linewidth',2,'linestyle','--')
xlim([3,10])
ylim([.2,.5])
legend({'electron','nuclear spin','$U(n)$','random $U(n)$'},...
    'interpreter','latex','location','best','color','none',...
    'edgecolor','none','NumColumns',2)
set(gca,'Ytick',[0.2,0.3,0.4,0.5],'Yticklabels',{'0.2','0.3','0.4','0.5'})

end


function properties_of_surf_plot(FNTsize,TITLE)

h.EdgeColor='none';    view(0,90)
colormap(brewermap([],'PuBu')); cbh=colorbar;  set(cbh,'YTick',[0:0.5:1]); 

shading interp
xlabel('\phi_0 (rad)');  ylabel('\phi_1 (rad)');

xlim([-pi,pi]); ylim([-pi,pi]);

set(gca,'Xticklabel',{'-\pi','0','\pi'},'Xtick',[-pi,0,pi])
set(gca,'Yticklabel',{'-\pi','0','\pi'},'Ytick',[-pi,0,pi])

title(TITLE,'interpreter','latex')

fig_defaults(FNTsize)

set(gca,'YMinorTick', 'off','XMinorTick', 'off')


end

function psi=AME62

ket0=[1;0];
ket1=[0;1];

kets={ket0,ket1};
psi=zeros(2^6,1);

for x=0:1
    for y=0:1
        for z=0:1
            
            psi =psi+ kron( kets{x+1}, kron(kets{y+1},kron(kets{z+1}  ,    GHZijk(x,y,z))));
            
        end
    end
end

end


function out=GHZijk(ii,jj,kk)

X=[0 1 ; 1 0];
Z=[1 0 ; 0 -1];

S={X,Z};
ket0=[1;0];
ket1=[0;1];

ket000=kron(ket0,kron(ket0,ket0));
ket111=kron(ket1,kron(ket1,ket1));

GHZ=1/sqrt(2)*(ket000+ket111);

if ii==jj && jj==kk
    
    c=-1;
    
else
    c=1;
    
    
end
    
out=c*kron(S{ii+1},kron(S{jj+1},S{kk+1}))*GHZ;


end


function Ep=epAME


ket0=[1;0];
ket1=[0;1];

kets={ket0,ket1};
psi=zeros(2^6,1);



for i1=0:1
    for i2=0:1
        for i3=0:1
            for j1=0:1
                for j2=0:1
                    for j3=0:1
    
    Uelem{i1+1,i2+1,i3+1,j1+1,j2+1,j3+1}=sqrt(2^3)*...
        (...
        kron(kets{i1+1}, kron(kets{i2+1},kron(kets{i3+1},...
        kron(kets{j1+1}, kron(kets{j2+1},kets{j3+1})))))...
        )'*AME62;
                    end
                end
            end
        end
    end
end
Uket=0;
for i1=0:1
    for i2=0:1
        for i3=0:1
            for j1=0:1
                for j2=0:1
                    for j3=0:1
    
    Uket=Uket+1/sqrt(2^3)*Uelem{i1+1,i2+1,i3+1,j1+1,j2+1,j3+1}*...
    kron(kets{i1+1}, kron(kets{i2+1},kron(kets{i3+1},...
        kron(kets{j1+1}, kron(kets{j2+1},kets{j3+1})))));
                    end
                end
            end
        end
    end
end



Uoper=reshape(Uket,[8,8]);

test=SubClass_Ent_and_Fid;

test.Uval=Uoper;

test=test.Ent_Power_Zanardi(1);
Ep=test.One_Tangle_Zanardi;


end


function Ep=Get_Ep_Random(n,Trials)


parfor ii=1:Trials
    
   A                 = rand(2^n);
   U                 = expm(-1i*(A+A'));
   Temp              = SubClass_Ent_and_Fid;
   Temp.Uval         = U;
   Temp              = Temp.one_tangles_General_U_byparts;
   Ep_Random(ii)     = max(Temp.one_tangles{1}); 
    
end

Ep=max(Ep_Random(:));


end




