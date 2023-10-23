function Fig1_n0n1VsN
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
close all
wL=314; A=80; B=25; s0=1/2; s1=-1/2; c    = 2*pi*1e-3;  kmax=4; Nmax = 400;

W0   = sqrt((wL+s0*A)^2+(s0*B)^2)*c;
W1   = sqrt((wL+s1*A)^2+(s1*B)^2)*c;
W    = W0+W1;

tk_one = @(k) 4*pi*(2*k-1)/W;
tk_two = @(k) 8*pi*(2*k-1)/W;


t_CPMG     = Get_Opt_Time("CPMG",tk_one,wL,A,B,s0,s1,kmax);
t_UDD3     = Get_Opt_Time("UDD3",tk_one,wL,A,B,s0,s1,kmax);
t1_UDD4    = Get_Opt_Time("UDD4",tk_one,wL,A,B,s0,s1,kmax);
t2_UDD4    = Get_Opt_Time("UDD4",tk_two,wL,A,B,s0,s1,kmax);

n0n1_CPMG   = Get_n0n1VsN("CPMG",wL,A,B,s0,s1,t_CPMG,kmax,Nmax);
n0n1_UDD3   = Get_n0n1VsN("UDD3",wL,A,B,s0,s1,t_UDD3,kmax,Nmax);
n0n1_UDD4_1 = Get_n0n1VsN("UDD4",wL,A,B,s0,s1,t1_UDD4,kmax,Nmax);
n0n1_UDD4_2 = Get_n0n1VsN("UDD4",wL,A,B,s0,s1,t2_UDD4,kmax,Nmax);


n0n1        = {n0n1_CPMG,n0n1_UDD3,n0n1_UDD4_1,n0n1_UDD4_2};


%======== Plot the results ================================================

Markers  = {'o','s','d','>'};
MarkerSz = [12,10,8,6];
Colors   = {[0 0.4470 0.7410],...
            [0.8500 0.3250 0.0980],...
            [0.9290 0.6940 0.1250],...
            [0.4940 0.1840 0.5560]};
max_Res  = [4,4,2,2]; %how many resonances to plot in each subplot

FntSize    = 22;
MarkerSz   = [8,8,8,8];
steps_of_N = [10,25,40,60];

for ii=1:length(n0n1)
    
    subplot(2,2,ii)
    
    if ii~=3 && ii~=4
    
        for jj=1:max_Res(ii)
            plot(1:steps_of_N(jj):Nmax,n0n1{ii}(1:steps_of_N(jj):end,jj),'linewidth',1,...
                'marker',Markers{jj},'markersize',MarkerSz(jj),...
                'markerfacecolor',Colors{jj});
            hold on
        end
    else
        for jj=1:max_Res(ii)
            plot(1:Nmax,n0n1{ii}(:,jj),'linewidth',1,...
                'marker',Markers{jj},'markersize',MarkerSz(jj),...
                'markerfacecolor',Colors{jj});
            hold on
        end
        
    end
    
    xlabel('$N$','interpreter','latex')
    ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
    ylim([-1,1])
    xlim([0,Nmax])
    fig_defaults(FntSize)
    
    if ii==1
    legend({'$k=1$','$k=2$','$k=3$','$k=4$'},'interpreter','latex','location','best',...
        'color','none','edgecolor','none','NumColumns',4)
    end
    
    set(gca,'Ytick',[-1,0,1],'Yticklabel',{'-1','0','1'},'YMinorTick', 'off'...
             ,'XMinorTick', 'off')
    
    
end


end

function tOpt = Get_Opt_Time(Option,tk,wL,A,B,s0,s1,kmax)


Nnuc = 1;

max_T_Values= 1e4;

tt=zeros(kmax,max_T_Values);


for kk=1:kmax
tt(kk,:) = linspace(tk(kk)-0.15,tk(kk)+0.15,max_T_Values);
end

n0n1 = zeros(max_T_Values,kmax);

if strcmp(Option,"CPMG")

parfor ii=1:length(tt)
    
        for kk=1:kmax
                    
            Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,1);   
            Spin = Spin.CPMG(tt(kk,ii));
            Spin = Spin.Rot_Angles_And_Axes;

            n0n1(ii,kk) = dot(Spin.axes{1},Spin.axes{2});

        end

end
    
elseif strcmp(Option,"UDD3")    
    
parfor ii=1:length(tt)
    
        for kk=1:kmax
                    
            Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,1);   
            Spin = Spin.UDD(3,tt(kk,ii));
            Spin = Spin.Rot_Angles_And_Axes;

            n0n1(ii,kk) = dot(Spin.axes{1},Spin.axes{2});

        end

end

elseif strcmp(Option,"UDD4")
 
parfor ii=1:length(tt)
    
        for kk=1:kmax
                    
            Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,1);   
            Spin = Spin.UDD(4,tt(kk,ii));
            Spin = Spin.Rot_Angles_And_Axes;

            n0n1(ii,kk) = dot(Spin.axes{1},Spin.axes{2});

        end

end
    
    
end




tOpt = zeros(1,kmax);
for kk=1:kmax
    
   [~,indx]=min(n0n1(:,kk),[],'all','linear') ;
   tOpt(kk) = tt(kk,indx);
    
end



end


function n0n1 = Get_n0n1VsN(Option,wL,A,B,s0,s1,times,kmax,Nmax)

n0n1=zeros(Nmax,kmax);
Nnuc=1;

if strcmp(Option,"CPMG")

parfor N=1:Nmax
    
   for kk=1:kmax
       
       Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,N);
       Spin = Spin.CPMG(times(kk));
       Spin = Spin.Rot_Angles_And_Axes;
       n0n1(N,kk) = dot(Spin.axes{1},Spin.axes{2});
       
       
   end
    
    
    
end
    
elseif strcmp(Option,"UDD3")   
    
parfor N=1:Nmax
    
   for kk=1:kmax
       
       Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,N);
       Spin = Spin.UDD(3,times(kk));
       Spin = Spin.Rot_Angles_And_Axes;
       n0n1(N,kk) = dot(Spin.axes{1},Spin.axes{2});
       
       
   end
    
    
    
end
     
elseif strcmp(Option,"UDD4")    
    
parfor N=1:Nmax
    
   for kk=1:kmax
       
       Spin = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,kk,N);
       Spin = Spin.UDD(4,times(kk));
       Spin = Spin.Rot_Angles_And_Axes;
       n0n1(N,kk) = dot(Spin.axes{1},Spin.axes{2});
       
       
   end
    
    
    
end
    
    
    
end



end





