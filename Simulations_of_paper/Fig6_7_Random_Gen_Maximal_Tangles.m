function Fig6_7_Random_Gen_Maximal_Tangles(Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

%Option can be Random to start a random generation of nuclei anew
%    or can be Data_From_Paper to reproduce Fig6-Fig7.

clc; close all;    s0=1/2; s1=-1/2;   wL=314;

switch Option
    
    case 'Random'
      
      disp('============= CPMG, k=1 ================')  
        
      k=1; tol_HF=25; Nnuc=2e3; 
      OUT_CPMG{1}=Maximize_Many_Tangles('CPMG',Nnuc,tol_HF,s0,s1,k,wL);
        
      disp('============= CPMG, k=2 ================')  
      
      k=2; tol_HF=20; Nnuc=2e3;
      OUT_CPMG{2}=Maximize_Many_Tangles('CPMG',Nnuc,tol_HF,s0,s1,k,wL);
      
      disp('============= UDD3, k=1 ================')  
      
      k=1; tol_HF=11; Nnuc=2e4;
      OUT_UDD3{1}=Maximize_Many_Tangles('UDD3',Nnuc,tol_HF,s0,s1,k,wL);
        
        
      disp('============= UDD3, k=3 ================')  
      
      k=3; tol_HF=11; Nnuc=2e4;
      OUT_UDD3{2}=Maximize_Many_Tangles('UDD3',Nnuc,tol_HF,s0,s1,k,wL);

      disp('============= UDD4, k=1 ================')  
      
      k=1; tol_HF=13; Nnuc=1e4;
      OUT_UDD4{1}=Maximize_Many_Tangles('UDD4',Nnuc,tol_HF,s0,s1,k,wL);
        
        
      disp('============= UDD4, k=2 ================')  
      
       k=2; tol_HF=13; Nnuc=1e4;
       OUT_UDD4{2}=Maximize_Many_Tangles('UDD4',Nnuc,tol_HF,s0,s1,k,wL);
      
       %=========== Gate errors ==========================================
       
        N=OUT_CPMG{2}.N; T=OUT_CPMG{2}.T; Ntarget=length(OUT_CPMG{2}.A);
        P_CPMG=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'CPMG'); %choose 2nd case, i.e. k=2.

        N=OUT_UDD3{2}.N; T=OUT_UDD3{2}.T; Ntarget=length(OUT_UDD3{2}.A);
        P_UDD3=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'UDD3'); %choose 2nd case, i.e. k=2.


        N=OUT_UDD4{2}.N; T=OUT_UDD4{2}.T; Ntarget=length(OUT_UDD4{2}.A);
        P_UDD4=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'UDD4'); %choose 2nd case, i.e. k=2.
        
    case 'Data_From_Paper'
        
        OUT_CPMG{1}=get_data_from_paper_CPMG('1',wL,s0,s1);
        OUT_CPMG{2}=get_data_from_paper_CPMG('2',wL,s0,s1);
        OUT_UDD3{1}=get_data_from_paper_UDD3('1',wL,s0,s1);
        OUT_UDD3{2}=get_data_from_paper_UDD3('2',wL,s0,s1);
        OUT_UDD4{1}=get_data_from_paper_UDD4('1',wL,s0,s1);
        OUT_UDD4{2}=get_data_from_paper_UDD4('2',wL,s0,s1);
        
        N=OUT_CPMG{2}.N; T=OUT_CPMG{2}.T; Ntarget=length(OUT_CPMG{2}.A);
        P_CPMG=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'CPMG'); %choose 2nd case, i.e. k=2.

        N=OUT_UDD3{2}.N; T=OUT_UDD3{2}.T; Ntarget=length(OUT_UDD3{2}.A);
        P_UDD3=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'UDD3'); %choose 2nd case, i.e. k=2.


        N=OUT_UDD4{2}.N; T=OUT_UDD4{2}.T; Ntarget=length(OUT_UDD4{2}.A);
        P_UDD4=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,'UDD4'); %choose 2nd case, i.e. k=2.
        
end

%===== Fig. 6 of Paper ====================================================

subplots_of_CPMG(OUT_CPMG,1)
subplots_of_UDD3(OUT_UDD3,2)
subplots_of_UDD4(OUT_UDD4,3)

%===== Fig. 7 of Paper ====================================================
%Note I am starting a new random generation of the unwanted spin bath.
%Expect different Fig.7, when 'Data_From_Paper' is selected.

%===== Plot the results: =================================================
figure(4)

subplot(2,2,1)
loglogplot(P_CPMG)
title('CPMG')

subplot(2,2,2)
kappa=1;
loglogTrivial(length(OUT_CPMG{2}.A),s0,s1,OUT_CPMG{2}.N,OUT_CPMG{2}.T,wL,kappa)
title('CPMG')

subplot(2,2,3)
loglogplot(P_UDD3)
title('UDD_3')

subplot(2,2,4)
loglogplot(P_UDD4)
title('UDD_4')
xlim([1e-3,1])

end

function loglogplot(S)

FntSize=22;
markers = {'o','d','*','^','>','+'};
listcolors ={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};

yData = S.Infid;

x_c   = S.center;
for kk=1:length(yData)
ErrBarStart(kk) =  S.ErrBar_start{kk};
ErrBarEnd(kk)   =  S.ErrBar_end{kk};
end   
       
for kk=1:length(yData)
   if ~isnan(yData{kk})
       Infid_1Spin(kk) =  yData{kk}(1);
       Infid_2Spin(kk) =  yData{kk}(2);
       Infid_3Spin(kk) =  yData{kk}(3);
       Infid_4Spin(kk) =  yData{kk}(4);
       Infid_5Spin(kk) =  yData{kk}(5);
       Infid_6Spin(kk) =  yData{kk}(6);
       
   else
       Infid_1Spin(kk) =  nan;
       Infid_2Spin(kk) =  nan;
       Infid_3Spin(kk) =  nan;
       Infid_4Spin(kk) =  nan;
       Infid_5Spin(kk) =  nan;
       Infid_6Spin(kk) =  nan;
       
       
   end
end

Infid=[Infid_1Spin;Infid_2Spin;Infid_3Spin;Infid_4Spin;Infid_5Spin;Infid_6Spin];

[row,~]=size(Infid);

for jj=1:row
    
   loglog(x_c,Infid(jj,:),'linewidth',2,'marker',markers{jj},'markerfacecolor',...
       listcolors{jj},'markeredgecolor',listcolors{jj},'color',listcolors{jj}) 
    hold on
    

end

line([1e-3,1e-3],[min(Infid(:)),max(Infid(:))],'color','k','linewidth',1,'linestyle','--')
line([min(x_c),max(x_c)],[1e-3,1e-3],'color','k','linewidth',1,'linestyle','--')
line([min(x_c),max(x_c)],[1e-1,1e-1],'color','k','linewidth',1,'linestyle','--')

hold on

    
errorbar(x_c,Infid(1,:),ErrBarStart,ErrBarEnd,'horizontal','linewidth',1,'color','k','linestyle','none')

  

fig_defaults(FntSize)
xlabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
ylabel('$1-F$','interpreter','latex')
legend({'1 spin','2 spins','3 spins','4 spins','5 spins','6 spins'},'color','none',...
    'edgecolor','none','location','best','NumColumns',1)

set(gca,'XMinorTick','on','YMinorTick','on')
end

function loglogTrivial(Ntarget,s0,s1,N,t,wL,kappa)

FntSize=22;
markers = {'o','d','*','^','>','+'};

Out=Get_Data_Trivial_Evol(kappa,t,N,s0,s1,wL,Ntarget);

for jj=1:length(Out.Infid)
loglog(Out.Ep(jj),Out.Infid(jj),'linewidth',2,'marker',markers{jj},'markerfacecolor',...
       'k','markeredgecolor','k','color','k','markersize',10,'linestyle','none') 
   hold on
end

fig_defaults(FntSize)
xlabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
ylabel('$1-F$','interpreter','latex')
legend({'1 spin','2 spins','3 spins','4 spins','5 spins','6 spins'},'color','none',...
    'edgecolor','none','location','best','NumColumns',1)

set(gca,'XMinorTick','on','YMinorTick','on')

end


function listcolors=COLORS

listcolors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','#0048BA',...
    '#B284BE','#F19CBB','#3DDC84','#CD9575','#915C83','#FBCEB1','#8F9779','#27346F','#9F8170',...
    '#3D2B1F','#CAE00D','#3D0C02','#BFAFB2','#1F75FE','#0D98BA','#8A2BE2',...
    '#006A4E','#D891EF','#5F9EA0','#EFBBCC','#C95A49','#6D9BC3','#232B2B','#D2691E','#6F4E37','#FBEC5D','#483D8B'};


end

function subplots_of_CPMG(S,fignum)
%S is the input struct
clrs = COLORS;

FntSize=20;

figure(fignum)

subplot(3,2,1) %(a)

for ii=1:length(S{1}.Ep)
stem(S{1}.phi0(ii)/(pi/2),S{1}.Ep(ii),"filled",'linewidth',2,...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
end
xlabel('$\phi_0/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.7,1])

subplot(3,2,2) %(b)

for ii=1:length(S{2}.Ep)
stem(S{2}.phi0(ii)/(pi/2),S{2}.Ep(ii),"filled",'linewidth',2,...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
end

xlabel('$\phi_0/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.7,1])

subplot(3,2,3) %(c)

for ii=1:length(S{1}.Ep)
plot(S{1}.tt,S{1}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(1,0,'$t=$3.1874','interpreter','latex','fontsize',FntSize)
line([S{1}.T,S{1}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))

fig_defaults(FntSize)

subplot(3,2,5) %(d)
clear legs
for ii=1:length(S{2}.Ep)
plot(S{2}.tt,S{2}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(7.5,0,'$t=$9.31','interpreter','latex','fontsize',FntSize)
line([S{2}.T,S{2}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))

fig_defaults(FntSize)




end

function subplots_of_UDD3(S,fignum)
%S is the input struct
clrs=COLORS;
FntSize=20;

figure(fignum)

subplot(3,2,1) %(a)

for ii=1:length(S{1}.Ep)
stem(S{1}.phi0(ii)/(pi/2),S{1}.Ep(ii),"filled",'linewidth',2,...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
end
xlabel('$\phi_0/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.98,1])

subplot(3,2,2) %(b)

for ii=1:length(S{2}.Ep)
stem(S{2}.phi0(ii)/(pi/2),S{2}.Ep(ii),"filled",'linewidth',2,...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
end

xlabel('$\phi_0/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.8,1])

subplot(3,2,3) %(c)

for ii=1:length(S{1}.Ep)
plot(S{1}.tt,S{1}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(2,0,'$t=$3.1877','interpreter','latex','fontsize',FntSize)
line([S{1}.T,S{1}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))
xlim([2,4])
fig_defaults(FntSize)

subplot(3,2,5) %(d)
clear legs
for ii=1:length(S{2}.Ep)
plot(S{2}.tt,S{2}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(14,-0.2,'$t=$15.9214','interpreter','latex','fontsize',FntSize)
line([S{2}.T,S{2}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))
xlim([13.5,18.5])
fig_defaults(FntSize)




end

function subplots_of_UDD4(S,fignum)
%S is the input struct
clrs=COLORS;
FntSize=20;

figure(fignum)

subplot(3,2,1) %(a)

for ii=1:length(S{1}.Ep)
stem(S{1}.phi0(ii)/(pi/2),S{1}.Ep(ii),"filled",'linewidth',2,'marker','o',...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
stem(S{1}.phi1(ii)/(pi/2),S{1}.Ep(ii),"filled",'linewidth',2,'marker','d',...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii},...
    'linestyle',':')
hold on
end
xlabel('$\phi_j/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.8,1])

subplot(3,2,2) %(b)

for ii=1:length(S{2}.Ep)
stem(S{2}.phi0(ii)/(pi/2),S{2}.Ep(ii),"filled",'linewidth',2,'marker','o',...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii})
hold on
stem(S{2}.phi1(ii)/(pi/2),S{2}.Ep(ii),"filled",'linewidth',2,'marker','d',...
    'color',clrs{ii},'markerfacecolor',clrs{ii},'markeredgecolor',clrs{ii},...
    'linestyle',':')
hold on
end

xlabel('$\phi_j/(\pi/2)$','interpreter','latex')
ylabel('$\epsilon_{p|q}/\epsilon_p^*$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.8,1])

subplot(3,2,3) %(c)

for ii=1:length(S{1}.Ep)
plot(S{1}.tt,S{1}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(3.23,0,'$t=$3.24','interpreter','latex','fontsize',FntSize)
line([S{1}.T,S{1}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))
xlim([3.22,3.26])
fig_defaults(FntSize)

subplot(3,2,5) %(d)
clear legs
for ii=1:length(S{2}.Ep)
plot(S{2}.tt,S{2}.n0n1(ii,:),'linewidth',2,...
    'color',clrs{ii})
legs{ii}=num2str(ii);
hold on
end
text(9.35,0.7,'$t\approx$9.35','interpreter','latex','fontsize',FntSize)
line([S{2}.T,S{2}.T],[-1,1],'color','k','linewidth',2,'linestyle','--')
xlabel('$t$ ($\mu$s)','interpreter','latex')
ylabel('$\textbf{n}_0\cdot\textbf{n}_1$','interpreter','latex')
legend(legs,'location','best','color','none','edgecolor','none','NumColumns',length(legs))
xlim([9.3,9.37])
fig_defaults(FntSize)




end

function [A,B]=Generate_distinct_HF(varargin) %tol_HF,Nnuc,max_HF

tol_HF = varargin{1};
Nnuc   = varargin{2};
if length(varargin)<3
   max_HF = 200;
else
   max_HF=varargin{3};
end
    
a = 10;  b = max_HF;  %Generate spins with HF parameters from 10 - max_HF kHz
A = a+(b-a).*rand(Nnuc,1);  %a+(b-a).*rand(N,1)
B = a+(b-a).*rand(Nnuc,1);
numA=length(A);


disp('Removing spins with similar HF params...')

cnt=0;
for ii=1:numA
    for jj=ii+1:numA

        if abs(A(ii)-A(jj))<tol_HF && abs(B(ii)-B(jj))<tol_HF
           
           cnt=cnt+1;
           indx_to_remove(cnt)=jj;
           
        end
    end
end

indx_to_remove = unique(indx_to_remove);
indx_to_keep   = setxor(1:length(A),indx_to_remove);

A=A(indx_to_keep); B=B(indx_to_keep);
disp(['# of distinct nuclei is:',num2str(length(A))])




end


function OUT=Maximize_Many_Tangles(Choose_Sequence,Nnuc,tol_HF,s0,s1,k,wL)

Pulse_Cutoff = 520;   maxNumber=30;
 
c=2*pi*1e-3;N=1;
[A,B]=Generate_distinct_HF(tol_HF,Nnuc);

W0 = @(A,B)  sqrt( (wL+s0*A)^2 +(s0*B)^2)*c;
W1 = @(A,B)  sqrt( (wL+s1*A)^2 +(s1*B)^2)*c;
W  = @(A,B)  W0(A,B)+W1(A,B);
  
time = @(A,B) 4*pi*(2*k-1)/(W(A,B));  
   
tt = linspace(time(A(1),B(1))-0.2,time(A(1),B(1))+0.2,1e4);

%Chose the 1st randomly generated spin and optimize the time a bit to hit 
%its resonance


switch Choose_Sequence
    
    case 'CPMG'
       
        for ii=1:length(tt)
            
            Spin = SubClass_U4Operations(wL,A(1),B(1),s0,s1,1,k,N);
            Spin = Spin.CPMG(tt(ii));
            Spin = Spin.Rot_Angles_And_Axes;
            n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
            
        end
        
        
    case 'UDD3'

        for ii=1:length(tt)
            
            Spin = SubClass_U4Operations(wL,A(1),B(1),s0,s1,1,k,N);
            Spin = Spin.UDD(3,tt(ii));
            Spin = Spin.Rot_Angles_And_Axes;
            n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
            
        end
    
    case 'UDD4'
        
        for ii=1:length(tt)
            
            Spin = SubClass_U4Operations(wL,A(1),B(1),s0,s1,1,k,N);
            Spin = Spin.UDD(4,tt(ii));
            Spin = Spin.Rot_Angles_And_Axes;
            n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
            
        end
        
end

[~,indx]=min(n0n1);
topt = tt(indx);


%Now, check the expected maxima of all spins:

switch Choose_Sequence
    
    case 'CPMG'
    
        for ii=1:length(A)

           Spin = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,N);
           Spin = Spin.CPMG(topt);
           Spin = Spin.Expected_Maxima(maxNumber);
           Nmax{ii} = Spin.Nmax;


        end

    case 'UDD3'

        for ii=1:length(A)

           Spin = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,N);
           Spin = Spin.UDD(3,topt);
           Spin = Spin.Expected_Maxima(maxNumber);
           Nmax{ii} = Spin.Nmax;


        end        
        
    case 'UDD4'
        Nall = 1:Pulse_Cutoff;
        parfor ii=1:length(A)
           
            for NN=1:Pulse_Cutoff
            
               Spin = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,NN);
               Spin = Spin.UDD(4,topt);
               Spin = Spin.Makhlin_Inv;
               ep(ii,NN) = Spin.Ep/(2/9);
              

            end
            
        end        
        
        for ii=1:length(A)
           
            Nmax{ii}=Nall(ep(ii,:)>0.98);
            
        end
        
        
end

for jj=1:length(Nmax)
    
    tempN = Nmax{jj} ;
    Nmax{jj}=tempN(tempN<=Pulse_Cutoff);
    
end

%Now, given Nmax, fix the first target set to be Nmax{1}

temp_Set = Nmax{1};

Iters=[];
cnt=0;
for jj=2:length(Nmax)
   
    
    Li = ismember(temp_Set,Nmax{jj});
    
    if any(Li)
        cnt=cnt+1;
        Iters = temp_Set(Li);
        store_spins(cnt)=jj; 
    end
    
    temp_Set=temp_Set(Li);
    
    if isempty(temp_Set) && ~isempty(Iters)
        temp_Set=Iters;
        
    elseif isempty(temp_Set) && isempty(Iters)
       
        error('Did not find common number of iterations for any spins.')
        
    end
    
    
    
    
end

Nopt=Iters;
Spins_For_Multipartite = [1,store_spins];
Acommon = A(Spins_For_Multipartite);
Bcommon = B(Spins_For_Multipartite);

disp(['The # of spins for multipartite gate is:',num2str(length(Acommon))])

%Finally, get nuclear one-tangles:

switch Choose_Sequence
    
    case 'CPMG'

        
        for ii=1:length(Acommon)

           Spin = SubClass_U4Operations(wL,Acommon(ii),Bcommon(ii),s0,s1,1,k,Nopt);
           Spin = Spin.CPMG(topt);
           Spin = Spin.Makhlin_Inv;
           Spin = Spin.Rot_Angles_And_Axes;
           n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
           phi0(ii) = Spin.angles{1};
           Ep(ii) =Spin.Ep/(2/9);

        end
         phi1=phi0;
    case 'UDD3'
        
        for ii=1:length(Acommon)

           Spin = SubClass_U4Operations(wL,Acommon(ii),Bcommon(ii),s0,s1,1,k,Nopt);
           Spin = Spin.UDD(3,topt);
           Spin = Spin.Makhlin_Inv;
           Spin = Spin.Rot_Angles_And_Axes;
           n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
           phi0(ii) = Spin.angles{1};
           Ep(ii) =Spin.Ep/(2/9);

        end
        phi1=phi0;
        
    case 'UDD4'

        for ii=1:length(Acommon)

           Spin = SubClass_U4Operations(wL,Acommon(ii),Bcommon(ii),s0,s1,1,k,Nopt);
           Spin = Spin.UDD(4,topt);
           Spin = Spin.Rot_Angles_And_Axes;
           n0n1(ii) = dot(Spin.axes{1},Spin.axes{2});
           phi0(ii) = Spin.angles{1};
           phi1(ii) = Spin.angles{2};
           Spin = Spin.Makhlin_Inv;
           Ep(ii) =Spin.Ep/(2/9);

        end

end

OUT.A = Acommon;
OUT.B = Bcommon;
OUT.N = Nopt;
OUT.T = topt;
OUT.Ep = Ep;
OUT.n0n1 = n0n1;
OUT.phi0 = phi0;
OUT.phi1 = phi1;


end


%======== Functions to reproduce results from paper ======================

function [A,B,T,N]=load_HF_from_paper(Sequence,Case_Option)

switch Sequence
    
    case 'CPMG'
        
        switch Case_Option
            
            case '1'
                
               A=[195.7761,27.7833,124.5269,100.0254,26.9260,65.7260,63.7672,106.8783,193.3990,144.1764];
               B=[49.6193,136.5072,128.3130,22.0722,181.3337,128.1540,74.9186,164.7410,122.0653,93.4147];
               T=3.1874; 
               N=56;
 
            case '2'
            
                A=[188.8471,56.3805,88.2939,56.5273,134.8014,82.9059,121.0999,10.2876];
                B=[131.3993,179.6982,109.6431,78.5107,150.6604,187.7117,73.4681,157.8233];
                T=9.3103; 
                N=8;
            
        end
        
            
    case 'UDD3'
        
       switch Case_Option
           
           case '1'
          
              A=[156.4447,140.3017,198.8740,66.0289,70.0823,123.2503,26.1352,41.6440,159.7784,61.7190,45.0811];
              B=[77.0341,86.0290,166.6354,49.3565,148.4476,121.9295,112.9580,103.9616,104.2800,76.8776,191.5879];
              T=3.1877; 
              N=487;
               
           case '2'
               
                A=[168.7800,82.9890,63.8160,136.9000,141.4600,142.1100,186.1000,199.6500];
                B=[12.8040,158.3000,88.1350,149.5200,99.4400,76.1910,56.7490,138.4300];
                T=15.9214; 
                N=93;
               
       end
        
       
    case 'UDD4'
    
    switch Case_Option
        
        case '1'
            

            A=[185.9721,66.7149,74.6909,142.1813,129.3018,176.7505,53.5991,22.8028,36.5410];
            B=[180.3179,101.6238,53.9081,92.3526,56.3928,56.9189,136.8595,92.3401,194.6969];
            T=3.2401; 
            N=252;

          
        case '2'

            A=[57.3009,83.4203,91.9723,167.8684,150.7603,81.2899,165.2526,179.0818];
            B=[157.2489,41.4065,183.3216,70.6491,190.5065,135.9573,99.1300,30.3376];
            T=9.3494; 
            N=41;
            
    end
    
end




end

function OUT=get_data_from_paper_CPMG(Case_Option,wL,s0,s1)

disp('=== Calculating parameters for multipartite gate for CPMG =========')

switch Case_Option
    
    case '1'
        
        [A,B,T,N]=load_HF_from_paper('CPMG','1');
   
    case '2'
        
        [A,B,T,N]=load_HF_from_paper('CPMG','2');        
        
end
k=1; %Doesnt matter because I am providing time.

parfor ii=1:length(A)

   Spin      = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,N) ;
   Spin      = Spin.CPMG(T);
   Spin      = Spin.Makhlin_Inv;
   Spin      = Spin.Rot_Angles_And_Axes;
   Ep(ii)    = Spin.Ep/(2/9);
   phi0(ii)  = Spin.angles{1};
   phi1(ii)  = Spin.angles{2};
   
end

tt   = linspace(T-3,T+3,1e3);
numA = length(A);

parfor ii=1:length(tt)
    for jj=1:numA
        Spin        = SubClass_U4Operations(wL,A(jj),B(jj),s0,s1,1,k,N) ;
        Spin        = Spin.CPMG(tt(ii));
        Spin        = Spin.Rot_Angles_And_Axes;
        n0n1(jj,ii) = dot(Spin.axes{1},Spin.axes{2});

    end
end


OUT.A    = A;
OUT.B    = B;
OUT.T    = T;
OUT.N    = N;
OUT.Ep   = Ep;
OUT.n0n1 = n0n1;
OUT.phi0 = phi0;
OUT.phi1 = phi1;
OUT.tt   = tt;

disp('=== Done =========')

end

function OUT=get_data_from_paper_UDD3(Case_Option,wL,s0,s1)

disp('=== Calculating parameters for multipartite gate for UDD3 =========')

switch Case_Option
    
    case '1'
        
        [A,B,T,N]=load_HF_from_paper('UDD3','1');
   
    case '2'
        
        [A,B,T,N]=load_HF_from_paper('UDD3','2');        
        
end
k=1; %Doesnt matter because I am providing time.

parfor ii=1:length(A)

   Spin      = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,N) ;
   Spin      = Spin.UDD(3,T);
   Spin      = Spin.Makhlin_Inv;
   Spin      = Spin.Rot_Angles_And_Axes;
   Ep(ii)    = Spin.Ep/(2/9);
   phi0(ii)  = Spin.angles{1};
   phi1(ii)  = Spin.angles{2};
   
end

tt   = linspace(T-3,T+3,1e3);
numA = length(A);

parfor ii=1:length(tt)
    for jj=1:numA
        Spin        = SubClass_U4Operations(wL,A(jj),B(jj),s0,s1,1,k,N) ;
        Spin        = Spin.UDD(3,tt(ii));
        Spin        = Spin.Rot_Angles_And_Axes;
        n0n1(jj,ii) = dot(Spin.axes{1},Spin.axes{2});

    end
end


OUT.A    = A;
OUT.B    = B;
OUT.T    = T;
OUT.N    = N;
OUT.Ep   = Ep;
OUT.n0n1 = n0n1;
OUT.phi0 = phi0;
OUT.phi1 = phi1;
OUT.tt   = tt;

disp('=== Done =========')

end

function OUT=get_data_from_paper_UDD4(Case_Option,wL,s0,s1)

disp('=== Calculating parameters for multipartite gate for UDD4 =========')

switch Case_Option
    
    case '1'
        
        [A,B,T,N]=load_HF_from_paper('UDD4','1');
   
    case '2'
        
        [A,B,T,N]=load_HF_from_paper('UDD4','2');        
        
end
k=1; %Doesnt matter because I am providing time.

parfor ii=1:length(A)

   Spin      = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,k,N) ;
   Spin      = Spin.UDD(4,T);
   Spin      = Spin.Makhlin_Inv;
   Spin      = Spin.Rot_Angles_And_Axes;
   Ep(ii)    = Spin.Ep/(2/9);
   phi0(ii)  = Spin.angles{1};
   phi1(ii)  = Spin.angles{2};
   
end

tt   = linspace(T-3,T+3,1e3);
numA = length(A);

parfor ii=1:length(tt)
    for jj=1:numA
        Spin        = SubClass_U4Operations(wL,A(jj),B(jj),s0,s1,1,k,N) ;
        Spin        = Spin.UDD(4,tt(ii));
        Spin        = Spin.Rot_Angles_And_Axes;
        n0n1(jj,ii) = dot(Spin.axes{1},Spin.axes{2});

    end
end


OUT.A    = A;
OUT.B    = B;
OUT.T    = T;
OUT.N    = N;
OUT.Ep   = Ep;
OUT.n0n1 = n0n1;
OUT.phi0 = phi0;
OUT.phi1 = phi1;
OUT.tt   = tt;

disp('=== Done =========')


end

function Out=Get_Gate_Errors(wL,s0,s1,N,T,Ntarget,Sequence)

tol_HF = 3; max_HF = 300;  Nnuc  = 1e5;  k=1;

[Aunw,Bunw]=Generate_distinct_HF(tol_HF,Nnuc,max_HF);

switch Sequence
    
    case 'CPMG'

        parfor ii=1:length(Aunw)

           Spin      = SubClass_U4Operations(wL,Aunw(ii),Bunw(ii),s0,s1,1,k,N);
           Spin      = Spin.CPMG(T);
           Spin      = Spin.Makhlin_Inv;
           Spin      = Spin.Rot_Angles_And_Axes;
           phi0(ii)  = Spin.angles{1};
           phi1(ii)  = Spin.angles{2};
           n0(ii,:)  = Spin.axes{1};
           n1(ii,:)  = Spin.axes{2};
           Ep(ii)    = Spin.Ep/(2/9);

        end        
        
    case 'UDD3'
        
        parfor ii=1:length(Aunw)

           Spin      = SubClass_U4Operations(wL,Aunw(ii),Bunw(ii),s0,s1,1,k,N);
           Spin      = Spin.UDD(3,T);
           Spin      = Spin.Makhlin_Inv;
           Spin      = Spin.Rot_Angles_And_Axes;
           phi0(ii)  = Spin.angles{1};
           phi1(ii)  = Spin.angles{2};
           n0(ii,:)  = Spin.axes{1};
           n1(ii,:)  = Spin.axes{2};
           Ep(ii)    = Spin.Ep/(2/9);

        end        
        
    case 'UDD4'
        
        parfor ii=1:length(Aunw)

           Spin      = SubClass_U4Operations(wL,Aunw(ii),Bunw(ii),s0,s1,1,k,N);
           Spin      = Spin.UDD(4,T);
           Spin      = Spin.Makhlin_Inv;
           Spin      = Spin.Rot_Angles_And_Axes;
           phi0(ii)  = Spin.angles{1};
           phi1(ii)  = Spin.angles{2};
           n0(ii,:)  = Spin.axes{1};
           n1(ii,:)  = Spin.axes{2};
           Ep(ii)    = Spin.Ep/(2/9);

        end        
        
end

Out=Interval_Gate_Error(Ntarget,Aunw,Bunw,n0,n1,phi0,phi1,Ep);


end

function intervals=Get_one_tangles_intervals


intervals={[0       ,  7  ]*1e-6,...   %1
           [0.70001 , 1.5 ]*1e-5,...   %2
           [1.50001 ,  4  ]*1e-5,...   %3
           [4.00001 ,  7  ]*1e-5,...   %4
           [7.00001 ,  9  ]*1e-5,...   %5
           [0.90001 ,  1.1]*1e-4,...   %6
           [1.10001 ,  4  ]*1e-4,...   %7
           [4.10001 ,  6  ]*1e-4,...   %8
           [6.10001 ,  8  ]*1e-4,...   %9
           [0.9 ,  1.1]*1e-3,...       %10
           [2 ,  4    ]*1e-3,...       %11
           [4 ,  6    ]*1e-3,...       %12
           [6 ,  8    ]*1e-3,...       %13
           [0.9 , 1.1 ]*1e-2,...       %14
           [2 , 4     ]*1e-2,...       %15
           [4 , 6     ]*1e-2,...       %16
           [6 , 8     ]*1e-2,...       %17
           [0.9 , 1.1 ]*1e-1,...       %18
           [1.4 , 1.6 ]*1e-1,...       %19
           [1.9 , 2.1 ]*1e-1,...       %20
           [2.4 , 2.6 ]*1e-1,...       %21
           [2.9 , 3.1 ]*1e-1,...       %22
           [3.4 , 3.6 ]*1e-1,...       %23
           [3.9 , 4.1 ]*1e-1,...       %24
           [4.4 , 4.6 ]*1e-1,...       %25
           [4.9 , 5.1 ]*1e-1,...       %26
           [5.4 , 5.6 ]*1e-1,...       %27
           [5.9 , 6.1 ]*1e-1,...       %28
           [6.4 , 6.6 ]*1e-1,...       %29
           [6.9 , 7.1 ]*1e-1,...       %30
           [7.4 , 7.6 ]*1e-1,...       %31
           };




end

function Out=Interval_Gate_Error(Ntarget,A,B,n0,n1,phi0,phi1,Ep)

intervals   = Get_one_tangles_intervals;
Spin_CutOff = 6;

[~,col1]=size(n0); [~,col2]=size(n1);

if col1~=3 || col2~=3
    error('n0 and n1 need to be Nnuc x 3 arrays.')
end

disp('============== Entering calculation of Gate errors =================')

for jj=1:length(intervals)
    
   indices       = find( (Ep>intervals{jj}(1)) & (Ep<intervals{jj}(2))   );
   L             = length(indices);
   
   if isempty(indices)
       
       disp('Did not find any unwanted spins for this interval.')
       Atemp=nan;
       Btemp=nan;
       Eptemp=nan;
       Infid=nan;
       
   else
       
       if L>Spin_CutOff
           
           indices=indices(1:Spin_CutOff);
           L=Spin_CutOff;
       end
       
       Atemp    = A(indices); 
       Btemp    = B(indices);
       n0temp   = n0(indices,:); 
       n1temp   = n1(indices,:);
       phi0temp = phi0(indices); 
       phi1temp = phi1(indices);
       Eptemp   = Ep(indices);
       
       for ll=1:Spin_CutOff %for each interval, get the gate error by gradually increasing the size
                  %of the unwanted spin bath.
          if ll<=L        
          TEMP      = SubClass_Ent_and_Fid;
          TEMP      = TEMP.Gate_Infid(ll+Ntarget,Ntarget,...
                      phi0temp(1:ll),phi1temp(1:ll),n0temp(1:ll,:),n1temp(1:ll,:));
          Infid(ll) = TEMP.Infid;
          else
              Infid(ll)=nan;
          end
       end
       
       
       
   end
   
    Out.A{jj}     = Atemp;
    Out.B{jj}     = Btemp;
    Out.Ep{jj}    = Eptemp;
    Out.Infid{jj} = Infid;
    
    
end

disp('========= Done. ===============================================')

%Get also the center of each interval and start/end for error bars of plots

for jj=1:length(intervals)
    
    Out.center(jj)       = 1/2*(intervals{jj}(1)+intervals{jj}(2));
    Out.ErrBar_start{jj} = abs(Out.center(jj) -intervals{jj}(1));
    Out.ErrBar_end{jj}   = abs(Out.center(jj) -intervals{jj}(2));
    
end
   


end

function Out=Get_Data_Trivial_Evol(kappa,t,N,s0,s1,wL,Ntarget)

MaxHF=300;

A = 2:1:300;
c=2*pi*1e-3;


%find B's:
cnt=0;
for ii=1:length(A)
        
       if  (8*kappa*pi/t)^2 > (wL*c+s0*A(ii)*c)^2 && 1/abs(s0)*sqrt(   (8*kappa*pi/t)^2 - (wL*c+s0*A(ii)*c)^2 )/c <= MaxHF
          
         cnt=cnt+1;
         Bsol(cnt)=  1/abs(s0)*sqrt(   (8*kappa*pi/(t))^2 - (wL*c+s0*A(ii)*c)^2 );  %B is in MHz
         BsolkHz(cnt)=  Bsol(cnt)/c;
         AsolkHz(cnt) = A(ii);
         
         
       end
       
    
end



cnt=0;
for ii=1:length(A)
        
       if  (8*kappa*pi/t)^2 > (wL*c+s1*A(ii)*c)^2 && 1/abs(s1)*sqrt(   (8*kappa*pi/t)^2 - (wL*c+s1*A(ii)*c)^2 )/c <= MaxHF
          
         cnt=cnt+1;
         Bsol2(cnt)=  1/abs(s1)*sqrt(   (8*kappa*pi/(t))^2 - (wL*c+s1*A(ii)*c)^2 );  %B is in MHz
         BsolkHz2(cnt)=  Bsol2(cnt)/c;
         AsolkHz2(cnt) = A(ii);
         
         
       end
       
    
end



if exist('AsolkHz2') && exist('AsolkHz')
A = [AsolKHz,AsolkHz2];
B = [BsolKHz,BsolkHz2];
elseif exist('AsolkHz2')
A = AsolkHz2;
B = BsolkHz2; 
elseif exist('AsolkHz')
A = AsolkHz;
B = BsolkHz;   
else
    return
end

A=A(1:6); %just get 6 nuclei that evolve trivially
B=B(1:6);

Nnuc=1;
k=1;

for ii=1:length(A)
    
    Spin = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,Nnuc,k,N);
    Spin = Spin.CPMG(t);
    Spin = Spin.Rot_Angles_And_Axes;
    Spin = Spin.Makhlin_Inv;
    
    phi0(ii)=Spin.angles{1};
    phi1(ii)=Spin.angles{2};
    n0(ii,:)=Spin.axes{1};
    n1(ii,:)=Spin.axes{2};
    Ep(ii)  =Spin.Ep/(2/9);

end

for ll=1:length(A)
    
    
    TEMP = SubClass_Ent_and_Fid;
    TEMP = TEMP.Gate_Infid(ll+Ntarget,Ntarget,...
                      phi0(1:ll),phi1(1:ll),n0(1:ll,:),n1(1:ll,:));
    Infid(ll) = TEMP.Infid;
    
end

Out.A=A;
Out.B=B;
Out.Ep=Ep;
Out.Infid=Infid;


end



