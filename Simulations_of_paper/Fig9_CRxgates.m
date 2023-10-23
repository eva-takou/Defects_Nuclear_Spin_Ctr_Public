function Fig9_CRxgates
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%Script to find the optimal CRx gates to compare with multi-spin scheme

CaseNum             = 23;
ResNum              = 3;
max_Res             = 12;
max_Time            = [1500,69,1500];
Unwanted_Tangle_Tol = [0.2,0.14,0.31]; %For C4: 0.2, For C5:

for ii=1:3

S{ii}=Optimal_Control(CaseNum,ResNum,ii,max_Res,max_Time(ii),Unwanted_Tangle_Tol(ii));
end


%Plot the rotation components

for ii=1:length(S)
Nuclear_Indx(ii) = S{ii}.NuclearIndex;
phi0(ii)         = S{ii}.phi0T;
phi1(ii)         = S{ii}.phi1T;
n0(ii,:)         = S{ii}.n0T;
n1(ii,:)         = S{ii}.n1T;
end

close all;
FntSize=22;
xlab = {num2str(Nuclear_Indx(1)),num2str(Nuclear_Indx(2)),num2str(Nuclear_Indx(3))};

figure(1) %Sequential



subplot(2,2,1)


stem(1:3,[n0(:,1),n1(:,1)],"filled",'linewidth',2,'marker','d')
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('Rot. axes')
legend({'$n_{x,0}$','$n_{x,1}$'},'interpreter','latex',...
    'location','best','color','none','edgecolor','none','NumColumns',1)
fig_defaults(FntSize)


subplot(2,2,2)
stem(1:3,[n0(:,3),n1(:,3)],"filled",'linewidth',2,'marker','>')
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('Rot. axes')
legend({'$n_{z,0}$','$n_{z,1}$'},'interpreter','latex',...
    'location','best','color','none','edgecolor','none','NumColumns',1)

fig_defaults(FntSize)

subplot(2,2,3)
stem(1:3,phi0/(pi/2),"filled",'linewidth',2)
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('$\phi_0/(\pi/2)$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.9,1.1])

figure(2)

%MultiSpin
[phi0,~,n0,n1]=Get_Rot_Components_MultiSpin(CaseNum,ResNum);



subplot(2,2,1)


stem(1:3,[n0(:,1),n1(:,1)],"filled",'linewidth',2,'marker','d')
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('Rot. axes')
legend({'$n_{x,0}$','$n_{x,1}$'},'interpreter','latex',...
    'location','best','color','none','edgecolor','none','NumColumns',1)
fig_defaults(FntSize)


subplot(2,2,2)
stem(1:3,[n0(:,3),n1(:,3)],"filled",'linewidth',2,'marker','>')
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('Rot. axes')
legend({'$n_{z,0}$','$n_{z,1}$'},'interpreter','latex',...
    'location','best','color','none','edgecolor','none','NumColumns',1)

fig_defaults(FntSize)

subplot(2,2,3)
stem(1:3,phi0/(pi/2),"filled",'linewidth',2)
xlabel('C nuclear index')
set(gca,'xtick',[1:3],'xticklabels',xlab)
ylabel('$\phi_0/(\pi/2)$','interpreter','latex')
fig_defaults(FntSize)
ylim([0.9,2])




end

function [phi0,phi1,n0,n1]=Get_Rot_Components_MultiSpin(CaseNum,ResNum)

B0 = 403; %in Gauss
B0 = B0*1e-4;
gamma_C13 = 6.728284 * 1e7; %rad/T*s
wL        = gamma_C13 * B0 ; 
wL        = wL*1e-6; %2.7 MHz
wL        = wL*1e3; %in KHz
wL        = wL/(2*pi);  %431.5484 kHz
s0        = 0;
s1        = -1;
k         = 1;
Nnuc      = 1;
load('27SpinsData.mat')

At   = Out{CaseNum,ResNum}.Atarget;
Bt   = Out{CaseNum,ResNum}.Btarget;
N    = Out{CaseNum,ResNum}.N;
t    = Out{CaseNum,ResNum}.t;

for ii=1:length(At)
    
    Spin       = SubClass_U4Operations(wL,At(ii),Bt(ii),s0,s1,Nnuc,k,N);
    Spin       = Spin.CPMG(t);
    Spin       = Spin.Rot_Angles_And_Axes;
    phi0(ii)   = Spin.angles{1};
    phi1(ii)   = Spin.angles{2};
    n0(ii,:)   = Spin.axes{1};
    n1(ii,:)   = Spin.axes{2};
    
end






end


function [A,B]=get_HF


A=[-20.72,...
   -23.22,...
   -31.25,...
   -14.07,...
   -11.346,...
   -48.58,...
   -8.32,...
   -9.79,...
   213.154,...
    17.643,...
    14.548,...
    20.569,...
    8.029,...
    -19.815,...
    -13.961,...
    -4.66,...
    -5.62,...
    -36.308,...
    24.399,...
    2.690,...
    1.212,...
    7.683,...
    -3.177,...
    -4.225,...
    -3.873,...
    -3.618,...
    -4.039];
B=[12,...
    13,...
    8,...
    13,...
    59.21,...
    9,...
    3,...
    5,...
    3,...
    8.6,...
    10,...
    41.51,...
    21,...
    5.3,...
    9,...
    7,...
    5,...
    26.62,...
    24.81,...
    11,...
    13,...
    4,...
    2,...
    0,...
    0,...
    0,...
    0];

A=A(1:23);
B=B(1:23);


end



function S=Optimal_Control(CaseNum,ResNum,Spin_Indx,max_Res,max_Time,Unwanted_Tangle_Tol)

%=================================================

load('27SpinsData.mat')

[A,B]=get_HF;

Target_Nuc = Out{CaseNum,ResNum}.TargetNuc(Spin_Indx);

At = A(Target_Nuc);
Bt = B(Target_Nuc);

Au = A(setxor(1:length(A),Target_Nuc));
Bu = B(setxor(1:length(A),Target_Nuc));



disp(['Finding optimal cases for C:',num2str(Target_Nuc)])

S=Optimal_CRx(At,Bt,max_Time,max_Res,Au,Bu,Unwanted_Tangle_Tol);
    
S.NuclearIndex = Target_Nuc;


end




function Out=Optimal_CRx(At,Bt,max_Time,max_Res,Au,Bu,Unwanted_Tangle_Tol)
%Find the optimal CRx of a single spin with minimal cross-talk.

B0 = 403; %in Gauss
B0 = B0*1e-4;
gamma_C13 = 6.728284 * 1e7; %rad/T*s
wL        = gamma_C13 * B0 ; 
wL        = wL*1e-6; %2.7 MHz
wL        = wL*1e3; %in KHz
wL        = wL/(2*pi);  %431.5484 kHz
wL        = 432;
s0        = 0;
s1        = -1;

c         = 2*pi*1e-3;
W0        = sqrt( (wL+s0*At)^2+(s0*Bt)^2 )*c;
W1        = sqrt( (wL+s1*At)^2+(s1*Bt)^2 )*c;
W         = W0+W1;
tk        = @(k) 4*pi*(2*k-1)/W; %resonance time

%For each resonance obtain the expected maxima

Nnuc = 1;
N    = 1; 
cnt  = 0;
maxN = 60;
for k=1:max_Res
     
    
    Spin     = SubClass_U4Operations(wL,At,Bt,s0,s1,Nnuc,k,N);
    Spin     = Spin.CPMG(tk(k));
    Spin     = Spin.Expected_Maxima(maxN);
    tempN    = Spin.Nmax;
    tempN    = tempN(tempN*tk(k)<max_Time); %satisfy time constraints

    %Arrange everything into cases:
    for ll=1:length(tempN)
       cnt=cnt+1; 
        
       Iters(cnt)=tempN(ll);
       times(cnt)=tk(k);
       Res(cnt)  =k;
       
    end
    
    
end

Max_Case = cnt;

%For each one of this cases find unwanted one tangles:

Ep = zeros(Max_Case,length(Au));
numA = length(Au);

disp(['Running evolution of ',num2str(numA),' unwanted spins... '])
parfor ii=1:Max_Case
    
   for jj=1:numA
       
      Spin      = SubClass_U4Operations(wL,Au(jj),Bu(jj),s0,s1,Nnuc,Res(ii),Iters(ii));
      Spin      = Spin.CPMG(times(ii)) ;
      Spin      = Spin.Makhlin_Inv;
      Ep(ii,jj) = Spin.Ep/(2/9);
       
   end
    
    
end
disp('Done.')
%Find when the unwanted one-tangles are minimal:

for ii=1:Max_Case
    
   temp =  Ep(ii,:);
   
   if any(temp>Unwanted_Tangle_Tol) %Not acceptable case
      
       Iters(ii)=nan;
       times(ii)=nan;
       Res(ii)  =nan;
      
   end
    
end
Ep    = Ep(~isnan(Iters),:);
Iters = Iters(~isnan(Iters));
times = times(~isnan(times));
Res   = Res(~isnan(Res));

if isempty(Iters)
    
    error(['Could not find cases that satisfy constraints. ' ...
    'Increase unwanted one-tangle tolerance or resonance # search.'])
elseif length(Iters)==1
    
    %Find also the target one tangle, and evolution of target spin.
    
    k=1; Nnuc=1;
    
    Spin  = SubClass_U4Operations(wL,At,Bt,s0,s1,Nnuc,k,Iters);
    Spin  = Spin.CPMG(times);
    Spin  = Spin.Makhlin_Inv;
    Spin  = Spin.Rot_Angles_And_Axes;
    
    EpT   = Spin.Ep/(2/9); 
    phi0T = Spin.angles{1};
    phi1T = Spin.angles{2};
    n0T   = Spin.axes{1};
    n1T   = Spin.axes{2};
    
    %Get the Gate error:
    
    for ii=1:length(Au)
        
       Spin      = SubClass_U4Operations(wL,Au(ii),Bu(ii),s0,s1,Nnuc,k,Iters);
       Spin      = Spin.CPMG(times);
       Spin      = Spin.Rot_Angles_And_Axes;
       phi0U(ii) = Spin.angles{1};
       phi1U(ii) = Spin.angles{2};
       n0U(ii,:) = Spin.axes{1};
       n1U(ii,:) = Spin.axes{2};
        
    end
    disp('Calculating gate error...')
    TEMP = SubClass_Ent_and_Fid;
    TEMP = TEMP.Gate_Infid(length(Au)+1,1,phi0U,phi1U,n0U,n1U);
    disp('Done.')
    
    Out.t       = times;
    Out.N       = Iters;
    Out.k       = Res;
    Out.Total_T = times.*Iters; 
    Out.EpU     = Ep;
    Out.EpT     = EpT;
    Out.Infid   = TEMP.Infid;
    Out.phi0T   = phi0T;
    Out.phi1T   = phi1T;
    Out.n0T     = n0T;
    Out.n1T     = n1T;
end





end

