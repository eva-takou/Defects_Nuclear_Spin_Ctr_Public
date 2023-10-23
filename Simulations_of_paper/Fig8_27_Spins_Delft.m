function Out=Fig8_27_Spins_Delft(Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

clc;

switch Option
    
    case 'Load'
        
        load('27SpinsData.mat')
        

    case 'Run_Again'
        
        Out=Get_Optimal_Cases;
end



%Collect the data into arrays:
cnt=0;
for k=1:5
    
    for Spin_Indx=1:27
        cnt=cnt+1;
        temp=Out{Spin_Indx,k};
        if isstruct(temp)
            
            EpT_mean(cnt) = mean(temp.EpT);
            EpU_mean(cnt) = mean(temp.EpUnw);
            Target_Spins(cnt) = length(temp.TargetNuc);
            Gate_Time(cnt)    = temp.Gate_Time;
            N(cnt)            = temp.N;
            Infid(cnt)        = temp.Infid;
            
        else
            EpT_mean(cnt)=nan;
            EpU_mean(cnt)=nan;
            Target_Spins(cnt)=nan;
            Gate_Time(cnt)=nan;
            N(cnt)=nan;
            Infid(cnt)=nan;
        end
        
    end
end

yData={EpT_mean,EpU_mean,Target_Spins,N,Gate_Time,Infid};

close all;
rows  = 3;
cols  = 2;
Markers  = {'o','s','v','>','<'};
FntSize  = 22;

ylims={[0.85,1],[0,0.04],[2,7],[0,400],[0,1500],[0,0.2]};
ylabels={'$\bar{\epsilon}_p^T$','$\bar{\epsilon}_p^U$',...
    '# of nuclei','# of pulses','Gate time (\mus)','$1-F$'};
interp = {'latex','latex','none','none','tex','latex'};


figure(1)

for ii=1:length(yData)
   
    subplot(rows,cols,ii)
    
    for k=1:5
        cnst = (k-1)*27;
        stem(1+cnst:27+cnst,yData{ii}(1+cnst:27+cnst),"filled",...
            'marker',Markers{k},'linewidth',2,'markersize',8)
        hold on
        
    
    end
    
        fig_defaults(FntSize)
        xlabel('Case #')
        ylabel(ylabels{ii},'interpreter',interp{ii})
        ylim(ylims{ii})
        xlim([1,135])
        if ii==1
            
            legend({'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$'},...
                'location','best','interpreter','latex','color','none',...
                'edgecolor','none','NumColumns',5)
            
        end
       
    
end






end

function [A,B,wL]=get_HF


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

B0 = 403; %in Gauss
B0 = B0*1e-4;

gamma_C13 = 6.728284 * 1e7; %rad/T*s
wL        = gamma_C13 * B0 ; 
wL        = wL*1e-6; %2.7 MHz
wL        = wL*1e3; %in KHz
wL        = wL/(2*pi);  %431.5484 kHz

end


function Out=Get_Optimal_Cases

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,1}=choose_starting_case(Spin_Indx,1);
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,2}=choose_starting_case(Spin_Indx,2);
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,3}=choose_starting_case(Spin_Indx,3);
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,4}=choose_starting_case(Spin_Indx,4);
   
end
toc

tic
parfor Spin_Indx=1:27
   
     Out{Spin_Indx,5}=choose_starting_case(Spin_Indx,5);
   
end
toc



end


function S=choose_starting_case(Spin_Indx,k)

Target_Tangle_Tol   = 0.8;
Unwanted_Tangle_Tol = 0.14;

s0 = 0;  s1 = -1;  c = 2*pi*1e-3; Tmax=1500;

[A,B,wL]=get_HF; 

Nnuc = length(A); 

W0 = @(A,B)  sqrt( (wL+s0*A)^2 +(s0*B)^2)*c;
W1 = @(A,B)  sqrt( (wL+s1*A)^2 +(s1*B)^2)*c;
W  = @(A,B) W0(A,B)+W1(A,B);

tk = @(A,B) 4*pi*(2*k-1)/(W(A,B));
T  = tk(A(Spin_Indx),B(Spin_Indx));

tt = linspace(T-0.25,T+0.25,50);
Nmax = round(Tmax/(T+0.25));

G1 = @(phi0,phi1,n0n1) ( cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2)  )^2;
ep = @(phi0,phi1,n0n1) 1-G1(phi0,phi1,n0n1);


phi0 = zeros(length(tt),Nnuc);    phi1 = zeros(length(tt),Nnuc);
n0   = zeros(length(tt),Nnuc,3);   n1  = zeros(length(tt),Nnuc,3);
N=1;
for jj=1:length(tt)
    
    for kk=1:Nnuc
    
        Spin         = SubClass_U4Operations(wL,A(kk),B(kk),s0,s1,1,k,N) ;
        Spin         = Spin.CPMG(tt(jj));
        Spin         = Spin.Rot_Angles_And_Axes;
        phi0(jj,kk)  = Spin.angles{1};
        phi1(jj,kk)  = Spin.angles{2};
        n0(jj,kk,:)  = Spin.axes{1};
        n1(jj,kk,:)  = Spin.axes{2};
        
    end
end

%For each time and each spin, get the entangling power as we vary N
phiN=@(phi,N) acos(cos(N*phi));

cnt=0; All_Nuc = 1:Nnuc;

for jj=1:length(tt)
    
    for N=1:Nmax
        
        Ep=zeros(1,Nnuc);
        
        for kk=1:Nnuc
        
        phi0N = phiN(phi0(jj,kk),N);
        phi1N = phiN(phi1(jj,kk),N);
        n0n1  = dot(n0(jj,kk,:),n1(jj,kk,:));
        Ep(kk) = ep(phi0N,phi1N,n0n1);
         
        end
        
        indx_Target = Ep>Target_Tangle_Tol;
        indx_Rest   = Ep<Target_Tangle_Tol;
        tempTarget  = Ep(indx_Target);
        tempRest    = Ep(indx_Rest);
      
        if length(tempTarget)==1 || isempty(tempTarget)
            continue
            
        else %found >=2 target spins.
            
            if all(tempRest<Unwanted_Tangle_Tol) %Acceptable Case
                
                cnt=cnt+1;
                Out{cnt}.TargetNuc        = All_Nuc(indx_Target);
                Out{cnt}.UnwantedNuc      = All_Nuc(indx_Rest);
                Out{cnt}.Atarget          = A(indx_Target);
                Out{cnt}.Btarget          = B(indx_Target);
                Out{cnt}.Aunw             = A(indx_Rest);
                Out{cnt}.Bunw             = B(indx_Rest);
                Out{cnt} .EpT             = tempTarget;
                Out{cnt}.EpUnw            = tempRest;
                Out{cnt}.t                = tt(jj);
                Out{cnt}.N                = N;
                Out{cnt}.Gate_Time        = N*tt(jj);               
                Number_Of_Target_Nuc(cnt) = length(All_Nuc(indx_Target));
                EpMean_T(cnt)             = mean(tempTarget);
           end
            
        end
        
    end
    
end

if cnt==0 %found no case
    S=nan;
    return
    
end

indices_Max_Nuclei=find(Number_Of_Target_Nuc==max(Number_Of_Target_Nuc));
    %Pick the case where we entangle the maximum possible number of nuclei with the
    %electron, while keeping mean(epT) the highest

EpMean_T=EpMean_T(indices_Max_Nuclei);
[~,indx] = max(EpMean_T,[],'all','linear');

indx = indices_Max_Nuclei(indx);
    
S=getfield(Out,{indx});
S=S{1};

k=1;    
for kk=1:length(S.UnwantedNuc)-4
       
  Spin        = SubClass_U4Operations(wL,S.Aunw(kk),S.Bunw(kk),s0,s1,1,k,S.N);
  Spin        = Spin.CPMG(S.t);
  Spin        = Spin.Rot_Angles_And_Axes;
  PH0(kk)     = Spin.angles{1};
  PH1(kk)     = Spin.angles{2};
  N0(kk,:)    = Spin.axes{1};
  N1(kk,:)    = Spin.axes{2};
      
end
   

%Get the gate error:
disp('==== Entering calculation of gate error =========================')

Nnuc=length(S.UnwantedNuc)-4+length(S.Atarget);

TEMP    = SubClass_Ent_and_Fid;
TEMP    = TEMP.Gate_Infid(Nnuc,length(S.Atarget),PH0,PH1,N0,N1);
S.Infid = TEMP.Infid;

disp('=========== Done ===============')

        



end




