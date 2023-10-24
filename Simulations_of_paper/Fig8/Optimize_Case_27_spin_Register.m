function S=Optimize_Case_27_spin_Register(Spin_Indx,k)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------

Target_Tangle_Tol   = 0.8;
Unwanted_Tangle_Tol = 0.14;
Tmax                = 1500;

s0 = 0;  
s1 = -1;  

[A,B,wL]=get_HF; 

Nnuc = length(A); 

T=Resonance_Time(A(Spin_Indx),B(Spin_Indx),wL,s0,s1,k,'Primal');

tt = linspace(T-0.25,T+0.25,50);
Nmax = round(Tmax/(T+0.25));

G1 = @(phi0,phi1,n0n1) ( cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2)  )^2;
ep = @(phi0,phi1,n0n1) 1-G1(phi0,phi1,n0n1);


phi0 = zeros(length(tt),Nnuc);    phi1 = zeros(length(tt),Nnuc);
n0   = zeros(length(tt),Nnuc,3);   n1  = zeros(length(tt),Nnuc,3);
N=1;

%--------- Get evolution for 1 unit ---------------------------------------
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

%---- For each time and each spin, get the entangling power as we vary N---
phiN=@(phi,N) acos(cos(N*phi));

cnt=0; 
All_Nuc = 1:Nnuc;

for jj=1:length(tt)
    
    for N=1:Nmax
        
        Ep=zeros(1,Nnuc);
        
        for kk=1:Nnuc
        
            phi0N  = phiN(phi0(jj,kk),N);
            phi1N  = phiN(phi1(jj,kk),N);
            n0n1   = dot(n0(jj,kk,:),n1(jj,kk,:));
            Ep(kk) = ep(phi0N,phi1N,n0n1); %One-tangle of k-th nuclear spin
         
        end
        
        indx_Target = Ep>Target_Tangle_Tol;
        indx_Rest   = Ep<Target_Tangle_Tol;
        tempTarget  = Ep(indx_Target);
        tempRest    = Ep(indx_Rest);
      
        if length(tempTarget)==1 || isempty(tempTarget) %Only 1 spin with maximal one-tangle (ignore)
            
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
                Out{cnt}.EpT              = tempTarget;
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

%------------ Post-select based on maximum mean(epT) ----------------------

indices_Max_Nuclei=find(Number_Of_Target_Nuc==max(Number_Of_Target_Nuc)); %Max possible # of nuclei entangled with electron


EpMean_T = EpMean_T(indices_Max_Nuclei);
[~,indx] = max(EpMean_T,[],'all','linear'); %Post-select maximum entanglement
indx     = indices_Max_Nuclei(indx);
    
S=getfield(Out,{indx});
S=S{1};

k=1;    
for kk=1:length(S.UnwantedNuc)-4 %Last 4 ignored because B=0, and they do not affect gate error.
       
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
