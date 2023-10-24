function [DeltaA,Ep_Opt]=subroutine_Fig13(wL,s0,s1,At,Bt,B)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to reproduce Fig. 13 of Ref. Phys. Rev X 13, 011004 (2023)
%



Time_Max=2000; %2000 \mus
Res_Max=12;

N=[];
t=[];

for k=1:Res_Max
    
   T      = Resonance_Time(At,Bt,wL,s0,s1,k,'Primal');
   
   Target = SubClass_U4Operations(wL,At,Bt,s0,s1,1,k,1);
   Target = Target.CPMG(T);
   Target = Target.Expected_Maxima(30);
   Iters  = Target.Nmax;
   Iters  = Iters(Iters*T<Time_Max);
   
   N = [Iters,N];
   t = [repmat(T,[1,length(Iters)]),t];
   
end

%--------------------------------------------------------------------------
%Make sure that for those iters and times we have maximal target one tangle

parfor ii=1:length(N)
    
   Target    = SubClass_U4Operations(wL,At,Bt,s0,s1,1,k,N(ii));
   Target    = Target.CPMG(t(ii));
   Target    = Target.Makhlin_Inv;
   Ep_Target = Target.Ep/(2/9);
    
   if Ep_Target<0.9
      error('The estimated maxima of Ep are wrong.') 
   end
   
end

%--------- Generate unwanted spins ----------------------------------------

A    = 10:0.25:200;
numA = length(A);

parfor ii=1:length(t)
    
   for jj=1:numA
       
      ThisSpin =  SubClass_U4Operations(wL,A(jj),B,s0,s1,1,k,N(ii));
      ThisSpin =  ThisSpin.CPMG(t(ii));
      ThisSpin =  ThisSpin.Makhlin_Inv;
      Ep(ii,jj) = ThisSpin.Ep/(2/9);
      
   end
    
end

t_Opt=zeros(1,numA);
N_Opt=zeros(1,numA);
Ep_Opt=zeros(1,numA);

for jj=1:numA
    
    temp = Ep(:,jj);
    
    [~,indx]=min(temp);
    
    t_Opt(jj)  = t(indx);
    N_Opt(jj)  = N(indx);
    Ep_Opt(jj) = min(temp);
    
end


[DeltaA,order]=sort(abs(At-A));
Ep_Opt = Ep_Opt(order);





end