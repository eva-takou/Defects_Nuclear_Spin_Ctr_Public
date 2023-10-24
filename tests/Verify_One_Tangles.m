%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------
%%
clc;
clear;
close all;

A=[60,50,100];
B=[30,20,10];
wL=314;
s0=0;
s1=-1;
k=1;
N=10;

t=4;

Nnuc=length(A);


temp=SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,k,N);
temp=temp.CPMG(t);

U=temp.Uevol;

ent_temp = SubClass_Ent_and_Fid;
ent_temp.Uval=U;

ent_temp=ent_temp.one_tangles_General_U_byparts;

one_tangle_el  = ent_temp.one_tangles{1}(1);  %Full numerical evaluation
one_tangle_nuc = ent_temp.one_tangles{1}(2:end);
%%


for jj=1:length(A)
    
   Spin=SubClass_U4Operations(wL,A(jj),B(jj),s0,s1,1,1,N);
   Spin=Spin.CPMG(t);
   Spin=Spin.Makhlin_Inv;
   
   G1(jj)=real(Spin.Makhlin{1}{1});
   ep_nuc(jj) = 2/9*(1-abs(G1(jj)));
   
end

ep_el_analy=1/3-1/3^(length(A)+1)*prod(1+2*G1);

if abs(ep_el_analy-one_tangle_el)>1e-9
    
    error('Analytical does not agree with numerical (one-tangle of electron)')
    
end

if any(abs(one_tangle_nuc-ep_nuc)>1e-9)
    
    error('Analytical does not agree with numerical (one-tangle of nuclei)')
    
end



