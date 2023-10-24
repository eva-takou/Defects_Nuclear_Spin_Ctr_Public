%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------
%%
clc;
clear;
close all;




a=10;
b=250;

Nnuc=7;
NT = 2;


A=a+(b-a)*rand(1,Nnuc);
B=a+(b-a)*rand(1,Nnuc);

Target_Spins=[];
while true
   
    indx = randi([2,Nnuc]);
    if ~any(Target_Spins==indx)
        
        Target_Spins=[Target_Spins,indx];
        
    end
    
    if length(Target_Spins)==NT
        break
    end
    
    
end


Unwanted_Spins = setxor(1:Nnuc,Target_Spins);

NT = length(Target_Spins);


At = A(Target_Spins);
Bt = B(Target_Spins);
wL=314;  s0=0;   s1=-1;  N=30;  t=5.4; k=1;
%% Get the target evolution

temp = SuperClass_Sequences(wL,At,Bt,s0,s1,NT,k,N);
temp = temp.CPMG(t);
U0   = temp.Uevol;

%% Get the total evolution

temp = SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,k,N);
temp = temp.CPMG(t);
U    = temp.Uevol;
%% Get the infidelity via Kraus
temp2 = SubClass_Ent_and_Fid;
Ek    = temp2.Get_Kraus(U,Nnuc,Target_Spins);


Infid_Brute = temp2.Get_Infid_From_Kraus(U0,Ek,Target_Spins)
%% Get the infidelity based on closed form expression

Aunw = A(Unwanted_Spins);
Bunw = B(Unwanted_Spins);

for jj=1:length(Aunw)
    
   Spin = SubClass_U4Operations(wL,Aunw(jj),Bunw(jj),s0,s1,1,1,N);
   Spin = Spin.CPMG(t);
   Spin = Spin.Rot_Angles_And_Axes;
   phi0(jj) = Spin.angles{1};
   phi1(jj) = Spin.angles{2};

   n0(jj,:) = Spin.axes{1};
   n1(jj,:) = Spin.axes{2};

end


test=SubClass_Ent_and_Fid;
Infid_Alt=test.Gate_Infid(Nnuc,NT,phi0,phi1,n0,n1).Infid




