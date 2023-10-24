function OUT=Fidelity_reduction_subroutine(tolHF,nn,Sequence_Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to calculate how the fidelity reduces due to correlations with
%unwanted nuclei. In each unwanted one-tangle interval we keep up to 6
%unwanted spins.
%
%Input:  tolHF: how different to be the parameters of the random generation
%          %nn: how many random spins to generate in total. 
%Output: A struct with the CPMG error, and other relevant parameters.
%--------------------------------------------------------------------------
%
%
%----------- Set some fixed parameters ------------------------------------
    
wL=314;     
s0=1/2;     
s1=-1/2;      

%Define target parameters
switch Sequence_Option
    
    case 'CPMG'

        %k=2;
        Nopt=8;    
        topt=9.31;
        Afinal=[188.85, 56.381, 88.294,56.527, 134.8, 82.906, 121.1, 10.288];
        %Bfinal=[131.4, 179.7, 109.64, 78.51,150.66, 187.71, 73.468,157.82];
        
    case 'UDD3'
        %k=3; 
        Nopt=93;    
        topt=15.9214;
        Afinal=[168.78, 82.989, 63.816,136.9, 141.46, 142.11, 186.1, 199.65];
        %Bfinal=[12.804, 158.3, 88.135, 149.52,99.44, 76.191, 56.749,138.43];        
        
    case 'UDD4'
        
        %k=2;
        Nopt=41;    
        topt=9.35;
        Afinal=[57.301, 83.42, 91.972,167.87, 150.76, 81.29, 165.25, 179.08];
        %Bfinal=[157.25, 41.407, 183.32, 70.649,190.51, 135.96, 99.13,30.338];  
        
end



Ntarget=length(Afinal);
%% Verify that the evolution operator of the target spins is correct:
%X   = [0 1 ; 1 0];   
%Y   = [0 -1i ; 1i 0];    
%Z   = [1 0 ; 0 -1];
% s00 = [1 0 ; 0 0];     
% s11 = [0 0 ; 0 1];
% 
% Rn=@(phi0,n) expm(-1i*phi0/2*( n(1)*X+n(2)*Y+n(3)*Z   ));
% 
% parfor ii=1:length(Afinal)
%     
%    TargetSpin = SubClass_U4Operations(wL,Afinal(ii),Bfinal(ii),s0,s1,1,k,Nopt);
%    TargetSpin = TargetSpin.CPMG(topt);
%    TargetSpin = TargetSpin.Rot_Angles_And_Axes;
%    phi0_Nopt(ii) = TargetSpin.angles{1};
%    phi1_Nopt(ii) = TargetSpin.angles{2};
%    n0(ii,:) = TargetSpin.axes{1};
%    n1(ii,:) = TargetSpin.axes{2};
%    
% end
% 
% test = SuperClass_Sequences(wL,Afinal,Bfinal,s0,s1,length(Afinal),k,Nopt);
% test = test.CPMG(topt);     
% 
% V0=1; V1=1;
% 
% for ii=1:length(Afinal)
%     
%     
%     V0 = kron(V0,Rn(phi0_Nopt(ii),n0(ii,:)));
%     V1 = kron(V1,Rn(phi1_Nopt(ii),n1(ii,:)));
% end
% 
% U0= -kron(s00,V0)-kron(s11,V1);
% 
% if abs(norm(U0-test.Uevol))<1e-5
% 
% else
%     error('Need to put some minus signs on some rot. angles to get the right U.')
% end

%% Now generate the other random spins, and find their one-tangles and at the end do a plot

min_HF = 10;
max_HF = 250;

[EpU,n0U,n1U,phi0U,phi1U,Arand,Brand]=unwanted_one_tangles_random_spins(wL,s0,s1,Nopt,topt,min_HF,max_HF,tolHF,nn,Sequence_Option);


%%  Find how tracing out some number of spins with particular Ep changes the Infid

intervals=get_unwanted_one_tangle_bins;

[Astore,Bstore,Epstore,Infid_U]=...
    arrange_one_tangles_in_intervals(Arand,Brand,EpU,phi0U,phi1U,n0U,n1U,Ntarget);

for kk=1:length(intervals)
    
    Trace1spin(kk) =Infid_U{kk}(1);
    Trace2spin(kk) =Infid_U{kk}(2);
    Trace3spin(kk) =Infid_U{kk}(3);
    Trace4spin(kk) =Infid_U{kk}(4);
    Trace5spin(kk) =Infid_U{kk}(5);
    Trace6spin(kk) =Infid_U{kk}(6);
    x_c(kk) =  (intervals{kk}(2)+intervals{kk}(1))/2;  %center of intervals
    
end

yData={Trace1spin,Trace2spin,Trace3spin,Trace4spin,Trace5spin,Trace6spin};
   
%Error-bars
start=zeros(1,length(intervals));
final=zeros(1,length(intervals));

for ii=1:length(intervals)
    start(ii)=abs(-intervals{ii}(1)+x_c(ii));
    final(ii)=abs(-intervals{ii}(2)+x_c(ii));
end

OUT.Infid1 = Trace1spin;
OUT.Infid2 = Trace2spin;
OUT.Infid3 = Trace3spin;
OUT.Infid4 = Trace4spin;
OUT.Infid5 = Trace5spin;
OUT.Infid6 = Trace6spin;
OUT.xData = x_c;
OUT.yData = yData;
OUT.ErrBarStart  = start;
OUT.ErrBarEnd  = final;

OUT.Interval_HFA = Astore;
OUT.Interval_HFB = Bstore;
OUT.Interval_EP  = Epstore;






end






