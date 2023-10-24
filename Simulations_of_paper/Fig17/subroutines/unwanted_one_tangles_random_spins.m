function [EpU,n0U,n1U,phi0U,phi1U,Arand,Brand]=unwanted_one_tangles_random_spins(wL,s0,s1,N,t,min_HF,max_HF,tolHF,nn,Sequence_Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to calculate the unwanted one-tangles for the CPMG sequence.
%
%Input: wL:Larmor frequency
%       s0/s1: Electron's spin projections
%       N: number of CMPG iterations
%       t: time of CPMG unit
%       min_HF/max_HF: min/max range of HF parameters in kHz
%       nn: # of random spins to generate
%       Sequence_Option: 'CPMG' or 'UDD3' or 'UDD4'
%Output: EpU: one-tangles of unwanted spins
%        n0U: Rot axis of unwanted spins when electron is in |0>
%        n1U: Rot axis of unwanted spins when electron is in |1>
%        phi0U: Rot angel of unwanted spins when electron is in |0>
%        phi1U: Rot angel of unwanted spins when electron is in |1>
    
Arand = min_HF+(max_HF-min_HF)*rand(1,nn);
Brand = min_HF+(max_HF-min_HF)*rand(1,nn);

indx=[];

%---- Remove HF parameters that are close up to a tolerance ---------------
for ii=1:length(Arand)
    
    for jj=ii+1:length(Arand)
        
        if abs(Arand(ii)-Arand(jj))<tolHF && abs(Brand(ii)-Brand(jj))<tolHF
            
            indx=[indx,jj];
            break
            
        end
        
    end
end

Arand(indx)=[];
Brand(indx)=[];

k=1;
Nnuc=1;

switch Sequence_Option
    
    case 'CPMG'
        
        parfor ii=1:length(Arand)

                OtherSpins = SubClass_U4Operations(wL,Arand(ii),Brand(ii),s0,s1,Nnuc,k,N);
                OtherSpins = OtherSpins.CPMG(t);
                OtherSpins = OtherSpins.Rot_Angles_And_Axes;

                n0U(ii,:) = OtherSpins.axes{1};
                n1U(ii,:) = OtherSpins.axes{2};
                phi0U(ii) = OtherSpins.angles{1};
                phi1U(ii) = OtherSpins.angles{2};

                OtherSpins = OtherSpins.Makhlin_Inv;
                EpU(ii)    = OtherSpins.Ep/(2/9);

        end
        
    case 'UDD3'
        
        parfor ii=1:length(Arand)

                OtherSpins = SubClass_U4Operations(wL,Arand(ii),Brand(ii),s0,s1,Nnuc,k,N);
                OtherSpins = OtherSpins.UDD(3,t);
                OtherSpins = OtherSpins.Rot_Angles_And_Axes;

                n0U(ii,:) = OtherSpins.axes{1};
                n1U(ii,:) = OtherSpins.axes{2};
                phi0U(ii) = OtherSpins.angles{1};
                phi1U(ii) = OtherSpins.angles{2};

                OtherSpins = OtherSpins.Makhlin_Inv;
                EpU(ii)    = OtherSpins.Ep/(2/9);

        end
        
    case 'UDD4'

        
        parfor ii=1:length(Arand)

                OtherSpins = SubClass_U4Operations(wL,Arand(ii),Brand(ii),s0,s1,Nnuc,k,N);
                OtherSpins = OtherSpins.UDD(4,t);
                OtherSpins = OtherSpins.Rot_Angles_And_Axes;

                n0U(ii,:) = OtherSpins.axes{1};
                n1U(ii,:) = OtherSpins.axes{2};
                phi0U(ii) = OtherSpins.angles{1};
                phi1U(ii) = OtherSpins.angles{2};

                OtherSpins = OtherSpins.Makhlin_Inv;
                EpU(ii)    = OtherSpins.Ep/(2/9);

        end
        
end

%Sort in ascending order
[EpU,Ind] = sort(EpU);
Arand     = Arand(Ind);
Brand     = Brand(Ind);
phi0U     = phi0U(Ind);
phi1U     = phi1U(Ind);
n0U       = n0U(Ind,:);
n1U       = n1U(Ind,:);



end