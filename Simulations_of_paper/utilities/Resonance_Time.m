function time=Resonance_Time(A,B,wL,s0,s1,k,Option)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

c=2*pi*1e-3;

W0 = sqrt((wL+s0*A)^2+(s0*B)^2)*c;
W1 = sqrt((wL+s1*A)^2+(s1*B)^2)*c;
W  = W0+W1;


switch Option
    
    
    case 'Primal'
        
        time =  4*pi*(2*k-1)/W;
        
        
    case 'Secondary'
        
        time =  8*pi*(2*k-1)/W;
        
end




end