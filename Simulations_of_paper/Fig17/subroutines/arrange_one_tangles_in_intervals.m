function [Astore,Bstore,Epstore,Infid_U]=arrange_one_tangles_in_intervals(A,B,Ep,phi0,phi1,n0,n1,Ntarget)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

intervals=get_unwanted_one_tangle_bins;


for jj=1:length(intervals)
       
            
    indices= find( (Ep>intervals{jj}(1)) & (Ep<intervals{jj}(2))   );

    
    if length(indices)>6
          
        indices=indices(1:6);
        
    end
            
    A_U   = A(indices); 
    B_U   = B(indices);   
    ph0_U = phi0(indices);
    ph1_U = phi1(indices);
    n0_U  = n0(indices,:);  
    n1_U  = n1(indices,:);
    Ep_U  = Ep(indices);

    Astore{jj}   = A_U;   
    Bstore{jj}   = B_U;    
    Epstore{jj}  = Ep_U;
%     ph0store{jj} = ph0_U;  
%     ph1store{jj} = ph1_U;  
%     n0store{jj}  = n0_U;   
%     n1store{jj}  = n1_U;
% 
    for ii=1:6       % give bound up to 5-6 spins?

        if ii<=length(A_U)

            TEMP          = SubClass_Ent_and_Fid;   
            TEMP          = TEMP.Gate_Infid(Ntarget+ii,Ntarget,ph0_U(1:ii),ph1_U(1:ii),n0_U(1:ii,:),n1_U(1:ii,:)); 
            infid(ii,jj)  = TEMP.Infid; 

        else

            infid(ii,jj)  = nan; 


        end

    end

    Infid_U{jj}=infid(:,jj);
            
end            




end