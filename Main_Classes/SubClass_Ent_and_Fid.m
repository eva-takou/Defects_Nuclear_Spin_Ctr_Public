classdef SubClass_Ent_and_Fid
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

  properties
      
      Uval     {mustBeNumeric}  
      Ek            %Kraus operators (could be empty)
      Fid           %Fidelity (as calculated from the total system using the Kraus operators)
      maxBoundPi    %Max bound of one-tangles for CR evolutions
      one_tangles   %One-tangles for a general U (or CR-type)
      Max_Bound     %Max bound of one-tangles for general U.
      epM_uni       %M-way entangling power (of target subspace)
      epM_nonuni    %M-way entangling power (of target subspace)
      epM_dephased  %M-way entangling power (of target subspace)
  end
  
    properties (Constant, Hidden)
       
        Q=1/sqrt(2)*[1  0  0  1i ; 
             0 1i  1  0  ;
             0 1i -1  0  ;
             1  0  0  -1i];
          
        
    end
  
    
  methods (Static)  
         
        function [Out]=GKF(U,U0,Nnuc,Target_Spins,option)
        %Output: Kraus operators corresponding to environmental spins.
        %        Gate fidelity due to cross-talk i.e., deviation of U
        %        from U0 due to presence of unwanted nuclei.
        %Input:  U: Total gate U. 
        %        U0: Target gate (of target subspace)
        %        Target_Spins: Ordered indices of target spins 
        %        Option: "Only_Kraus" or "And_Fid" to output only 
        %        Kraus operators or gate fidelity as well.
        
            %Give two target gates U0 and two m indices.
            
            %U0{1} is the target for single spin + electron.
            %U0{2} is the target for the "target spins"+ electron.
            
            %First make sure that arguments were passed correctly.
            
            if isempty(U0{2}) && ~isempty(Target_Spins)
                error('You provided target spins but not a target gate for the calculation of the fidelity.')
            end
            
            if ~isempty(U0{2}) && isempty(Target_Spins)
                error('You provided target gate for many spins, but not the indices of the target spins.')
            end
            
            if ~strcmp(option,"And_Fid") && ~strcmp(option,"Only_Kraus")
                error('Option can be either "And_Fid" or "Only_Kraus".')
            end
            
            %Part 1: get Kraus for the evolution of electron + single
            %spins.
            
            if isempty(Target_Spins)
                   [ww,ss]=size(U0{1});
                if ss~=4 || ww~=4
                  error('The input U0{1} is not 4x4.')   
                end
                
            end
            
            Nenvir = Nnuc-1;
            dim_E  = 2^Nenvir;
            id_envir = eye(dim_E);  
            ket00 = zeros(dim_E,1); %column vector
            ket00(1) = 1;   %assume the state |0000000000.....>
            e0= kron(eye(4),ket00); 
            
            %P2 switches the position between 2 nuclear spins.
            %It is a permutation matrix.
            %Repeatedly we switch the position of the i-th spin
            %To find the Kraus of the i-th spin with the electron.
            
            
            P2 = [1  0  0  0 ;...
                  0  0  1  0 ;...
                  0  1  0  0 ;...
                  0  0  0  1 ];

                %Kraus of 1st spin.
                Ek1 = zeros(4,4,dim_E);
                
                for ii=1:dim_E
                    
                    ek = kron(eye(4),id_envir(ii,:));
                    Ek1(:,:,ii) = ek*U*e0;
                    
                end            
                
                Identities = repmat({eye(2)},[1 1 Nnuc-1]);
                Identities{1}=P2;                
                
            %Kraus of 2nd spin till last.
            
            Ek_ALL = zeros(4,4,Nnuc,dim_E);
            
            for Ni=2:Nnuc

                if Ni==2
                    temp=Identities;
                    P = superkron(temp);
                    
                else

                    for ii=1:Ni-1

                        temp  = Identities;
                        temp = circshift(temp,(Ni-1)-1); %moves P2 by some positions to the right
                        P=superkron(temp);

                    if isequal(temp{1},P2)

                        continue %go to the next for loop

                    else
                        for jj=1:Ni-2
                        
                            temp = circshift(temp,-1);
                            P = superkron(temp)*P;
                            
                        end

                    end

                    end
                end

                P = kron(eye(2),P);

                    for ll=1:dim_E
                    ek = kron(eye(4),id_envir(ll,:));

                    Ek_ALL(:,:,Ni,ll) = ek*P*U*P'*e0; %Note this starts from the 2nd position
                    end

            end          
            
            for ll=1:dim_E
            
                Ek_ALL(:,:,1,ll)= Ek1(:,:,ll);
                
            end
            
            %The above are the Kraus when we consider single spin + electron.
            
            %Whether to output Kraus only or Infidelity as well.
            switch option
                
                case "Only_Kraus"
                    Out.Kraus=Ek_ALL;
                    
                case "And_Fid"
                    
                   Infid=zeros(1,Nnuc);
                   m=4;
                    for jj=1:Nnuc
                        
                        expr1=0;   expr2=0;  
                        
                        for ii=1:dim_E
                        
                        Mk=U0{1}'*Ek_ALL(:,:,jj,ii);  
                        

                        expr1= expr1 + Mk'*Mk;
                        expr2= expr2 + abs(trace(Mk))^2;

                        end
                        
                        Fid=1/(m*(m+1))* (   trace(expr1) + expr2     );
                        Infid(jj) = 1-Fid;
                    end           
           
                    Out.Kraus=Ek_ALL;
                    Out.Infid = Infid;
                    
            end
            
            
            %Now, we might also want to find Kraus and Fid, for some
            %particular spins. (Might want to keep more than single spin +
            %electron).
            
            %Give number of subsystems we want to find Kraus and Fid.
            
            
            if ~isempty(Target_Spins)
            
            %Need to rearrange unwanted spins at the last position
            %Use setxor, to find which spins are the unwanted.
            
            
            Ntarget = numel(Target_Spins);
            dim_TE = 2^(Ntarget+1); %dimension of target+electron
            All_Spins = 1 : Nnuc;
            
            Rest_Spins = setxor(Target_Spins,All_Spins); %This gives us the position of unwanted spins
            Rest_Spins = sort(Rest_Spins,'descend');  %Rest = Unwanted
            
            
            Nrest = Nnuc - Ntarget;
            dim_E = 2^Nrest;
            id_envir = eye(dim_E);
            
            ket00 = zeros(dim_E,1); %column vector
            ket00(1) = 1;   %assume the state |0000000000.....>
            e0= kron(eye(dim_TE),ket00); 
             
            %To push the i-th spin to the end
            %Apply P2 to the ith position
            %The rest will be identities
            %Keep on applying P2 till it is pushed all the
            %way to the last element of Identities
            
            %We choose first the spins that are already close
            %to the end
            
            %Special case: if all the target spins are in the 1st positions
            %there is no need for re-ordering and then projecting.
            
            
            if Ntarget<2 
                error('Need to have at least 2 target spins.')
            end
            
           
            
            
            Identities = repmat({eye(2)},[1 1 Nnuc-1]);
            
            %Could happen that I give Rest_Spins=end.
            %So if this is the case skip everything
            
            Unew = U;
            
            for ii = 1 : length(Rest_Spins)
                
                temp = Identities;
                
                if Rest_Spins(ii)==All_Spins(end)
                    break
                    
                end
                temp{Rest_Spins(ii)}=P2; %Put P2 to the i-th position.
                P=superkron(temp);
                
                
                for jj = 1 : Nnuc - Rest_Spins(ii)
                    
                    if isequal(temp{end},P2)
                        
                        break %it breaks the inner loop (stops j) & 
                        %continues running the same i-th loop until we go to next i
                    else
                    temp = circshift(temp,1); %move P2 to the right
                    P = superkron(temp)*P;
                    end
                    
                end
                
                %do a check here:
                if ~isequal(temp{end},P2)
                    error('The last entry of the list temp is not P2!')
                end
                
                P = kron(eye(2),P);
                Unew=P*Unew*P'; %I think it is now fine. Then, we just need to "trace out" the unwanted spins
                
                
                
            end
            
            
            
            %Now trace out the unwanted spins:
            
            Ek_Target_Spins = zeros(dim_TE,dim_TE,dim_E);
             for ll=1:dim_E
                    ek = kron(eye(dim_TE),id_envir(ll,:));
                    Ek_Target_Spins(:,:,ll) = ek*Unew*e0; 
             end
            
            
            %Finally, get the fidelity as well:
            
                   
                   m=dim_TE;
                   [ww,ss]=size(U0{2});
                   
                   if ww~=dim_TE || ss~=dim_TE
                       error('U0{2} is not 2^(Ntarget+1)x2^(Ntarget+1)')
                   end
                        expr1=0;   expr2=0;   
                        for ll=1:dim_E
                        
                        Mk=U0{2}'*Ek_Target_Spins(:,:,ll);  
                        

                        expr1= expr1 + Mk'*Mk;
                        expr2= expr2 + abs(trace(Mk))^2;

                        
                        
                        end
                        
                        Fid=1/(m*(m+1))* (   trace(expr1) + expr2     );
                        Infid_MultiQ = 1-Fid;
                              
                        Out.Kraus_MultiQ = Ek_Target_Spins;
                        Out.Infid_MultiQ = Infid_MultiQ;
             
            
            
            else
                return
                
            end
            
            
            
            
            
            
                
        end

        function [Out]=GKFnew(U,U0,Nnuc,Target_Spins,option)
            
            %Give two target gates U0 and two m indices.
            
            %U0{1} is the target for single spin + electron.
            %U0{2} is the target for the "target spins"+ electron.
            
            %First make sure that arguments were passed correctly.
            
            if isempty(U0{2}) && ~isempty(Target_Spins)
                error('You provided target spins but not a target gate for the calculation of the fidelity.')
            end
            
            if ~isempty(U0{2}) && isempty(Target_Spins)
                error('You provided target gate for many spins, but not the indices of the target spins.')
            end
            
            if ~strcmp(option,"And_Fid") && ~strcmp(option,"Only_Kraus")
                error('Option can be either "And_Fid" or "Only_Kraus".')
            end
            
            %Part 1: get Kraus for the evolution of electron + single
            %spins.
            
            if isempty(Target_Spins)
                   [ww,ss]=size(U0{1});
                if ss~=4 || ww~=4
                  error('The input U0{1} is not 4x4.')   
                end
                
            end
            
            
            Nenvir = Nnuc-1;
            dim_E  = 2^Nenvir;
            id_envir = eye(dim_E);  
            ket00 = zeros(dim_E,1); %column vector
            ket00(1) = 1;   %assume the state |0000000000.....>
            e0= kron(eye(4),ket00); 
            
            %P2 switches the position between 2 nuclear spins.
            %It is a permutation matrix.
            %Repeatedly we switch the position of the i-th spin
            %To find the Kraus of the i-th spin with the electron.
            
            
            P2 = [1  0  0  0 ;...
                  0  0  1  0 ;...
                  0  1  0  0 ;...
                  0  0  0  1 ];

                %Kraus of 1st spin.
                
                Ek1 = zeros(4,4,dim_E);
                
                for ii=1:dim_E
                    ek = kron(eye(4),id_envir(ii,:));
                    Ek1(:,:,ii) = ek*U*e0;
                end            
                
                Identities = repmat({eye(2)},[1 1 Nnuc-1]);
                Identities{1}=P2;                
                
            %Kraus of 2nd spin till last.
            
            Ek_ALL = zeros(4,4,Nnuc,dim_E);
            
            for Ni=2:Nnuc

                if Ni==2
                    temp=Identities;
                    P = superkron(temp);
                    
                else

                    for ii=1:Ni-1

                        temp  = Identities;
                        temp = circshift(temp,(Ni-1)-1); %moves P2 by some positions to the right
                        P=superkron(temp);

                    if isequal(temp{1},P2)

                        continue %go to the next for loop

                    else
                        for jj=1:Ni-2
                        temp = circshift(temp,-1);
                        P = superkron(temp)*P;
                        end

                    end


                    end
                end

                P = kron(eye(2),P);

                    for ll=1:dim_E
                    ek = kron(eye(4),id_envir(ll,:));

                    Ek_ALL(:,:,Ni,ll) = ek*P*U*P'*e0; %Note this starts from the 2nd position
                    end


            end          
            
            
            for ll=1:dim_E
            Ek_ALL(:,:,1,ll)= Ek1(:,:,ll);
            end
            
            %The above are the Kraus when we consider single spin + electron.
            
            
            %Whether to output Kraus only or Infidelity as well.
            switch option
                
                case "Only_Kraus"
                    Out.Kraus=Ek_ALL;
                    
                case "And_Fid"
                    
                    
                   
                   Infid=zeros(1,Nnuc);
                   m=4;
                    for jj=1:Nnuc
                        
                        expr1=0;   expr2=0;  
                        
                        for ii=1:dim_E
                        
                        Mk=U0{1}'*Ek_ALL(:,:,jj,ii);  
                        

                        expr1= expr1 + Mk'*Mk;
                        expr2= expr2 + abs(trace(Mk))^2;

                        
                        
                        end
                        
                        Fid=1/(m*(m+1))* (   trace(expr1) + expr2     );
                        Infid(jj) = 1-Fid;
                    end           
           
                    
          
                    
                    Out.Kraus=Ek_ALL;
                    Out.Infid = Infid;
                    
                    
                    
            end
            
            
            %Now, we might also want to find Kraus and Fid, for some
            %particular spins. (Might want to keep more than single spin +
            %electron).
            
            %Give number of subsystems we want to find Kraus and Fid.
            
            
            if ~isempty(Target_Spins)
            
            %Need to rearrange unwanted spins at the last position
            %Use setxor, to find which spins are the unwanted.
            
            
            Ntarget = numel(Target_Spins);
            dim_TE = 2^(Ntarget+1); %dimension of target+electron
            All_Spins = 1 : Nnuc;
            
            Rest_Spins = setxor(Target_Spins,All_Spins); %This gives us the position of unwanted spins
            Rest_Spins = sort(Rest_Spins,'descend');  %Rest = Unwanted
            
            
            Nrest = Nnuc - Ntarget;
            dim_E = 2^Nrest;
            id_envir = eye(dim_E);
            
            ket00 = zeros(dim_E,1); %column vector
            ket00(1) = 1;   %assume the state |0000000000.....>
            e0= kron(eye(dim_TE),ket00); 
             
            %To push the i-th spin to the end
            %Apply P2 to the ith position
            %The rest will be identities
            %Keep on applying P2 till it is pushed all the
            %way to the last element of Identities
            
            %We choose first the spins that are already close
            %to the end
            
            %Special case: if all the target spins are in the 1st positions
            %there is no need for re-ordering and then projecting.
            
            
            if Ntarget<2 
                error('Need to have at least 2 target spins.')
            end
            
           
            
            
            Identities = repmat({eye(2)},[1 1 Nnuc-1]);
            
            %Could happen that I give Rest_Spins=end.
            %So if this is the case skip everything
            
            Unew = U;
            
            for ii = 1 : length(Rest_Spins)
                
                temp = Identities;
                
                if Rest_Spins(ii)==All_Spins(end)
                    break
                    
                end
                temp{Rest_Spins(ii)}=P2; %Put P2 to the i-th position.
                P=superkron(temp);
                
                
                for jj = 1 : Nnuc - Rest_Spins(ii)
                    
                    if isequal(temp{end},P2)
                        
                        break %it breaks the inner loop (stops j) & 
                        %continues running the same i-th loop until we go to next i
                    else
                    temp = circshift(temp,1); %move P2 to the right
                    P = superkron(temp)*P;
                    end
                    
                end
                
                %do a check here:
                if ~isequal(temp{end},P2)
                    error('The last entry of the list temp is not P2!')
                end
                
                P = kron(eye(2),P);
                Unew=P*Unew*P'; %I think it is now fine. Then, we just need to "trace out" the unwanted spins
                
                
                
            end
            
            
            
            %Now trace out the unwanted spins:
            
            Ek_Target_Spins = zeros(dim_TE,dim_TE,dim_E);
             for ll=1:dim_E
                    ek = kron(eye(dim_TE),id_envir(ll,:));
                    Ek_Target_Spins(:,:,ll) = ek*Unew*e0; 
             end
            
            
            %Finally, get the fidelity as well:
            
                   
                   m=dim_TE;
                   [ww,ss]=size(U0{2});
                   
                   if ww~=dim_TE || ss~=dim_TE
                       error('U0{2} is not 2^(Ntarget+1)x2^(Ntarget+1)')
                   end
                        expr1=0;   expr2=0;   
                        for ll=1:dim_E
                        
                        Mk=U0{2}'*Ek_Target_Spins(:,:,ll);  
                        

                        expr1= expr1 + Mk'*Mk;
                        expr2= expr2 + abs(trace(Mk))^2;

                        
                        
                        end
                        
                        Fid=1/(m*(m+1))* (   trace(expr1) + expr2     );
                        Infid_MultiQ = 1-Fid;
                              
                        Out.Kraus_MultiQ = Ek_Target_Spins;
                        Out.Infid_MultiQ = Infid_MultiQ;
             
            
            
            else
                return
                
            end
            
            
            
            
            
            
                
        end
        
        function [Out]=Gate_Infid(Nnuc,Ntarget,phi0,phi1,n0,n1)
        %Input: Nnuc: total # of nuclear spins
        %       Ntarget: # of nuclear spins of target subspace
        %       phi0,phi1: Rot angles of unwanted nuclei
        %       n0,n1: Rot axes of unwanted nuclei
        %Output: Gate infidelity.
        
        %Analytical evaluation of gate infidelity which expresses deviation
        %of target gate (of target subsoace) from ideal due to the presence 
        %of unwanted nuclei.
        
         dTarget = 2^(Ntarget+1);   
         Nunw    = Nnuc-Ntarget;
         
         if Nunw~=length(phi0)
             
             error('Inconsistent # of unwanted spins and rot angles argument.')
             
         end
         
        %NumberOfKet ranges from 0 till 2^{L-K}-1 where L#of all nuclei
        %                                               K#of target ones
        
        %loop through all these combinations, read the m-th bit
        %and if the m-th position is 0, we are supposed to get
        %the factor cos(phi_0/2)-1i*nz0*sin(phi_0/2) for when the electron
        %is in 0 (and similarly for 1).
        %if the m-th position is 1 then we get the other coefficient
        
        expr=0;
        
        for NumberOfKet = 0:2^(Nunw)-1 
        
            CurrentKet=dec2bin(NumberOfKet,Nunw)-'0';
            
            %To get the logical indices do this:
            
            Nuc_in0 = CurrentKet==0;
            Nuc_in1 = CurrentKet==1;
            
            ang0_el0    = phi0(Nuc_in0); ang0_el1    = phi1(Nuc_in0);
            ang1_el0    = phi0(Nuc_in1); ang1_el1    = phi1(Nuc_in1);

            %for nuclei in |0>
            f0_ThisKraus_el0 = prod(cos(ang0_el0/2)-1i*(n0(Nuc_in0,3).').*sin(ang0_el0/2));
            f0_ThisKraus_el1 = prod(cos(ang0_el1/2)-1i*(n1(Nuc_in0,3).').*sin(ang0_el1/2));
            
            %for nuclei in |1>
            f1_ThisKraus_el0 = prod((-1i)*sin(ang1_el0/2).*(n0(Nuc_in1,1)+1i*n0(Nuc_in1,2)).');
            f1_ThisKraus_el1 = prod((-1i)*sin(ang1_el1/2).*(n1(Nuc_in1,1)+1i*n1(Nuc_in1,2)).');
             
            
            if  isempty(ang0_el0) && ~isempty(ang1_el0) %all nuclei in |1>
                 
                f0_ThisKraus_el0=1;
                f0_ThisKraus_el1=1;
                
            elseif ~isempty(ang0_el0) && isempty(ang1_el0) %all nuclei in |0>
                
                f1_ThisKraus_el0=1;
                f1_ThisKraus_el1=1;
                
            end
            
             expr = expr + abs(f0_ThisKraus_el0*f1_ThisKraus_el0+...
                               f0_ThisKraus_el1*f1_ThisKraus_el1)^2;
             
        end
        
       
        F          = 1/(dTarget+1)*(1+2^(Ntarget-1)*expr); 
        Out.Infid  = 1-F;   
       
        end
        
        
  end
  
  methods (Static)
      
      function obj=one_tangles_General_U_Max_Bound(D)
      %To calculate the one-tangles of general U, based on the partition
      %method of the paper on multipartite entanglement.
      %This bound is not tight! (Only when AME(2n,d) exists).
          dQ = 2;
          
          c  = 1;
          ab = D-1;
        
          expr_bound=0;
        
          for ii=1:D
           
              multip  = nchoosek(D,ii);
              x_prime = ii; 
              y_prime = D-ii;
            
              dim_abx_prime = 2^(ab+x_prime);
              dim_cy_prime  = 2^(c+y_prime);
              min_dim       = min(dim_abx_prime,dim_cy_prime);
              expr_bound    = expr_bound + multip*1/min_dim;
            
          end
          
          x_prime = 0; 
          y_prime = D;
        
          dim_abx_prime = 2^(ab+x_prime);
          dim_cy_prime  = 2^(c+y_prime);
        
          min_dim       = min(dim_abx_prime,dim_cy_prime);
          expr_bound    = expr_bound + 1/min_dim;
          expr_bound    = dQ^D/(dQ+1)^D *expr_bound;
          expr_bound    = 1-expr_bound;
          obj.Max_Bound = expr_bound;
           
      end
      
  end
  
  methods %One-tangles

      
      function obj=one_tangles_CR(obj,Nnuc,phi0,phi1,n0,n1)
      %Calculate one-tangles for CR-type evolution operator.
      %Input: Nnuc: # of nuclear spins
      %       phi0,phi1: Rot angles
      %       n0,n1: Rot axes
      
          D  = Nnuc+1;
          G1 = @(phi0,phi1,n0n1) (cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2))^2;
          Fs = zeros(1,Nnuc);
          
          for ii=1:Nnuc
          
              Fs(ii)=G1(phi0(ii),phi1(ii),dot(n0(jj,:),n1(jj,1)));
        
          end
          
          epNuc                 = 2/9*(1-Fs);   
          ep_electron           = 1/3 - 1/3^(D)*prod(1+2*Fs);
          obj.one_tangles       = {epNuc,ep_electron};
          
          epmax_Nuc      = 2/9;
          epmax_Electron = 1/3-1/(3^D);
        
          obj.maxBoundPi = {epmax_Nuc,epmax_Electron,[]};
          
      end
      
      function obj=one_tangles_General_U_Zanardi(obj)
      %To calculate one-tangles of arbitrary U. We follow a method based on
      %Zanardi.
      
          U     = obj.Uval;
          [~,d] = size(U);
          d     = log(d)/log(2);
          Omp0  = 1;
          dQ    = 2;
          
          for indx=1:d
              
             Pplus            = 1/2*(eye(2^(2*d))+ArbDimSWAP(indx,indx+d,2*d));
             Omp0             = Omp0*Pplus;
             
          end
          
          Omp0 = Omp0* 1/( dQ+1 )^d;
          
          ep      = zeros(1,d);
          epNames = cell(1,d);
          
          for indx=1:d
              
              Pminus        = 1/2*(eye(2^(2*d))-ArbDimSWAP(indx,indx+d,2*d));
              ep(indx)      = 2 * trace(  kron(U,U) * Omp0 * kron(U',U') * Pminus) ;
              epNames{indx} = strcat(num2str(indx),'|rest');       
        
              
          end

          obj.one_tangles = {real(ep),epNames};
            

      end
      
      function obj=one_tangles_General_U_byparts(obj)
      %To calculate the one-tangles of general U, based on the partition
      %method of the paper on multipartite entanglement.
      
          U  = obj.Uval;
          n  = log(size(U,2))/log(2);
          U  = U(:);
          U  = U/norm(U);
          dQ = 2;
           
          byparts_primal    = {1:n};               %Consider only separating 1 subsystem from rest.
          byparts_secondary = all_byparts(n);      %Get all bypartitions of the secondary system
           
          expr   = zeros(1,n);
          labels = cell(1,n);
          cnt    = 0;
          
          for jj=1:length(byparts_primal)  %Loop over primal system bypartitions

              primal_partitions = byparts_primal{jj}; 
              [~,cols_prime]    = size(primal_partitions); 
          
              for kk=1:cols_prime  
              
                  cnt = cnt+1;
                  p   = primal_partitions(:,kk);  

                  %Label the p|q: 
                  labels{cnt} = strcat(num2str(p),'|rest');
              
                  for ii=2:(n+1) %Will add the ii=1 by hand later %loop over secondary system bipartitions
                  
                      xprime       = byparts_secondary{ii};
                      [~,cols_sec] = size(xprime); %dimension of secondary bipartition

                      for ll=1:cols_sec
                      
                          expr(cnt) = expr(cnt) +...
                                      trace(   ptrace( U*U',   [p.' (xprime(:,ll)+n).'] , repmat(2,[1,2*n])  )^2 );  %tracing the p and the x'

                      end
                      
                  end
                  
                  %Add the empty bipartition of the secondary system
                  
                  expr(cnt) = expr(cnt) + trace( ptrace(U*U', p.' , repmat(2,[1,2*n]) )^2  );
                  expr(cnt) = ( 1- dQ^n/(dQ+1)^n * expr(cnt)  );  %Notice here: we drop the 2 factor, following Zanardi's convention.

              end
              
          end
          
          obj.one_tangles={expr,labels};
          
      end
      
  end
  
    
end