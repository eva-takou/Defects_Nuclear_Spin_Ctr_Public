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
         
      function Ek=Get_Kraus(U,Nnuc,Target_Spins)
          
          %Target_Spins: Subsystem indices for nuclear spins to keep.
          
          if any(Target_Spins>Nnuc)
              
              error('The nuclear spin indices exceed the # of total nuclei.')
              
          end
          
          K    = length(Target_Spins); %# of subsystems in target nuclear subspace
          Kall = K+1;                  %# of susbystems in target subspace (including electron)
          
          %------- Parameters for the environment -------------------------
          
          Nenv      = Nnuc-K;                   %nuclei in environment
          dim_E     = 2^Nenv;                   %dimension of environment
          id_E      = eye(dim_E);               %Identity
          ket0_E    = zeros(dim_E,1);          
          ket0_E(1) = 1;                        %initial state of environment
          e0_E      = kron(eye(2^Kall),ket0_E); %Promote to total space
          
          %--- Bring all target spins to first positions ------------------
          
          Target_Spins = sort(Target_Spins)+1; %Add 1 because 1st subsystem is electron
         
          if any(Target_Spins>Kall) %Need to permute

              for jj=2:length(Target_Spins)+1

                  if Target_Spins(jj-1)~=jj

                      SWAP = ArbDimSWAP(jj,Target_Spins(jj-1),Nnuc+1); 
                      U    = SWAP*U*SWAP'; %Permute to jj-th position

                  end

              end              
              
          end
          
          %"Trace-out" all the final subsystems
          
          Ek = zeros(2^(Kall),2^(Kall),dim_E);
          
          for k=1:dim_E
              
              ek = kron(eye(2^(Kall)),id_E(k,:)); %For systems we do not trace out, we need to put Id.
              Ek(:,:,k)=ek*U*e0_E;
              
          end
          
          
          %Check completeness of Kraus:
          compl=0;
          
          for k1=1:dim_E
              
              
              compl = compl + Ek(:,:,k1)'*Ek(:,:,k1);
                  
              
          end
          
          imag_compl = imag(compl);
          
          if all(imag_compl<1e-9)
              
              compl=real(compl);
              
          end
          
          if norm(compl-eye(2^(Kall)))>1e-9
              
              error('Completeness is not satisfied')
              
          end
          
          
      end
      
      function Infid=Get_Infid_From_Kraus(U0,Ek,Target_Spins)
          
          Kall = length(Target_Spins)+1;
          
          if log2(length(U0))~=Kall
              
             error('Dimension of target gate does not match the # of target susbystems.') 
             
          end
          
          expr1=0;
          expr2=0;
          
          for k=1:size(Ek,3)
             
              Mk=U0'*Ek(:,:,k);
              
              expr1=expr1+Mk'*Mk;
              expr2=expr2+abs(trace(Mk))^2;
              
              
          end
          
          m=2^Kall;
          
          Infid = 1 - 1/(m*(m+1))*(trace(expr1)+expr2);
          
          
          
          
          
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