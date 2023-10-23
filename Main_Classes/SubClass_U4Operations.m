classdef SubClass_U4Operations < SuperClass_Sequences
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
    
%Note 1: This subclass does operations only on unitaries U(4). It is sort of a 
%        characterization class for any unitary U. It borrows all the sequences from the
%        SuperClass_Sequences and it can decompose the U operator into the
%        local and non-local part. It can further draw the circuit.

%This class can only do operations on U(4) matrices

%Note 1: This is a subclass that inherits properties and methods from the Sequences superclass

%Note 2: What we include here: 
%       1)KAK decomp (only for unitary 4 x 4),  
%       3)GET_KRAUS_AND_FID algorithm
%       5)Entangling measures (unitary and non-unitary on 4x4)

%Note 3: To be updated: Write a script that generalizes the one-tangles to arbitrary dimensions
%Note 4: To be updated: Write a script that simplifies the KAK circuit
%Note 5: To apply a particular method of a superclass onto a subclass
%        we do it by superMethod@MySuperClass(obj,superMethodArgs)


   properties (GetAccess=public, SetAccess=protected) 
       
      Ep           %Entangling power of 4x4 unitary
      Makhlin      %Makhlin invariants, G1,G2, and gj,
      cj           %cj-Weyl chamber coefficients
      Local_Gates  %U(2) local gates (from KAK decomposition) 
      EntClass     %EntClass (the middle non-local part from KAK decomposition)
      angles       %Rotation angles for CR-type sjj\otimes Rnj(\phi_j) 
      axes         %Rotation axes for CR-type sjj\otimes Rnj(\phi_j) 
      Nmax         %Iterations of sequence where Ep becomes maximal.
      
   end
   
   properties (Constant,Hidden)
       
        %There are many conventions for this tranformation into the
        %Bell-basis.
        %I am following the one from Ref: Nonlocal properties of two-qubit gates and mixed states
        %and optimization of quantum computations
        
      QQ=1/sqrt(2)*[1   0  0   1i ; ...
                    0  1i  1   0  ; ...
                    0  1i -1   0  ; ...
                    1   0  0  -1i];
      
      x=[0 1 ; 1 0];
      
      y=[0 -1i ; 1i 0];
      
      z=[1 0 ; 0 -1];
      
      xx= [   0     0     0     1;...
              0     0     1     0;...
              0     1     0     0;...
              1     0     0     0];
      
      yy=[ 0     0     0    -1;...
           0     0     1     0;...
           0     1     0     0;...
          -1     0     0     0];
      
      
      zz=[1     0     0     0;...
          0    -1     0     0;...
          0     0    -1     0;...
          0     0     0     1];
       
       
       
   end

   
   methods  %Constructor
       
       %Call superclass constructor
       function obj = SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N)
           
           obj = obj@SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,k,N);
           
       end
        
   end
   
   %A subclass method can call superclass method only if both have the same
   %name. From the subclass, reference the superclass method as:
   %SuperMethod@SuperClass(obj,superMethod Args)
   
   methods %Inherit here all the from SuperClass
              
       function obj=CPMG(obj,varargin)
          
           if nargin==2
                
               time = varargin{1};
               obj  = CPMG@SuperClass_Sequences(obj,time);
               
           elseif nargin==3
               
               times = varargin{1};
               iters = varargin{2};
               
               obj=CPMG@SuperClass_Sequences(obj,times,iters);
               
           end
           
       end

       function obj=UDD(obj,n,time)
           
          obj=UDD@SuperClass_Sequences(obj,n,time);
          
       end
       
       
       function obj=XY2(obj,varargin)
           
           if nargin==2
          
               time = varargin{1};
               obj  = XY2@SuperClass_Sequences(obj,time);
               
           elseif nargin==3
               
               times = varargin{1};
               iters = varargin{2};
               obj   = XY2@SuperClass_Sequences(obj,times,iters);
               
           end
           
       end

       
   end
   
   methods (Sealed = true,Static) %Scripts for KAK decomposition
       
            %Simultaneous Diagonalization via Givens Rotations
            function [U,A]=cyclic_Jacobi_Givens_Diag(A)
             %Requires: input A in cell (A{1} and A{2}). 
             %Outputs: Simultaneous diagonalization matrix and diagonal
             %matrices A in cell array again.

                        epsilon = sqrt(eps);

                        %cost function of off-diagonal terms
                        off = @(X) norm( X{1},'fro')^2 - sum(diag(X{1}).^2) + norm(X{2},'fro')^2 - sum(diag(X{2}).^2);

                        h = @(X,p,q) [X(p,p)-X(q,q) , X(p,q)+X(q,p)];

                        h_comp = @(h,indx) h(indx);

                        [dim,~] = size(A{1}); 

                        In = eye( dim );
                        proj = @(i,j)  In(:,i)*In(:,j)'; %|i><j|

                        %Givens matrix:

                        %Its dim, has to be of the dim of A{1} or A{2}.

                        G = @(p,q,c,s) In + (c-1)*proj(p,p) + (c-1)*proj(q,q) ...
                            + (s)*proj(p,q) - s * proj(q,p);

                        U = eye(dim);

                        err = off(A);


                        while err>epsilon

                            for p=1:dim-1 %p<q
                                for q=p+1:dim

                                    ton = h_comp( h(A{1},p,q),1  ) + h_comp( h(A{2},p,q),1  );

                                    toff = h_comp( h(A{1},p,q),2  ) +  h_comp( h(A{2},p,q), 2  );
                                    theta = 0.5*atan2(toff, ton+sqrt(ton*ton+toff*toff));

                                    cd_star = cos(theta);
                                    sd_star = sin(theta);

                                    %update A's
                                    A{1}= G(p,q,cd_star,sd_star)*A{1}*G(p,q,cd_star,sd_star)';
                                    A{2}= G(p,q,cd_star,sd_star)*A{2}*G(p,q,cd_star,sd_star)';
                                    %update U
                                    U = G(p,q,cd_star,sd_star)*U;
                                    %update the error
                                    err=off(A);
                                end
                            end
                        end


                        U = U';
                        %Output: the diagonalizing matrix U,
                        %and the diagonalized version A{1} and A{2}


                        %Pass it as U=U' because if we dont, then it is %P*G*P' = Gdiag
                        %but we want P'*G*P, to be consistent with Eckart's convention.



                        




            end
            
            %Diagonalization via SO4 svd matrices
            function [U,Dtot,V]=SO4_diag(UB)
             %Diagonalizes the real and imaginary part Bell-Basis matrix with SO4 e-vecs
             %using Eckart's theorem. It calls the Jacobi-Givens angles

                        %Take the real and imaginary parts:
                        A = (UB+conj(UB))/2;
                        B = (UB-conj(UB))/(2i);

                        %Do SVD on the real part:

                        [LA,S,RA] = svd(A);   %A=U*S*V': Guaranteed to produce SO(4) matrices since we are working with real matrix


                        [~,dim_U]=size(UB);

                        %Find the rank of A

                        tol = dim_U*S(1,1)*eps;
                        rank_A =1;
                        j=2;
                        done='false';


                        while j<dim_U+1 && strcmp(done,'true')~=1

                           if S(j,j)>tol
                               rank_A=rank_A+1;
                           else
                               done='true';
                           end
                           j=j+1;


                        end


                        %Now, if the n=rank(A)<dim_U, from Eckarts theorem:

                        %D = [D_1{n,n} D_2{dim_U-n,n};...
                        %     D_3{n,dim_U-n} D_4{dim_U-n,dim_U-n}]= > D = [D1 O_2 ; O_3 O_4]

                        %And it further says that the application of the SVD of A onto B,
                        %i.e. LA' B RA = B_transf, gives:
                        %B_transf=  [G K ;
                        %            L H]


                        %If A is full-rank then D=D1_{n,n}
                        %Then B_transf  = G

                        %If A is not full rank, then we have the above form for K and L
                        %where K needs to be 




                        Btransf = LA'* B * RA;

                        if rank_A<dim_U

                            D = S(1:rank_A,1:rank_A);
                            D2 = zeros( dim_U - rank_A ,  rank_A);
                            D3 = zeros( rank_A         ,  dim_U-rank_A);
                            D4 = zeros(  dim_U-rank_A   ,  dim_U-rank_A);


                            G = Btransf(1:rank_A,1:rank_A);          %Now, G should be of rank(A)
                            %note that if rank_A is not dim_U
                            %then we need a givens rotation of smaller dim


                            K = Btransf(1 : dim_U - rank_A , rank_A+1 : dim_U);  %this is upper right (Should be 0.)
                            L = Btransf(rank_A+1 : dim_U , 1 : dim_U-rank_A );   %this is lower left  (Should be 0.)
                            H = Btransf(rank_A+1 : dim_U   , rank_A+ 1 : dim_U); %non-zero



                        elseif rank_A == dim_U

                            D=S;
                            G = LA' * B * RA;

                            %K, L , H are empty
                            %so in that case U =  LA*P;
                            % and            V =  RA*P;
                            %because,  LA' * A * VA = D
                            %          P'*(LA'*A*RA)*P= D -> U = LA*P and V=RA*P

                            %and LA'*B*VA = Btransf->  P' * LA' * B * VA * P = G
                            %                             U' * B * V = G



                        end



                        %First, do further checks: D and G have to commute and they have to be
                        %orthogonal matrices



                        if norm(G*D-D*G)>1e-9
                            error('the input matrices do not commute.')
                        end

                        if norm(G-G.')>1e-9
                            error('the matrix G is not orthogonal.')
                        end
                        if norm(D-D.')>1e-9
                            error('the matrix D1 is not orthogonal.')
                        end

                        %Store the det(G) and det(D1) (if they are not special orthogonal we might
                        %need to do something.

                        %det_G = det(G);
                        %det_D1 = det(D1);


                        %Now, next step according to Eckart, is to diagonalize simultaneously G and
                        %D1.

                        %Call the matrix that does the simultaneous diagonalization P

                        [P,Ds]=SubClass_U4Operations.cyclic_Jacobi_Givens_Diag({G,D});  %P*G*P' = Gdiag

                        Gdiag = Ds{1};
                        D = Ds{2};

                        %Next, Eckart says, if we have full rank, then from Btransf
                        %K is empty, L is empty and H is empty.




                        if rank_A == dim_U

                            U = LA*P;
                            V = RA*P;

                            Dtot = D + 1i  * Gdiag;


                        elseif rank_A<dim_U
                            %now H will not be diagonal
                            %so we further need to diagonalize it
                            %find its SVD with matrices Q, R: Q'*H*R= Hdiag

                            [Q,Hdiag,R]= svd(H);  %H=Q*S*R' - > Q' H R  = Hdiag

                            %and then, U = [P,0;0 Q] and V = [P, 0;0 R]
                            %times the SVD LA and RA

                            U = LA*[P , D2; D3 , Q];
                            V = RA*[P , D2; D3 , R];

                            %And the diagonalized will be:
                            B_transf_diag = [Gdiag, D2; D3, Hdiag];    
                            Dtot = [D D2 ; D3 D4]+1i*B_transf_diag;
                        end




                        if det(U)<1e-3
                            U(:,1)=-U(:,1);
                            Dtot(1,1)=-Dtot(1,1);
                        end


                        %some fucking bug here

                        if det(V)<0

                            V(:,1)=-V(:,1);
                            Dtot(1,1)=-Dtot(1,1);

                        end









                        end
                        
            %SO4 to SU2 x SU2 decompositon
            function [A,B]=SO4ToSU2xSU2(q)
            %SO4 to SU2 x SU2 decompositon of the single qubit gates. They
            %are returned to the callback function in the computational
            %basis.
            
            
                        
                        M=SubClass_U4Operations.QQ;
                       


                        %Check if Orthogonal matrix was passed
                        if norm(q*q.'-eye(4))>1e-8
                            error('input matrix for SO(4)->SU(2)xSU(2) is not orthogonal')
                        end
                        %Check if Special orthogonal matrix was passed
                        if abs(det(q)-1)>1e-8
                           error('input matrix for SO(4)->SU(2)xSU(2) is not special') 
                        end


                        U = M*q*M'; %inverse transformation



                        %Start by assuming that 
                        %U =
                        %[a11*B, a12*B]
                        %[-a12'*B, a11'*B]
                        %Will check validity of this assumption at the end


                        A = zeros(2,2);
                        




                        a11_sq = (U(1:2,1:2)*U(3:4,3:4)');
                        a11_sq=a11_sq(1,1);
                        a12_sq = -(U(1:2,3:4)*U(3:4,1:2)');
                        a12_sq = a12_sq(1,1);

                        a11_a12H = (U(1:2,1:2)*U(1:2,3:4)');
                        a11_a12H = a11_a12H(1,1);
                        A(1,1) = sqrt(a11_sq);
                        A(1,2) = sqrt(a12_sq);


                        %sqrt ambiguous in its sign, 
                        %so the values just
                        %assigned to A(1,1) and
                        %to A(1,2) may have the wrong sign.
                        %Fix relative sign:
                        if(abs(A(1,1)*A(1,2)'-a11_a12H)>1e-9)
                            A(1,2)=-A(1,2);
                        end


                        A(2,1) = -A(1,2)';
                        A(2,2) = A(1,1)';
                        if(abs(A(1,1))>abs(A(1,2)))
                            B = U(1:2,1:2)/A(1,1);
                        else
                            B = U(1:2,3:4)/A(1,2); 
                        end


                        %the moment of truth

                        if norm(kron(A,B)-U)>1e-9
                            error("the error in SO(4)->SU(2)xSU(2) is too large")
                        end


            end
            
            
            function b = diag_to_class_vec(dd)
            %Find the cj-Weyl chamber (canonical) coefficients of the input matrix.
            %first components are cx,cy,cz and last component is the phase.


                        %Here we take the Dx + i Dy part that we found from SO4_diag, i.e. the Dtot
                        %part.


                        %find the entangling class vector from the diagonal components of the
                        %diagonal matrix
                        %and the global phase stored in the 4th component of cvec

                        dd = diag(dd);

                        t1 = angle(dd(1));
                        t2 = angle(dd(2));
                        t3 = angle(dd(3));
                        t4 = angle(dd(4));


                        Gamma = [1 1 -1 1;...
                                 1 1  1 -1;...
                                 1 -1 -1 -1;...
                                 1 -1 1 1];

                        k = inv(Gamma)*[t1 ; t2; t3 ; t4];

                        k0 = k(1);  %0
                        kx = k(2);  %pi/4
                        ky = k(3);  %0
                        kz = k(4);  %0



                        cvec=[kx;ky;kz;k0];

                        %Now, make sure we get the right entangling class


                        %b will satisfy
                        %pi/2> bx>=by>=bz>=0
                        %pi/2>=bx +by
                        %if bz=0 then pi/4>=bx











                        b = cvec;
                        phi=pi/2;


                        %Shift repeatedly kx such that kx will be in pi/2
                        %Similarly, shift ky and kz into [0 pi/2)



                        for jj=1:3 
                            if b(jj)<0

                                while b(jj)<0
                                    b(jj)=b(jj)+phi;
                                end
                            end

                            if b(jj)>=phi
                                while b(jj)>=phi
                                    b(jj)=b(jj)-phi;
                                end
                            end



                        end

                        %make kx>=ky>=kz
                        b(1:3) = sort(b(1:3),'descend') ; 

                        %If kx+ky>pi/2
                        %Transform k into (pi/2-ky,pi/2-kx,kz)

                        %This can be done by 1 swapping kx with ky
                        %then reverse kx -> -kx and ky->-ky
                        %then shift kx by pi/2 and then shift ky by pi/2

                        if b(1)+b(2)>phi
                            old_b1=b(1);

                            b(1)=b(2);
                            b(2)=old_b1;

                            b(1)=phi-b(1);
                            b(2)=phi-b(2);



                            %at this point b(1)>b(2) but b(3) may be larger than b(1) || b(2)


                            if b(3)>b(2) || b(3)>b(1)
                              b(1:3)=sort(b(1:3),'descend');  %sort again in descending order such that kx>=ky>=kz

                            end

                        end

                        if b(3)<1e-14 && b(1)>pi/4
                            b(1)=phi-b(1);
                        end







                        end
       
   end
   
   methods %Here add the extra capabilities
            
            %This is a purely numerical evaluation.
            %Produces right rotation angles and axes;
            %For the UDD4 however, we might get the gate
            %U=s00 x Rn(phi0) - s11 x Rn(phi1)
            function obj = Rot_Angles_And_Axes(obj)
                %Get rotation angles and axes of target spin. (phi0,phi1,n0,n1)
                X=[0 1; 1 0];
                Z=[1 0; 0 -1];
                Y=[0 -1i;1i 0];
                s00=[1 0 ; 0 0];
                s11=[0 0 ; 0 1];
                
                
                %it should take the obj.Uval to find the angles numerically.
                U=obj.Uevol;  %note this will calculate the rot. angle for N pulses
                
                U12  = U(1:2,1:2);
                U34  = U(3:4,3:4);
                phi0 = real(2*acos( trace(U12/2 ) ));
                phi1 = real(2*acos( trace(U34/2 ) ));
                
                %After N iterations, we might get a rot. angle>pi
                
                flag1=0; 
                flag2=0;
                
                %If we exceed the range [0,pi] do this:
                if phi0>pi
                    phi0=phi0-2*pi;
                    
                end
                
                if phi1>pi
                    phi1=phi1-2*pi;
                    
                    
                end
                
                %I am allowed to do the above, since there is nothing wrong
                %with shifting by multiples of 2pi.
                
                %The dot product would then become:
                %sin(phj/2-pi)=sin(phj/2)cos(pi) = - sin(phj/2)
                
                %Now, if phj<0 I would get a (-)
                %But if I shift phj->-phj such that now phj>0
                %then i would get sin(phj/2) w/o a minus sign.
                %so to take this into account, I would have to 
                %make n0 -> -n0 and n1 ->-n1
                
                
                if phi0<0
                    phi0=-phi0;  
                    flag1=1;
                end
                
                if phi1<0
                    phi1=-phi1;
                    flag2=1;
                end
                
                
                        Opers ={X,Y,Z};
                        
                        nj =@(oper,UU,phj) -imag( 1/(2*sin(phj/2))*trace(oper*UU)          );
                        
                       n0=zeros(1,3);
                       n1=n0;
                        for ii=1:3
                           
                            n0(ii)=nj(Opers{ii},U12,phi0);
                            n1(ii)=nj(Opers{ii},U34,phi1);
                            
                        end
                        
                  if flag1==1    
                  n0=-n0; 
                  end
                  
                  if flag2==1
                      n1=-n1;
                  end
                  
                   
                  obj.axes={n0,n1};
                  
                  %do a test here to make sure that we get the right
                  %evolution operator.
                  
                  R = @(a,nx,ny,nz) expm(-1i*a/2*(nx*X+ny*Y+nz*Z) );
                  
                  Utest1 = kron(s00,R(phi0,n0(1),n0(2),n0(3)))+kron(s11,R(phi1,n1(1),n1(2),n1(3)));
                  Utest2 = kron(s00,R(phi0,n0(1),n0(2),n0(3)))-kron(s11,R(phi1,n1(1),n1(2),n1(3)));
                  
                  %You might obtain the U or -U gate (but global phase
                  %does not matter). Consider the following test condition:
                  
                  
                  Fid=@(U,U0,n) 1/(n*(n+1))*...
                  ( trace( (U0'*U) *(U0'*U)'   ) +  abs(  trace( U0'*U )   )^2);
                  
                    test1=Fid(Utest1,U,4);
                    test2=Fid(Utest2,U,4);
                  
              
                    if abs(test1-1)>1e-5 && abs(test2-1)>1e-5
                       error('Did not produce the right decomposition.') 
                    else
                        obj.angles={phi0,phi1};
                        
                    end
              
                  

                        
            end
            

            function obj = Expected_Maxima0(obj,maxNum)
               %This function is based on the fact that the time is always chosen such that n0n1~-1. 
               %Whatever the input object, it carries the info for the
               %rotation angle at N iterations.
               %To find all the maxima we need the rotation angle for 1
               %iteration.
               
               
               obj = obj.Rot_Angles_And_Axes;
               
               phi=obj.angles;
               phi0=phi{1}; phi1=phi{2};
               
               %now it might be the case that phi0 and phi1 are very close to 2pi
               
               %This is something i am not 100% sure I should code it like
               %that.
               
               
               if (2*pi-phi0)>0 &&  (2*pi-phi0)<0.2
                   phi0=2*pi-phi0;
               end
               
               if (2*pi-phi1)>0 &&  (2*pi-phi1)<0.2
                   phi1=2*pi-phi1;
               end

               
               Nmaxima = zeros(1,maxNum);
               
               if obj.N==1
                  %the easy case
                  for ll=1:maxNum
                  Nmaxima(ll) = round( (2*ll-1)*pi/(phi0+phi1));
                  end
                  
               else
                  %we know that phi0(N) = acos(cos(Nphi0^(1)))
                  % cos(phi0(N)) = cos(N*phi0^(1))
                  % 1/N*acos(cos(phi0(N))=phi0^(1)
                  
                  phi0_1iter = 1/N*acos(cos(phi0));
                  phi1_1iter = 1/N*acos(cos(phi1));
                  
                  for ll=1:maxNum
                  Nmaxima(ll) = round( (2*ll-1)*pi/(phi0_1iter+phi1_1iter));
                  end

                   
               end
               
               flag=0;
               
               
               %check for instance 5 maxima to make sure that we get the
               %right periodicity; alteratively fit a function
               
               %This is not entirely wrong, but in some cases it is
               %unecessary and we produce exactly the same maxima.
               
               
               if obj.N==1 %it holds only in this case
                        G1=zeros(1,4);
                       for ii=1:4
                           %get the evolution operator in 1 iteration
                           test = SubClass_U4Operations(obj.wL,obj.A,obj.B,obj.s0,obj.s1,1,obj.k,1);
                           test.Uevol = obj.Uevol^Nmaxima(ii);

                           test = test.Makhlin_Inv;

                           G1(ii) = test.Makhlin{1}{1};

                       end
                            
                       G1=real(G1);
                           %If at least one G1 has real part >tol
                           %then proceed with the other case.
                           
                           if abs(G1(1))>5e-3 || abs(G1(2))>5e-3 || abs(G1(3))>5e-3 || abs(G1(4))>5e-3
                               
                               %If we want to display the warning: 
                               %warning('on')
                               
                               warning('The maxima of the entangling power are not correct. Proceeding with the fitting function.') 
                               flag=1;
                               
                           else
                                obj.Nmax = Nmaxima;
                               return %No need to proceed to the rest of the code.

                           end

                       

                   
               end
               
               
               
                      %take up to 200 pulses, and find the entangling power
                      
                      NN = 200;
                      ENT_Power = zeros(1,NN);
                      for ii=1:NN
                          
                          U1 = obj.Uevol;
                          
                          UN = U1^ii;
                          
                          Empty = SubClass_U4Operations(obj.wL,obj.A,obj.B,obj.s0,obj.s1,1,obj.k,ii);
                          Empty.Uevol = UN;
                          
                          Empty = Empty.Makhlin_Inv;
                          
                          ENT_Power(ii) = Empty.Ep;
                          
                          
                          
                      end
                      ENT_Power=ENT_Power/(2/9);
                      NumPulses = 1:NN;
                      
                      ss=fit(NumPulses.', ENT_Power.','fourier1');
                      
                      
                      %It has the form
                      %a0 + a1*cos(w*NumPulses)+b1*sin(w*NumPulses)
                        w = ss.w;
                        a0 = ss.a0;
                        a1 = ss.a1;
                        b1 = ss.b1;
                        
                        %This holds true if b1=0.
                        Maxima_Fitted_Alt=zeros(1,maxNum);
                        if abs(b1)<1e-3
                            for ll=1:maxNum

                                Maxima_Fitted_Alt(ll) = round( 2*ll*pi/w - 1/w*acos((1-a0)/a1));
                            end                        
                        else
                            
                            for ll=1:maxNum
                               Maxima_Fitted_Alt(ll) = round( (2*ll-1)*pi/w - 1/w*atan(b1/a1) ); 
                                
                            end
                            
                            %error('b1~=0. The true maxima could not be obtained.')
                        end
                        
                        
                  if flag==1
                      
                      obj.Nmax = Maxima_Fitted_Alt;
                  else
                      error('There is some bug in the code.')
                      
                  end
                        
                        
            end
            
            
            
            function obj = Expected_Maxima(obj,maxNum)
               %This function is based on the fact that the time is always chosen such that n0n1~-1. 
               %Whatever the input object, it carries the info for the
               %rotation angle at N iterations.
               %To find all the maxima we need the rotation angle for 1
               %iteration.
               
               obj = obj.Rot_Angles_And_Axes;
               
               phi=obj.angles;
               phi0=phi{1}; phi1=phi{2};
               n0n1=dot(obj.axes{1},obj.axes{2});
               
               if obj.N~=1
                  error('I need the angle in one iteration.') 
                   
               end
               
               Iters=zeros(1,2*maxNum);
               
               for jj=1:maxNum
               
               Iters(jj)        = round(1/phi0*(2*jj*pi-2*atan(sqrt(-1/n0n1)))   );
               Iters(jj+maxNum) = round(1/phi0*(2*(jj-1)*pi+2*atan(sqrt(-1/n0n1)))     );
 
               end
                        
               obj.Nmax = sort(Iters);
               
               
            end
            
            %Get only G and gj from here
            function obj=Makhlin_Inv(obj)
             %Outputs the G, gj and EP of the input object   
                   if isempty(obj.Uevol)
                       
                       error('You need first to provide a value for Uevol.')
                       
                       
                   end
                       
                   
                   U=obj.Uevol;  U = U/det(U)^(1/4); 
            
                   Q = SubClass_U4Operations.QQ;
                   
                   
                   UB = Q' *U *Q;  m  = UB.'*UB;
                
                   G1 = (trace(m)^2/(16*det(U))) ;

                   G2 = 1/(4*det(U))*(trace(m)^2 - trace(m*m)  );
                    
                   g1=real(G1); g2 = imag(G1); g3=G2;
                   
                   
                   
                   Gj=[g1,g2,g3];
                   
                   Ent_Power = 2/9*(1-abs(G1));
                   
                
                   obj.Makhlin={{G1,G2,'Gj'},{Gj,'gj'}};
                
                   obj.Ep = Ent_Power;
                
                
                
                
                
                
                
            end
        
            function obj=KAK(Uinput)
            %Uses the operations from Ref: https://www.ams.org/journals/bull/1939-45-02/S0002-9904-1939-06910-3/home.html    
            %And Ref: “Orthonormalapproximate joint block-diagonalization” (Telecom Paris) by C ́edric F ́evotte Fabian J. Theis.    
            %And Ref: https://arxiv.org/abs/quant-ph/0507171
            %Returns: The entangling class (note not actually what corresponds to the cj, but the decomposition we find
            %by the algorithm-sublety here). It returns also the cj coefficients and the local gates.        
                    Uinput=Uinput.Uevol;
                    XX= SubClass_U4Operations.xx; YY= SubClass_U4Operations.yy; ZZ= SubClass_U4Operations.zz;
                    M=SubClass_U4Operations.QQ;
                    
                    
                    phase=det(Uinput)^(1/4);
                    [m,n]=size(Uinput);


                    
                    %Raise errors:
                    if m~=n && m~=4
                        error('input matrix is not 4 by 4.')
                    else
                        dim_U=4;

                    end

                    if norm(Uinput*Uinput'-eye(4))>1e-8
                        error('input matrix is not unitary.')
                    end

                    



                            %Bell-Basis transformation:

                            Uinput =  Uinput /phase;  %SU4 reduction
                            UB     =  M' * Uinput * M;

                            [U,DD,V]=SubClass_U4Operations.SO4_diag(UB);  %U*DD*V' = UB 


                            %So most part is done.
                            %Now, make the checks

                            if norm(U*DD*V'-UB)>1e-4
                                error('Problem with U*D*V^dagger-UB.')
                            elseif norm(U*U.'-eye(4))>1e-6
                                error('U was not found orthogonal')
                            elseif norm(V*V.'-eye(4))>1e-6
                                error('V was not found orthogonal')
                            elseif norm(DD*DD'-eye(4))>1e-6
                                error('DD not unitary')
                            end


                            %Now, get the class of DD

                            %What if instead, we get the class of DDnew = M * DD * M'?, no it doesnt
                            %give the correct thing.

                            b = SubClass_U4Operations.diag_to_class_vec(DD);

                            b=2*b;

                            Ent_Class = expm(1i/2* (b(1)*XX+b(2)*YY+b(3)*ZZ))*expm(1i*eye(4)*b(4));


                            %And now, do the S04 to SU2xSU2


                            [A0,B0]=SubClass_U4Operations.SO4ToSU2xSU2(U);  %(A0 x B0) = M * U * M' %U and V are returned back to the computational basis
                            [A1,B1]=SubClass_U4Operations.SO4ToSU2xSU2(V);  %(A1 x B1) = M * V * M' 


                            %This:  (M'*kron(A0,B0)*M) *  DD * (M'*kron(A1,B1)*M)' = UB
                            %holds all the time do the construction.

                            %Since UB = M' * Uinput * M
                            %it should then hold:

                            % (M'*kron(A0,B0)*M) *  DD * (M'*kron(A1,B1)*M)' = M' *U * M

                            %So this should also hold
                            % kron(A0,B0)*M*DD*M'*kron(A1,B1)' = Uinput


                            %But, then, we need to know in general the relation between M * DD * M'
                            %and the entangling class we find



                            % kron(A0,B0) * M * DD * M' * kron(A1',B1') = Uinput





                            if norm(M'*kron(A0,B0)*M-U)>1e-7
                                error('M^dag* (A0 x B0) M does not give U')
                            elseif norm(M'*kron(A1,B1)*M-V)>1e-7
                                error('M^dag* (A1 x B1) M does not give V')
                            end



                            err = norm(kron(A0,B0)*(M*DD*M')*(kron(A1,B1))' -Uinput );

                            if err>1e-4

                                error('the multiplication of (A0 x B0)*Ent_Class*(A1^dag x B1^dag) does not give the correct result.')


                            end


                            
                            obj.cj = b(1:3);
                            obj.EntClass = M*DD*M';
                            obj.Local_Gates = {{A0,B0},{A1,B1}};


                
            end
            
            function Weyl_Chamber_C(obj)


                    %we want     0<=c3<=c2<=c1<=pi/2 &  pi/2<c1<pi ,  0<=c3<=c2<pi-c1

                    %Thus, we have a tetrahedron with the following vertices
                    
                    
                    
                    
                    figure(2)
                    
                    
                    P =  [0 0 0 ; ...  %Identity, or A x B
                          pi 0 0 ; ... %Identity or A x B
                          pi/2 pi/2 0; ... %DCNOT
                          pi/2 pi/2 pi/2]; %SWAP

                      %Also scale wrt pi

                      P =P/pi;

                      OtherP = [1/2 1/4 0 ;...  %B(Bgate)
                                1/2  0  0 ;...  %L(CNOT)
                                1/4  1/4 1/4;... %P(sqrt(SWAP))
                                1/4  1/4 0;...   %Q
                                3/4  1/4 0;... %M
                                3/4   1/4 1/4;... %N
                                1/2  1/4   1/4]; %R 

                      names={'B','L(CNOT)','P(sqSWAP)','Q','M','N','R'};

                    for ii=1:4
                        if ii~=3
                    scatter3(P(ii,1),P(ii,2),P(ii,3),80,'filled','MarkerFacecolor',[0 0.45 0.74])
                    hold on
                        end

                    end

                    for ii=1:7
                        scatter3(OtherP(ii,1),OtherP(ii,2),OtherP(ii,3),80,'filled','MarkerFaceColor',[1 1 0.07])
                    hold on
                    end

                    hold on
                    scatter3(P(3,1),P(3,2),P(3,3),80,'filled','MarkerFacecolor',[1 1 0.07])
                    hold on
                    %Now we want to draw lines between the points

                    for ii=1:4
                        for jj=1:4
                            if ii~=jj
                                hold on
                                line([P(ii,1),P(jj,1)],[P(ii,2),P(jj,2)],[P(ii,3),P(jj,3)],'linewidth',2)
                                hold on
                            end
                        end
                    end

                    text(P(1,1),P(1,2),P(1,3),'O','FontSize',22,'color',[0 0.45 0.74]) %Identity
                    text(P(2,1),P(2,2),P(2,3),'A1','FontSize',22,'color',[0 0.45 0.74]) %Identity
                    text(P(3,1),P(3,2),P(3,3),'A2(DCNOT)','FontSize',22,'color',[1 1 0.07]) %DCNOT
                    text(P(4,1),P(4,2),P(4,3),'A3(SWAP)','FontSize',22,'color',[0 0.45 0.74]) %SWAP


                    for ii=1:7

                        text(OtherP(ii,1),OtherP(ii,2),OtherP(ii,3),names{ii},'FontSize',22,'color',[1 1 0.07]) %Identity

                    end

                    set(gcf,'color','w')
                    set(gca,'fontsize',22,'fontname','Microsoft Sans Serif')
                    xlabel('c_1/\pi','color','w')
                    ylabel('c_2/\pi')
                    zlabel('c_3/\pi')
                    set(gca,'XColor',[1 1 1]);
                    set(gca,'YColor',[1 1 1]);
                    set(gca,'ZColor',[1 1 1]);
                    set(gcf,'color','k')

                    hold on

                    line([0 0],[0 0],[0 1],'color','k','linewidth',4)
                    hold on
                    line([0 0],[0 1],[0 0],'color','k','linewidth',4)
                    hold on
                    line([0 1],[0 0],[0 0],'color','k','linewidth',4)

                    %line between DCNOT and CNOT
                    hold on
                    line([1/2 1/2],[0 1/2],[0 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between DCNOT and sq(SWAP)
                    hold on
                    line([1/2 1/4],[1/2 1/4],[0 1/4],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between Q and sq(SWAP)
                    hold on
                    line([1/4 1/4],[1/4 1/4],[0 1/4],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between Q and DCNOT
                    hold on
                    line([1/4 1/2],[1/4 1/2],[0 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between Q and L(CNOT)
                    hold on
                    line([1/4 1/2],[1/4 0],[0 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between DCNOT and M
                    hold on
                    line([1/2 3/4],[1/2 1/4],[0 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between CNOT and M
                    hold on
                    line([1/2 3/4],[0 1/4],[0 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between N and M
                    hold on
                    line([3/4 3/4],[1/4 1/4],[1/4 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between N and L(CNOT)
                    hold on
                    line([3/4 1/2],[1/4 0],[1/4 0],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')
                    %line between N and P(sqSWAP)
                    hold on
                    line([3/4 1/4],[1/4 1/4],[1/4 1/4],'color',[1 1 0.07],'linewidth',2,'linestyle','-.')

                    zlim([0,0.6])
                    ylim([0,0.6])
                    axis off
            
                    
                    
                    c_space_Coeffs = obj.cj/pi;
                    
                   
                    
                    
                   
                            
                            
                            scatter3(c_space_Coeffs(1),...
                                c_space_Coeffs(2),c_space_Coeffs(3),200,'filled','markerfacecolor','r')
                            
                       
                        
                       set(gca,'color','k')
                       
                       axis on
                       
                       
                       XXX = [ 0 0 0  ;...
                              1/2 1/2 0 ;...
                              1 0 0;...
                              1/2 1/2 1/2];
                       T=[1 2 3 4];
                       tetramesh(T,XXX,'FaceAlpha',0.3,'facecolor',[0 0.45 0.74])
                       % [0 0.45 0.74],[0 0.45 0.74],[0 0.45 0.74],[0 0.45 0.74]
            end

            
   end
   
   
   
   
   
   
end
  