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
                %Method followed based on Ref:
                %Makhlin, Y. Nonlocal Properties of Two-Qubit Gates and Mixed States, 
                %and the Optimization of Quantum Computations. Quantum Information Processing 1, 243–252 (2002).
             
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
        

            
   end
   
   
   
   
   
   
end
  