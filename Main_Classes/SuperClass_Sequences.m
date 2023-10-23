classdef SuperClass_Sequences  
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

%Note 1: The Hamiltonian assumes no nuclear-nuclear interactions for now.

%Note: If we do not specify a default value for a property with attributes, the
%matlab initializes as that property as []. So it raises an error if we do
%not provide a default value, e.g. {mustBePositive} is not []
%You can check this if you do isnumeric([])? and matlab raises 1
%but ispositive([]) gives nothing. So it doesnt know if it is positive or
%not and raises the error.


    properties (Constant,Hidden)
        
        c=2*pi*1e-3;
        
    end
    
    
    
    %Attributes of properties
    %public access: unrestricted access
    %private access: access by a class member only (not subclass) 
    %protected access: access by class or subclasses 
    
    %GetAccess -> show (public), do not show to anyone (private) do not
    %show to public but show to class and subclass (protected)
    
    %SetAccess -> allow modification to all (public), 
    %            do not allow modification to anyone, meaning read-only (private or protected)
    
    
    %                  Attributes of methods
    %    Sealed = true (if you try to write a function with the same name in a subclass throw an error)
    %    Static L if it doesn't depend on the object
    
    
    
    properties (GetAccess=public,SetAccess=public)  %it both shows the properties to the public and is writable
        
        N_nuc (1,1) {mustBePositive} = 1 % # of Nuclear spins
        wL (1,1) {mustBeNumeric,mustBePositive}  = 314 %Larmor frequency
        A                                    %parallel HF parameter (needs to be in kHz)
        B                                    %perpendicular HF paramemter (needs to be in kHz)
        s0 (1,1) {mustBeNumeric}    %1st electron spin projection
        s1 (1,1) {mustBeNumeric}    %2nd electron spin projections
        N  (1,1) {mustBePositive}  =10 %number of pulses
        k  (1,1)  {mustBePositive} = 1 %resonance #
        Uevol                       %evolution operator (At the N-th iteration)
        
    end
    
    methods 
        
        %Constructor, to initialize the class
        %This needs to have the same name with the class
        %NOTE: Conversion constant to rad*MHz is implemented internally.
        %Externally need to give values in kHz.
        function  obj  = SuperClass_Sequences(wwL,AA,BB,S0,S1,Nnuc,K,NN)
        
                const = SuperClass_Sequences.c;

                if nargin==0

                      obj.A=60*const; obj.B=30*const;
                      obj.s0=0;
                      obj.s1=-1;
                      
                        
                else
                    obj.A = AA*const;
                    obj.B = BB*const;
                    obj.wL=wwL*const;
                    obj.s0 = S0;
                    obj.s1 = S1;
                    obj.N_nuc = Nnuc;
                    obj.k=K;
                    obj.N=NN;

                end
        
        end
        
    end
    
    
    methods 

        function [Hdiag,T,Tinv] = Get_Ham(obj)
         %Hamiltonian for Nnuc nuclear spins  + the electron (no n-n interactions)
         S0 = obj.s0;  S1 = obj.s1;
         
         
         
         wwL = obj.wL; AA = obj.A; BB = obj.B;
         Nnuc=obj.N_nuc;
         
            Zn=[1 0 ; 0 -1]; X=[0 1 ; 1 0]; Ze=[S0 0 ; 0 S1];

            L1=repmat({eye(2)},[1 1 Nnuc+1]); %this is Zn x Id x Id
            L2=repmat({eye(2)},[1 1 Nnuc]);  %this is Z x Id
            L3=L2;  %this is X x Id
            
            L1{2}=Zn;   L2{1}=Zn;   L3{1}=X;
            
            
            f=@(cellA,k) circshift(cellA,k);  %Shifts element by k positions to the right

                H0=0;
                
                          for ii=1:Nnuc
                              H0=H0+wwL/2*superkron(f(L1,ii-1))+...    
                                  AA(ii)/2*kron(Ze,superkron(f(L2,ii-1)))+...
                                  BB(ii)/2*kron(Ze,superkron(f(L3,ii-1)));
                          end


                [T,~]=eig(H0); % V^(-1)*H*V=H_diagonal
                Tinv=inv(T);
                Hdiag=Tinv*H0*T;


               
        end
      
    end
    
    
    methods %Define standard sequences: CPMG, XY2, UDD
        
        
        function obj = CPMG(obj,varargin)
            %UDDn sequence, Syntax: UDD(n,time(s),iters(optional))
            %If single time is given, then we have
            %one CPMG sequence of N iterations.
            %If multiple times are given, then we are composing CPMG
            %sequences with different times and iterations per CPMG
            %sequence.
           
            Nnuc = obj.N_nuc;
            
            if nargin==2  %we also provide the time
                
                t  = varargin{1};
                NN = obj.N;
              
            elseif nargin==3 %We provide times and iterations
                
                t     = varargin{1};
                NN    = varargin{2};
                
                if length(NN)~=length(t)
                   error('Attempted to compose CPMG sequences, but length of times is not equal with length of iters.') 
                end                
                
            
            end
            
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0];  XI=kron(X,L);  U_pi=expm(-1i*pi/2*XI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal = 1;
            
            for jj=1:length(t)
                
                Ut     = UT(t(jj)/4);
                UU     = (T*Ut*Tinv) * U_pi * (T*Ut*Ut*Tinv) * U_pi * (T*Ut*Tinv);
                Utotal = UU^NN(jj)*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end
        
        function obj = UDD(obj,varargin)
            %UDDn sequence, Syntax: UDD(n,time)
            
            wwL = obj.wL;  AA = obj.A;     BB = obj.B;  S0 = obj.s0;   S1 = obj.s1;  NN = obj.N;
            
            Nnuc = obj.N_nuc;  
            
            
            kk = obj.k;
            
            if nargin==3  %we have provided both n and time
                n = varargin{1};
                t = varargin{2};
            elseif nargin==2 %we have provided only n
                n = varargin{1};
                
            
                 if ~isempty(kk)
              Om0 = sqrt(  (wwL+S0*AA(1))^2 +(S0*BB(1))^2 );
              Om1 = sqrt(  (wwL+S1*AA(1))^2 +(S1*BB(1))^2 );
              Om = Om0 + Om1;
              
              t= 4*pi*(2*kk-1)/(Om);  
                 else 
                     error('You did not provide a resonance time. Either call UDD(obj,n,time) or provide a resonance time.')
                 end
            end
            
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0];  XI=kron(X,L);  U_pi=expm(-1i*pi/2*XI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
                   
                   product=1;
                   expon2= ( sin( pi*(n+1)/(2*n+2))  )^2 - ( sin(pi*(n)/(2*n+2) ) )^2 ;
                   
                      for qq=1:n

                          expon1=( sin( pi*qq/(2*n+2) ) )^2 - ( sin( pi*(qq-1)/(2*n+2)) )^2;

                               if (-1)^n==1   %even
                                   U= T* UT(t*expon1) *Tinv  *U_pi;
                               else  %odd
                                   U= T* UT(t/2 *expon1) *Tinv  *U_pi;
                               end

                          product=product*U;
                      end

                      
                    if (-1)^n==1 %n is even
                          Uend=T* UT(t*expon2) *Tinv;
                          UU=product*Uend;
                    else  %n is odd
                          Uend= T* UT(t/2*expon2) *Tinv;
                          product=product*Uend;      
                          UU=product^2;    %one more time here 
                         
                    end
                    
                     obj.Uevol=UU^NN;
            
        end
        
        function obj = XY2(obj,varargin)
            %CPMG sequence, Syntax: CPMG(time(s),iters (optional))
            Nnuc = obj.N_nuc;
            
            if nargin==2  %we also provide the time
                
                t  = varargin{1};
                NN = obj.N;
              
            elseif nargin==3
                
                t  = varargin{1};
                NN = varargin{2};
                
                if length(NN)~=length(t)
                   error('Attempted to compose CPMG sequences, but length of times is not equal with length of iters.') 
                end                
                
            end
            
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X     = [0 1; 1 0]; 
            Y     = [0 -1i; 1i 0]; 
            YI    = kron(Y,L);  
            XI    = kron(X,L);  
            U_piX = expm(-1i*pi/2*XI); 
            U_piY = expm(-1i*pi/2*YI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal=1;
            for jj=1:length(t)
                
                Ut      = UT(t(jj)/4);
                UU      = (T*Ut*Tinv) * U_piY * (T*Ut*Ut*Tinv) * U_piX * (T*Ut*Tinv);
                Utotal  = UU^(NN(jj))*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end

        function obj = CPMG_error(obj,varargin)
            %CPMG sequence, Syntax: CPMG(time(s),iters (optional),systematic_pulse_error)
            
            Nnuc = obj.N_nuc;
            
            if nargin==3  %we also provide the time
                
                t         = varargin{1};
                NN        = obj.N;
                pulse_err = varargin{2};
              
            elseif nargin==4
                
                t         = varargin{1};
                NN        = varargin{2};
                pulse_err = varargin{3};
                
                if length(NN)~=length(t)
                   error('Attempted to compose CPMG sequences, but length of times is not equal with length of iters.') 
                end
              
            end
            
             
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0];  XI=kron(X,L);  U_pi=expm(-1i*(pi+pulse_err*pi)/2*XI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal=1;
            for jj=1:length(t)
            
                Ut = UT(t(jj)/4);
                UU = (T*Ut*Tinv) * U_pi * (T*Ut*Ut*Tinv) * U_pi * (T*Ut*Tinv);
                
                Utotal  =  UU^(NN(jj))*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end
        
        function obj = XY2_error(obj,varargin)
            %XY2 sequence, Syntax: XY2(time(s),iters(optional),systematic_pulse_error) 
            Nnuc = obj.N_nuc;
            
            if nargin==3  %We also provide the time
                
                t         = varargin{1};
                NN        = obj.N;
                pulse_err = varargin{2};
                
            elseif nargin==4
                
                t         = varargin{1};
                NN        = varargin{2};
                pulse_err = varargin{3};
              
            end
             
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0]; Y=[0 -1i; 1i 0]; YI=kron(Y,L);  XI=kron(X,L);  
            U_piX=expm(-1i*(pi+pulse_err*pi)/2*XI); U_piY=expm(-1i*(pi+pulse_err*pi)/2*YI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal=1;
            for jj=1:length(t)
                
                Ut=UT(t(jj)/4);
                UU= (T*Ut*Tinv) * U_piY * (T*Ut*Ut*Tinv) * U_piX * (T*Ut*Tinv);
                
                Utotal  = UU^(NN(jj))*Utotal;
            end
            

            obj.Uevol = Utotal;
            
        end
        
        function obj = CPMG_error_variance(obj,varargin)
            %CPMG sequence, Syntax: CPMG(time(s),iters(optional),standard_deviation)
            
            Nnuc = obj.N_nuc;
            
            if nargin==3  %we also provide the time
                
              
                t        = varargin{1};
                NN       = obj.N;
                variance = varargin{2};
                
            elseif nargin==4
                
                t        = varargin{1};
                NN       = varargin{2};
                variance = varargin{3};
                
            end
            
            
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0];  XI=kron(X,L);  U_pi=@(r) expm(-1i*(pi+r*pi)/2*XI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            
            Utotal=1;
            
            for kk=1:length(t)
                
                Ut=UT(t(kk)/4);
                Utemp=@(r) (T*Ut*Tinv) * U_pi(r(2)) * (T*Ut*Ut*Tinv) * U_pi(r(1)) * (T*Ut*Tinv);
                UU=1;
                
                for jj=1:NN(kk)
                
                    r=normrnd(0,variance,1,2);
                    UU=Utemp(r)*UU;
                
                end              
                
                Utotal =  UU*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end
        
        function obj = XY2_error_variance(obj,varargin)
            %XY2 sequence, Syntax: XY2(time(s),iters (optional),standard_deviation) 
            
            Nnuc = obj.N_nuc;
            
            if nargin==3  %We also provide the time
                
                t        = varargin{1};
                NN       = obj.N;
                variance = varargin{2};
              
            elseif nargin==4

                t        = varargin{1};
                NN       = varargin{2};
                variance = varargin{3};
                
            end
             
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0]; Y=[0 -1i; 1i 0]; YI=kron(Y,L);  XI=kron(X,L);  
            U_piX=@(r) expm(-1i*(pi+r*pi)/2*XI); U_piY=@(r) expm(-1i*(pi+r*pi)/2*YI);
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal=1;
            
            for kk=1:length(t)
                
            
                Ut=UT(t(kk)/4);
                
                Utemp=@(r) (T*Ut*Tinv) * U_piY(r(2)) * (T*Ut*Ut*Tinv) * U_piX(r(1)) * (T*Ut*Tinv);
                UU=1;
                
                
                for jj=1:NN(kk)

                    r=normrnd(0,variance,1,2);
                    UU=Utemp(r)*UU;

                end
                
                Utotal = UU*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end

        function obj = XY2_error_FULL(obj,varargin)
            %Composition of XY2 sequences, Syntax: XY2(time(s),iters (optional),pulse_error)
            %The errors are fixed according to errors reported in deLange
            %(Delft thesis)
            
            Nnuc = obj.N_nuc;
            
            x_error = 0.005;
            z_error = 0.05;
            
            if nargin==3 
                
              t         = varargin{1};
              NN        = obj.N;
              pulse_err = varargin{2};
              
            elseif nargin==4

              t         = varargin{1};
              NN        = varargin{2};
              pulse_err = varargin{3};
            
            end
             
            [Hdiag,T,Tinv] = obj.Get_Ham;
            
            L = repmat({eye(2)},[1, 1, Nnuc]);  L = superkron(L);
           
            X=[0 1; 1 0]; Y=[0 -1i; 1i 0]; YI=kron(Y,L);  XI=kron(X,L);  
            Z=[1 0 ;0 -1];
            
            ZI=kron(Z,L);
            
            U_piX= expm(-1i*(pi+pulse_err)/2*...
                                       (sqrt(1-z_error^2)*XI+z_error*ZI) );
                                   
            U_piY= expm(-1i*(pi+pulse_err)/2*...
                                       (x_error*XI+sqrt(1-z_error^2-x_error^2)*YI+z_error*ZI) );
            
           
            
            UT =@(x)  expm(-1i*Hdiag*x);
            
            Utotal=1;
            
            for kk=1:length(t)
                
                Ut     =  UT(t(kk)/4);
                UU     = (T*Ut*Tinv) * U_piY * (T*Ut*Ut*Tinv) * U_piX * (T*Ut*Tinv);
                Utotal = UU^(NN(kk))*Utotal;
                
            end
            
            obj.Uevol = Utotal;
            
        end
         
    end
    
    
    
end