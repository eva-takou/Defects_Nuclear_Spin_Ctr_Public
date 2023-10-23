function out=ArbDimSWAP(sys1,sys2,NumSys)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Generate an arbitrary SWAP gate between sys1 and sys2 for a total of
%NumSys qubit systems.
%The usual SWAP= [1 0 0 0 
%                 0 0 1 0
%                 0 1 0 0 
%                 0 0 0 1] for example is ArbDimSWAP(1,2,2).
%--------------------------------------------------------------------------

if sys1>NumSys || sys2>NumSys
    
    error('Cannot swap because the susbsystem index is greater than the # of total qubit systems.')
    
end


%to swap sys1 with sys 2, check when the indx of the ket1 \oplus ket2
% for those subsystems is equal to 1. then store that indx. Find the next
% one and store that index as well. then given the ordering just flip
% those position. Then call identity(ordering,:) and this creates the SWAP


%Generate all possible combination of bitstrings 

bitCombs = dec2bin(0:2^NumSys-1) - '0'; %Each row has all the possible bitstrings.  

[dim1,~]=size(bitCombs);

%Take the column of sys1 and sys2 and find when i\oplus j=1
cnt=0;

indx = {};

if NumSys==2
    
   SWAP = [1 0 0 0 ; 0  0 1 0 ; 0 1 0 0 ; 0 0 0 1] ;
   
   out=SWAP;
return


end





for ii=1:dim1
    
    temp = bitxor(bitCombs(ii,sys1),bitCombs(ii,sys2));
    
    if temp==1
        
        %keep searching the next ket that satisfies this:
        
        for jj=ii+1:dim1
            
            if jj~=ii
                if bitxor(bitCombs(jj,sys1),bitCombs(jj,sys2))==1
                    
                    
          %need to make sure that when bitxor = 1 all the rest bits are the same
        
        check1 = bitCombs(ii,:);
        check1(:,[sys1 sys2]) = []; %remove the columns of sys1 and sys2
        
        check2 = bitCombs(jj,:);
        check2(:,[sys1 sys2]) = [];
                  
                    
        if check1==check2
            cnt=cnt+1;
        
            indx{cnt}=[ii,jj];    %store the (i,j) pair that we need to swap
        end
                  
                    
                    
                    
                    
                end
            end
        end
            
        
        
        
    end
   
    
end


Id = eye(2^NumSys);

%generate an ordered array

arrtemp = 1:2^NumSys;    
arr=arrtemp;
%now flip the positions as the indx cell dictates
for ii=1:length(indx) 
    
    arr(indx{ii})=arr(flip(indx{ii}));
    
end



out = Id(arr,:);














end