function [byparts]=distinct_byparts(n)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Generate distinct bypartitions for order n.
%e.g. if n=3 we can have 1|23, 2|13 and 3|12 
%
%Output: byparts: a cell array from which we need to access the partitions
%        columnwise.

b=1:n;
byparts=cell(1,floor(n/2));

for ii=1:floor(n/2)

    byparts{ii}=nchoosek(b,ii).';  %In some cases it introduces more than we need, so run the next section as well.

end        


for ii=1:floor(n/2)

   temp = byparts{ii};

   for jj=1:length(temp)

      rest{ii,jj}=setdiff(b,temp(:,jj));

   end

end            

%number of distinct permutations for each ii
%if n=even and ii=n/2 we have number/2 distinct permutations

number=@(ii) factorial(n)/ (factorial(ii)* ( factorial(n-ii) )  );


if (-1)^n==1 %n is even

temp2=byparts;


                    for ii=1:floor(n/2)


                      if ii==n/2

                        tempnew = temp2{ii};

                        if ii==1
                            col=n;
                        else
                            col = number(ii);

                        end


                        for jj=1:(col/2)

                            for kk=1:col

                             if tempnew(:,jj)==rest{ii,kk}.' 


                               tempnew(:,jj)=[];

                             end



                            end


                        end

                        byparts{ii}=tempnew;
                      else


                          continue

                      end





                    end

end


    %now do the tests;
    %for each ii we need to have either n!/(ii!(n-ii)!) or n!/(ii!(n-ii)!)*1/2
    %distinct bypartitions (for the even in the latter case)



    for ii=1:floor(n/2)

       [~,col]=size(byparts{ii});




       if (-1)^n==1

           if ii==n/2

           ways=factorial(n)/( factorial(ii) * factorial(n-ii) )*1/2;
           else
               ways=factorial(n)/( factorial(ii) * factorial(n-ii) );
           end


       else
           ways=factorial(n)/( factorial(ii) * factorial(n-ii) );
       end

       if col~=ways
          error('Extra bipartitions were found that are repeated.') 

       end

    end




end
