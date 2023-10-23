function [byparts]=all_byparts(n)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Generate all bipartitions for order n. 
%e.g. for n=3 we can have *|123, 1|23, 2|13, 3|12, 12|3, 13|2, 23|1, 123|*.
%(* means empty set)
%Output: The bypartitions in a cell array. They need to be accessed
%        collumnwise from each cell.
%--------------------------------------------------------------------------

b=1:n;

byparts = cell(1,n+1); %but 2^n in total when selected columnwise
for ii=1:n

   byparts{ii}=nchoosek(b,ii).'; 

end
%and add+1 which will be the empty set

byparts{n+1}=[];           
byparts=circshift(byparts,1); %bring the empty set in the beginning



end
