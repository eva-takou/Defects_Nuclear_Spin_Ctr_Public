function M = superkron(varargin)
%Jean-Daniel (2023). Superkron (https://www.mathworks.com/matlabcentral/fileexchange/32578-superkron), MATLAB Central File Exchange. Retrieved October 23, 2023.
% M = superkron(M1, M2, [M3, ...]);
% or
% M = superkron({M1, M2, ...})
%
% Computes the kronecker product between several arbitrarily-dimensional
% matrices.
% Written by : Jean-Daniel Bancal
% Last modified : 4.11.2011
% We check the input
if nargin < 2
    if iscell(varargin{1}) % We also accept input in the form superkron({M1,M2,...})
        M = varargin{1}{1};
        for i = 2:length(varargin{1})
            M = superkron(M, varargin{1}{i});
        end
        return;
    else
        M = varargin{1};
        return
    end
end
% Computes the kronecker product between the two first arguments
% First, we care for the dimensions of the matrices
ss1 = size(varargin{1});
s2 = size(varargin{2});
nbDim = max(length(ss1), length(s2));
if length(s2) < nbDim
    s2 = [s2 ones(1,length(ss1)-length(s2))]; % This line includes a suggestion by Mauro Faccin.
elseif length(ss1) < nbDim
    ss1 = [ss1 ones(1,length(s2)-length(ss1))]; % This line includes a suggestion by Mauro Faccin.
end
% We do the product
M = reshape(varargin{1},numel(varargin{1}),1) * reshape(varargin{2},1,numel(varargin{2}));
% And put the numbers back into their place
M = reshape(M,[ss1 s2]);
M = permute(M,reshape([nbDim+1:2*nbDim;1:nbDim],1,2*nbDim));
M = reshape(M, ss1.*s2);
%If the input contains more than two matrices, we iterate the process
for i = 3:nargin
    M = SubClass_Ent_and_Fid.superkron(M, varargin{i});
end


%obj.MM = M;
%dbstop in superkron
end
