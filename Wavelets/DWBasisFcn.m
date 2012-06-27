function [F, Coeffs] = DWBasisFcn(Tree, Level, Node, Index)
% function F  = DWBasisFcn(Tree, Level, Node, Index)
%
% DWBASISFCN returns a representation of a particular wavelet packet basis
% function or functions with respect to the initial basis.  There are three
% different possible syntices:
%
% DWBasisFcn(Tree, [1 2], [3 4], [1 2]) returns a Mx2 array with the
% function from level 1, node 3, index 1 in the first column and the
% function from level 2, node 4, index 2 in the second.
%
% DWBasisFcn(Tree, 1, 2, [1 2 3]), on the other hand, returns the functions
% with index 1, 2, and 3 from level 1, node 2.
%
% DWBasisFcn(Tree, 1, 2) returns all of the basis functions from level 1,
% node 2.
%
% In:
%    Tree  = the diffusion wavelet packet tree
%    Level = a scalar or vector giving the level or levels of the basis
%            functions to extract
%    Node  = a scalar or vector giving the index of the node or nodes
%    Index = indices of the basis functions
%
% Out:
%    F     = an MxN array specifying N basis functions in R^M
%
% Version History:
%   jcb        2/2006         adapted from old version and modified to
%                             allow for multiple functions

F = [];

% initialize coefficient cell array
 Coeffs = cell(size(Tree));

if ~exist('Index')
   if length(Level) > 1 | length(Node) > 1
      fprintf('DWBasisFcn: invalid syntax\n');
      return;
   end

   Index = 1:size(Tree{Level(1),Node(1)}.Basis,2);
end

% two different possible syntaxes
if length(Index) >= 1 && length(Level) == 1 && length(Node)==1

   % many indices, all in the same node
   M = size(Tree{Level, Node}.Basis, 2);
   N = length(Index);

   Coeffs{Level, Node} = sparse(Index, 1:N, ones(1,N), M, N, N);

   % reconstruct the basis function
   F = DWRecon(Tree, Coeffs);

elseif length(Index) > 1 && length(Level) > 1 && length(Node) > 1

   N = length(Index);
   R = size(Tree{1,1}.Basis,1);

   F = zeros(R, N);

   % many different indices all over the place ... we have to do this
   % function by function
   for j=1:N
      Coeffs = cell(size(Tree));
      M = size(Tree{Level(j), Node(j)}.Basis, 2);
      Coeffs{Level(j), Node(j)} = sparse(Index(j), 1, 1.0, M, 1, 1);

      F(:,j) = DWRecon(Tree, Coeffs);
   end

else

   fprintf('DWBasisFcn.m: unrecognized syntax\n');
end

