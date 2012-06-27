function CoeffTree = DWCoeffs( Tree, Fcns)
% function CoeffTree = DWCoeffs( Tree, Fcns)
%
% Compute the coefficients of the given functions in each of the subspaces
% represented in the given diffusion wavelet/wavelet packet tree.
%
% In:
%    Tree      = a diffusion wavelet tree returned by dwtree or dwptree
%    Fcns      = an mxn matrix of n functions in R^n, each represented in the delta
%                basis
% Out:
%    CoeffTree = a cell array with the same structure as Tree containing the
%                coefficient for each subspace
%
% Dependencies:
%     none
%
% Version History:
%     jcb      2/2006         created from old ComputeCoeffs.m
%
% (c) Copyright Yale University, James C Bremer Jr. and Mauro Maggioni
%

Levels = size(Tree, 1);    % number of levels in the tree
%MaxNodes  = 2^floor(Levels/2); % maximum number of nodes at each level
MaxNodes = size(Tree,2);

% initialize the coefficient tree to have the same structure as Tree
CoeffTree = cell(size(Tree));

% compute the coefficients at the first level
for r=1:MaxNodes
   if ~isempty(Tree{1,r})
      ljout(sprintf('%s', DWNodeName(1, r)), 10);
      
      if isempty(Tree{1,r}.Basis),
         CoeffTree{1,r} = 0;
      else
         CoeffTree{1,r} = Tree{1,r}.Basis'*Fcns;
      end
      fprintf('\b\b\b\b\b\b\b\b\b\b');
   end
end


for j=2:Levels
   for r=1:MaxNodes
      if ~isempty(Tree{j,r})
         ljout(sprintf('%s', DWNodeName(j, r)), 10);
         
         if isempty(Tree{j,r}.Basis),
            CoeffTree{j,r} = 0;
         else
            % represent Fcns in the new basis
            ParentNode = ceil(r/2);
            CoeffTree{j, r} = Tree{j, r}.Basis'*CoeffTree{j-1,ParentNode};
         end
         
         fprintf('\b\b\b\b\b\b\b\b\b\b');
      end
   end

end


