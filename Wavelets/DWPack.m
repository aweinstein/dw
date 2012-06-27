function Basis = DWPack(Tree, CoeffList)

% function Basis2 = DWPack(Tree, CoeffList)
%
%
%

[Levels Nodes] = size(Tree);
Basis = cell(size(Tree));

for Level=1:size(Tree,1)
   for Node=1:size(Tree,2)

       idxs = find(CoeffList(:,1)==Level & CoeffList(:,2)==Node);

       if ~isempty(idxs)
         k = length(idxs);


         Basis{Level, Node} = sparse(CoeffList(idxs, 3), ones(k,1), CoeffList(idxs, 4), size(Tree{Level,Node}.Basis,2), 1, k);
      end


   end
end


% now fill in the coefficients
%for j=1:size(CoeffList,1)
%   Level = CoeffList(j, 1);
%   Node  = CoeffList(j, 2);
%   Index = CoeffList(j, 3);
%   Coeff = CoeffList(j, 4);

%   Basis2{Level, Node}(Index, 1) = Coeff;
%end




