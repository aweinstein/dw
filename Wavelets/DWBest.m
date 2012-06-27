function [Basis, Costs, Marks] = DWBest(Coeffs, CostFcn)
% function Basis = DWBest(Coeffs, [CostFcn])
%
% Find the representation of the function which minimizes the given additive
% cost function.
%
% In:
%    Coeffs  = cell array containing wavelet packet coefficients, returned by
%              DWCoeffs.m
%    CostFcn = an optional argument specifying a cost fcn.
%
% Out:
%    Basis   = a cell array containing the
%

Basis = [];

if nargin < 2
   fprintf('dwbest.m: using default L^1 cost function\n');
   CostFcn = @(x) sum(abs(x));
end

Costs  = zeros(size(Coeffs));
Marks  = zeros(size(Coeffs));

Levels   = size(Coeffs, 1);
MaxNodes = size(Coeffs, 2);

%{
% compute all the costs and mark all nodes with no children
%}
for j=1:Levels
   for r=1:MaxNodes
      if ~isempty(Coeffs{j,r})
         Costs(j,r) = sum(feval(CostFcn, Coeffs{j,r}));
         if ~HasChildren(Coeffs, j, r)
            Marks(j,r) = 1;
         end
      end
   end
end

%{
% mark the nodes that will be used in the best basis
%}
for Level=Levels:-1:1
   for Node=1:MaxNodes

      if ~isempty(Coeffs{Level, Node}) & HasChildren(Coeffs, Level, Node)
         LeftNode  = 2*Node-1;
         RightNode = 2*Node;
         ChildrenCost = Costs(Level+1, LeftNode) + Costs(Level+1, RightNode);

         if ChildrenCost >= Costs(Level, Node)
            Marks = MarkSubtree(Marks, Level, Node, 0);
            Marks(Level, Node) = 1;
         else
            % the children cost less ... we don't need to change any marks, but
            % we need to update the cost
            Costs(Level, Node) = ChildrenCost;
         end
      end
   end
end

%{
% now setup the basis cell array
%}
Basis = cell(size(Coeffs));
for Level=1:Levels
   for Node=1:MaxNodes
      if Marks(Level, Node)
         Basis{Level, Node} = Coeffs{Level, Node};
      end
   end
end

fprintf('DWBest.m: cost of best basis = %g\n', Costs(1,1)+Costs(1,2));

function Children = HasChildren(Coeffs, Level, Node)

Children = 0;

Level = Level+1;
if size(Coeffs,1) >= Level
   LeftNode  = 2*Node-1;
   RightNode = 2*Node;

%   fprintf('(%d, %d) (%d, %d)\n', Level, LeftNode, Level, RightNode);

   if size(Coeffs,2) >= RightNode
      if ~isempty(Coeffs{Level, RightNode}) & ~isempty(Coeffs{Level, LeftNode})
         Children = 1;
      end
   end

end


function Marks = MarkSubtree(Marks, Level, Node, Value)
% function Marks = MarkSubtree
%
% Set the mark value for all nodes lying in the subtree with parent (Level, Node)

if Level <= size(Marks,1) && Node <= size(Marks,2)
   % first set the value for the children
   Marks = MarkSubtree(Marks, Level+1, 2*Node-1, Value);
   Marks = MarkSubtree(Marks, Level+1, 2*Node, Value);

   % now mark this node
   Marks(Level, Node) = Value;
end